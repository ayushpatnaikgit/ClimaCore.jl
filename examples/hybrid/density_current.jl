push!(LOAD_PATH, joinpath(@__DIR__, "..", ".."))

using Test
using StaticArrays, IntervalSets, LinearAlgebra, UnPack

import ClimaCore:
    ClimaCore,
    slab,
    Spaces,
    Domains,
    Meshes,
    Geometry,
    Topologies,
    Spaces,
    Fields,
    Operators
import ClimaCore.Domains.Geometry: Cartesian2DPoint
using ClimaCore.Geometry
using JLD2

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

# set up function space
function hvspace_2D(
    xlim = (-π, π),
    zlim = (0, 4π),
    helem = 10,
    velem = 50,
    npoly = 4,
)
    FT = Float64
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_tags = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, nelems = velem)
    vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)
    horzdomain = Domains.RectangleDomain(
        Geometry.XPoint{FT}(xlim[1])..Geometry.XPoint{FT}(xlim[2]),
        Geometry.YPoint{FT}(-0)..Geometry.YPoint{FT}(0),
        x1periodic = true,
        x2boundary = (:a, :b),
    )
    horzmesh = Meshes.EquispacedRectangleMesh(horzdomain, helem, 1)
    horztopology = Topologies.GridTopology(horzmesh)

    quad = Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    hv_center_space =
        Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)
    return (hv_center_space, hv_face_space)
end

# set up rhs!
#hv_center_space, hv_face_space = hvspace_2D((-25600, 25600), (0, 6400), 128, 64,4)
hv_center_space, hv_face_space = hvspace_2D((-25600, 25600), (0, 6400), 64, 32,4)
# Δx = 100m
# Δz = 100m

const MSLP = 1e5 # mean sea level pressure
const grav = 9.8 # gravitational constant
const R_d = 287.058 # R dry (gas constant / mol mass dry air)
const γ = 1.4 # heat capacity ratio
const C_p = R_d * γ / (γ - 1) # heat capacity at constant pressure
const C_v = R_d / (γ - 1) # heat capacity at constant volume
const R_m = R_d # moist R, assumed to be dry

function pressure(ρθ)
    if ρθ >= 0
        return MSLP * (R_d * ρθ / MSLP)^γ
    else
        return NaN
    end
end

Φ(z) = grav * z

# Reference: https://journals.ametsoc.org/view/journals/mwre/140/4/mwr-d-10-05073.1.xml, Section 5a
function init_dry_rising_current_2d(x, z)
    x_c = 0.0
    z_c = 3000.0
    θ_b = 300.0
    θ_c = -16.632 # Ullrich & Jablonowski apply the perturbation to T rather than to θ, hence this magnitude
    cp_d = C_p
    cv_d = C_v
    p_0 = MSLP
    g = grav
    x_r = 4000
    z_r = 2000

    # auxiliary quantities
    r = sqrt((x - x_c)^2/x_r^2 + (z - z_c)^2/z_r^2)
    θ_p = r <= 1.0 ? 0.5 * θ_c * (1.0 + cospi(r)) : 0.0 # potential temperature perturbation

    θ = θ_b # potential temperature
    π_exn = 1.0 - g * z / cp_d / θ # exner function
    T = π_exn * θ # temperature
    p = p_0 * π_exn^(cp_d / R_d) # pressure
    ρ = p / R_d / T # density
    ρθ = ρ * (θ_b + θ_p) # potential temperature density

    return (ρ = ρ, ρθ = ρθ, ρuₕ = ρ * Geometry.Cartesian1Vector(0.0))
end

# initial conditions
coords = Fields.coordinate_field(hv_center_space);
face_coords = Fields.coordinate_field(hv_face_space);

Yc = map(coords) do coord
    current = init_dry_rising_current_2d(coord.x, coord.z)
    current
end;

ρw = map(face_coords) do coord
    Geometry.Cartesian3Vector(0.0)
end;

Y = Fields.FieldVector(Yc = Yc, ρw = ρw);

function rhs!(dY, Y, _, t)
    ρw = Y.ρw
    Yc = Y.Yc
    dYc = dY.Yc
    dρw = dY.ρw

    # spectral horizontal operators
    hdiv = Operators.Divergence()
    hgrad = Operators.Gradient()
    hwdiv = Operators.WeakDivergence()
    hwgrad = Operators.WeakGradient()

    # vertical FD operators with BC's
    vdivf2c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
        top = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
    )
    vvdivc2f = Operators.DivergenceC2F(
        bottom = Operators.SetDivergence(Geometry.Cartesian3Vector(0.0)),
        top = Operators.SetDivergence(Geometry.Cartesian3Vector(0.0)),
    )
    uvdivf2c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(
            Geometry.Cartesian3Vector(0.0) ⊗ Geometry.Cartesian1Vector(0.0),
        ),
        top = Operators.SetValue(
            Geometry.Cartesian3Vector(0.0) ⊗ Geometry.Cartesian1Vector(0.0),
        ),
    )
    If_bc = Operators.InterpolateC2F(
        bottom = Operators.SetValue(Geometry.Cartesian1Vector(0.0)),
        top = Operators.SetValue(Geometry.Cartesian1Vector(0.0)),
    )
    If = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    Ic = Operators.InterpolateF2C()
    ∂ = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
        top = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
    )
    ∂f = Operators.GradientC2F()
    ∂c = Operators.GradientF2C()
    B = Operators.SetBoundaryOperator(
        bottom = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
        top = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
    )

    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    fcf = Operators.FluxCorrectionF2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )



    # 1) compute hyperviscosity coefficients
    #
    @. dYc.ρ = hdiv(hgrad(Yc.ρ))
    @. dYc.ρθ = hdiv(hgrad(Yc.ρθ))
    @. dYc.ρuₕ = hdiv(hgrad(Yc.ρuₕ))
    @. dρw = hdiv(hgrad(ρw))
    Spaces.weighted_dss!(dYc)

    κ = 0.0
    @. dYc.ρ = κ * hdiv(hgrad(dYc.ρ))
    @. dYc.ρθ = κ * hdiv(hgrad(dYc.ρθ))
    @. dYc.ρuₕ = κ * hdiv(hgrad(dYc.ρuₕ))
    @. dρw = κ * hdiv(hgrad(dρw))

    uₕ = @. Yc.ρuₕ / Yc.ρ
    w = @. ρw / If(Yc.ρ)
    wc = @. Ic(ρw) / Yc.ρ
    p = @. pressure(Yc.ρθ)

    # density
    @. dYc.ρ += -∂(ρw) + fcc(w, Yc.ρ)
    @. dYc.ρ -= hdiv(Yc.ρuₕ)

    # potential temperature
    @. dYc.ρθ += -(∂(ρw * If(Yc.ρθ / Yc.ρ))) + fcc(w, Yc.ρθ)
    @. dYc.ρθ -= hdiv(uₕ * Yc.ρθ)

    # horizontal momentum
    Ih = Ref(
        Geometry.Axis2Tensor(
            (Geometry.Cartesian1Axis(), Geometry.Cartesian1Axis()),
            @SMatrix [1.0]
        ),
    )
    @. dYc.ρuₕ += -uvdivf2c(ρw ⊗ If_bc(uₕ)) + fcc(w, Yc.ρuₕ)
    @. dYc.ρuₕ -= hdiv(Yc.ρuₕ ⊗ uₕ + p * Ih)


    # vertical momentum
    @. dρw +=
        B(
            Geometry.transform(
                Geometry.Cartesian3Axis(),
                -(∂f(p)) - If(Yc.ρ) * ∂f(Φ(coords.z)),
            ) - vvdivc2f(Ic(ρw ⊗ w)),
        ) + fcf(wc, ρw)
    uₕf = @. If_bc(Yc.ρuₕ / Yc.ρ) # requires boundary conditions
    @. dρw -= hdiv(uₕf ⊗ ρw)

    # 3) diffusion
    

    κ = 75.0 # m^2/s

    Yfρ = @. If(Yc.ρ) # Face interpolation

 #  hgrad defined on spectral element nodes
 #  ∂ = vgrad ... = defined on centers
    ∇ₕuₕ = @. hgrad(Yc.ρuₕ / Yc.ρ)
    ∇ᵥuₕ = @. ∂f(Yc.ρuₕ / Yc.ρ)
    ∇ₕw = @. hgrad(ρw / Yfρ)    
    ∇ᵥw = @. ∂c(ρw / Yfρ)    

    #  1a) horizontal div of horizontal grad of horiz momentun
    ρκc = @. Yc.ρ * κ 
    ρκf = @. Yfρ * κ
    @. dYc.ρuₕ += hdiv(ρκc * (hgrad(Yc.ρuₕ / Yc.ρ)))

    #  1b) vertical div of vertical grad of horiz momentun
    @. dYc.ρuₕ += uvdivf2c(ρκf * (∂f(Yc.ρuₕ / Yc.ρ)))

    #  1c) horizontal div of horizontal grad of vert momentum
    @. dρw += hdiv(ρκf * (hgrad(ρw / Yfρ)))

    #  1d) vertical div of vertical grad of vert momentun
    @. dρw += vvdivc2f(ρκc *  (∂c(ρw / Yfρ)))

    # 2a) horizontal div of horizontal grad of potential temperature
    @. dYc.ρθ += hdiv(ρκc * (hgrad(Yc.ρθ / Yc.ρ)))

    # 2b) vertical div of vertial grad of potential temperature
    @. dYc.ρθ += ∂(ρκf * (∂f(Yc.ρθ / Yc.ρ)))

    Spaces.weighted_dss!(dYc)
    Spaces.weighted_dss!(dρw)
    return dY
end

dYdt = similar(Y);
rhs!(dYdt, Y, nothing, 0.0)


# run!
using OrdinaryDiffEq
Δt = 0.2
prob = ODEProblem(rhs!, Y, (0.0, 900.0))
sol = solve(
    prob,
    SSPRK33(),
    dt = Δt,
    saveat = 50.0,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);
@save "current_output.jld2" sol

ENV["GKSwstype"] = "nul"
import Plots
Plots.GRBackend()

dirname = "current"
path = joinpath(@__DIR__, "output", dirname)
mkpath(path)

# post-processing
import Plots

anim = Plots.@animate for u in sol.u
    Plots.plot(u.Yc.ρθ ./ u.Yc.ρ)
    #Plots.plot(Operators.matrix_interpolate(u.Yc.ρθ ./ u.Yc.ρ),4)
end
Plots.mp4(anim, joinpath(path, "current_theta.mp4"), fps = 20)

If = Operators.InterpolateC2F(
   bottom = Operators.Extrapolate(),
   top = Operators.Extrapolate(),
)

anim = Plots.@animate for u in sol.u
    Plots.plot(@. u.ρw ./ If(u.Yc.ρ))
end
Plots.mp4(anim, joinpath(path, "current_w.mp4"), fps = 20)

anim = Plots.@animate for u in sol.u
    Plots.plot(Operators.matrix_interpolate(u.Yc.ρθ ./ u.Yc.ρ),4)
end
Plots.mp4(anim, joinpath(path, "current_theta_itp.mp4"), fps = 20)
