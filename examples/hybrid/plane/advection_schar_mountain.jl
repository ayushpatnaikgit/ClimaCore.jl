using Test
using LinearAlgebra, StaticArrays

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
import ClimaCore.Geometry: ⊗

import Logging
import TerminalLoggers
Logging.global_logger(TerminalLoggers.TerminalLogger())

# set up function space
function hvspace_2D(
    xlim = (-π, π),
    zlim = (0, 4π),
    helem = 64,
    velem = 32,
    npoly = 4,
)
    FT = Float64
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_names = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, nelems = velem)
    vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = Domains.IntervalDomain(
        Geometry.XPoint{FT}(xlim[1]),
        Geometry.XPoint{FT}(xlim[2]),
        periodic = true,
    )
    horzmesh = Meshes.IntervalMesh(horzdomain; nelems = helem)
    horztopology = Topologies.IntervalTopology(horzmesh)

    quad = Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    hv_center_space =
        Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)
    return (hv_center_space, hv_face_space)
end

# set up rhs!
hv_center_space, hv_face_space = hvspace_2D((-30000, 30000), (0, 25000))

const MSLP = 1e5 # mean sea level pressure
const grav = 9.8 # gravitational constant
const R_d = 287.058 # R dry (gas constant / mol mass dry air)
const γ_ratio = 1.4 # heat capacity ratio
const C_p = R_d * γ_ratio / (γ_ratio - 1) # heat capacity at constant pressure
const C_v = R_d / (γ_ratio - 1) # heat capacity at constant volume
const R_m = R_d # moist R, assumed to be dry
const uᵣ = 10.0

function pressure(ρθ)
    if ρθ >= 0
        return MSLP * (R_d * ρθ / MSLP)^γ_ratio
    else
        return NaN
    end
end
function rayleigh_sponge(z;
                         z_sponge=15000.0,
                         z_max=25000.0,
                         α = 0.5,  # Relaxation timescale
                         τ = 0.5,
                         γ = 2.0)
    if z >= z_sponge
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α * sinpi(τ * r)^γ
        return β_sponge
    else
        return eltype(z)(0)
    end
end

function rayleigh_sponge_x(x;
                         Δ_sponge=20000.0,
                         α = 0.5,  # Relaxation timescale
                         τ = 0.5,
                         γ = 2.0)
    if x >= 20000.0 ||  x <= -20000.0
        r = (abs(x)-Δ_sponge)  / (10000.0)
        β_sponge = α * sinpi(τ * r)^γ
        return β_sponge
    else
        return eltype(x)(0)
    end
end

Φ(z) = grav * z

# Reference: https://journals.ametsoc.org/view/journals/mwre/140/4/mwr-d-10-05073.1.xml, Section 5a
function init_advection_schar_2d(x, z)
    θ₀ = 280.0
    cp_d = C_p
    cv_d = C_v
    p_0 = MSLP
    𝒩 = 0.01
    θ = θ₀ * exp(𝒩^2 * z / grav)
    π_exn = 1.0 + grav^2/cp_d/θ₀/𝒩^2 * (exp(-𝒩^2*z/grav))
    T = π_exn * θ # temperature
    p = p_0 * π_exn^(cp_d / R_d) # pressure
    ρ = p / R_d / T # density
    ρθ = ρ * θ # potential temperature density
    return (ρ = ρ, ρθ = ρθ)
end

# initial conditions
coords = Fields.coordinate_field(hv_center_space)
face_coords = Fields.coordinate_field(hv_face_space)

Yc = map(coords) do coord
    advection_over_mountain = init_advection_schar_2d(coord.x, coord.z)
    advection_over_mountain
end

ρw = map(face_coords) do coord
    Geometry.WVector(0.0)
end;

Y = Fields.FieldVector(
    Yc = Yc,
    ρuₕ = Yc.ρ .* Ref(Geometry.UVector(uᵣ)),
    ρw = ρw,
)

function energy(Yc, ρu, z)
    ρ = Yc.ρ
    ρθ = Yc.ρθ
    u = ρu / ρ
    kinetic = ρ * norm(u)^2 / 2
    potential = z * grav * ρ
    internal = C_v * pressure(ρθ) / R_d
    return kinetic + potential + internal
end
function combine_momentum(ρuₕ, ρw)
    Geometry.transform(Geometry.UWAxis(), ρuₕ) +
    Geometry.transform(Geometry.UWAxis(), ρw)
end
function center_momentum(Y)
    If2c = Operators.InterpolateF2C()
    combine_momentum.(Y.ρuₕ, If2c.(Y.ρw))
end
function total_energy(Y)
    ρ = Y.Yc.ρ
    ρu = center_momentum(Y)
    ρθ = Y.Yc.ρθ
    z = Fields.coordinate_field(axes(ρ)).z
    sum(energy.(Yc, ρu, z))
end

θ_0 = sum(Yc.ρθ)
mass_0 = sum(Yc.ρ) # Computes ∫ρ∂Ω such that quadrature weighting is accounted for.

function rhs!(dY, Y, _, t)
    ρuₕ = Y.ρuₕ
    ρw = Y.ρw
    Yc = Y.Yc
    dYc = dY.Yc
    dρuₕ = dY.ρuₕ
    dρw = dY.ρw
    ρ = Yc.ρ
    ρθ = Yc.ρθ
    dρθ = dYc.ρθ
    dρ = dYc.ρ

    # spectral horizontal operators
    hdiv = Operators.Divergence()
    hgrad = Operators.Gradient()
    hwdiv = Operators.WeakDivergence()
    hwgrad = Operators.WeakGradient()

    # vertical FD operators with BC's
    vdivf2c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.WVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0)),
    )
    vvdivc2f = Operators.DivergenceC2F(
        bottom = Operators.SetDivergence(Geometry.WVector(0.0)),
        top = Operators.SetDivergence(Geometry.WVector(0.0)),
    )
    uvdivf2c = Operators.DivergenceF2C(
        bottom = Operators.SetValue(
            Geometry.WVector(0.0) ⊗ Geometry.UVector(0.0),
        ),
        top = Operators.SetValue(Geometry.WVector(0.0) ⊗ Geometry.UVector(0.0)),
    )
    If = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    Ic = Operators.InterpolateF2C()
    ∂ = Operators.DivergenceF2C(
        bottom = Operators.SetValue(Geometry.WVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0)),
    )
    ∂f = Operators.GradientC2F()
    ∂c = Operators.GradientF2C()
    B = Operators.SetBoundaryOperator(
        bottom = Operators.SetValue(Geometry.WVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0)),
    )

    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    fcf = Operators.FluxCorrectionF2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )

    uₕ = @. ρuₕ / ρ
    w = @. ρw / If(ρ)
    wc = @. Ic(ρw) / ρ
    p = @. pressure(ρθ)
    θ = @. ρθ / ρ
    Yfρ = @. If(ρ)

    ### HYPERVISCOSITY
    # 1) compute hyperviscosity coefficients
    @. dρθ = hwdiv(hgrad(θ))
    @. dρuₕ = hwdiv(hgrad(uₕ))
    @. dρw = hwdiv(hgrad(w))
    Spaces.weighted_dss!(dYc)
    Spaces.weighted_dss!(dρuₕ)
    Spaces.weighted_dss!(dρw)

    κ₄ = 0.0 # m^4/s
    @. dρθ = -κ₄ * hwdiv(ρ * hgrad(dρθ))
    @. dρuₕ = -κ₄ * hwdiv(ρ * hgrad(dρuₕ))
    @. dρw = -κ₄ * hwdiv(Yfρ * hgrad(dρw))

    # density
    @. dρ = -∂(ρw)
    @. dρ -= hdiv(ρuₕ)

    # potential temperature
    @. dρθ += -(∂(ρw * If(ρθ / ρ)))
    @. dρθ -= hdiv(uₕ * ρθ)

    # horizontal momentum
    Ih = Ref(
        Geometry.Axis2Tensor(
            (Geometry.UAxis(), Geometry.UAxis()),
            @SMatrix [1.0]
        ),
    )
    @. dρuₕ += -uvdivf2c(ρw ⊗ If(uₕ))
    @. dρuₕ -= hdiv(ρuₕ ⊗ uₕ + p * Ih)

    # vertical momentum
    z = coords.z
    @. dρw += B(
        Geometry.transform(Geometry.WAxis(), -(∂f(p)) - If(ρ) * ∂f(Φ(z))) -
        vvdivc2f(Ic(ρw ⊗ w)),
    )
    uₕf = @. If(ρuₕ / ρ) # requires boundary conditions
    @. dρw -= hdiv(uₕf ⊗ ρw)

    ### DIFFUSION
    κ₂ = 75.0 # m^2/s
    #  1a) horizontal div of horizontal grad of horiz momentun
    @. dρuₕ += hwdiv(κ₂ * (ρ * hgrad(ρuₕ / ρ)))
    #  1b) vertical div of vertical grad of horiz momentun
    @. dρuₕ += uvdivf2c(κ₂ * (Yfρ * ∂f(ρuₕ / ρ)))

    #  1c) horizontal div of horizontal grad of vert momentum
    @. dρw += hwdiv(κ₂ * (Yfρ * hgrad(ρw / Yfρ)))
    #  1d) vertical div of vertical grad of vert momentun
    @. dρw += vvdivc2f(κ₂ * (ρ * ∂c(ρw / Yfρ)))

    #  2a) horizontal div of horizontal grad of potential temperature
    @. dρθ += hwdiv(κ₂ * (ρ * hgrad(ρθ / ρ)))
    #  2b) vertical div of vertial grad of potential temperature
    @. dρθ += ∂(κ₂ * (Yfρ * ∂f(ρθ / ρ)))

    # Application of Sponge [Lateral + Top-Boundary Sponge]
    @. dρuₕ -= (rayleigh_sponge_x(coords.x) * (ρuₕ - Geometry.UVector(uᵣ)) + rayleigh_sponge(coords.z) * (ρuₕ - Geometry.UVector(uᵣ)))
    @. dρw -= (rayleigh_sponge_x(face_coords.x) * ρw + rayleigh_sponge(face_coords.z) * ρw)

    Spaces.weighted_dss!(dYc)
    Spaces.weighted_dss!(dρuₕ)
    Spaces.weighted_dss!(dρw)
    return dY
end

dYdt = similar(Y);
rhs!(dYdt, Y, nothing, 0.0);


# run!
using OrdinaryDiffEq
Δt = 0.1
prob = ODEProblem(rhs!, Y, (0.0, 900.0))

integrator = OrdinaryDiffEq.init(
    prob,
    SSPRK33(),
    dt = Δt,
    saveat = 50.0,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
);

if haskey(ENV, "CI_PERF_SKIP_RUN") # for performance analysis
    throw(:exit_profile)
end

sol = @timev OrdinaryDiffEq.solve!(integrator)

ENV["GKSwstype"] = "nul"
using ClimaCorePlots, Plots
Plots.GRBackend()

dir = "schar_mountain"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

# post-processing
using ClimaCorePlots, Plots
anim = Plots.@animate for u in sol.u
    Plots.plot(u.Yc.ρθ ./ u.Yc.ρ)
end
Plots.mp4(anim, joinpath(path, "theta.mp4"), fps = 20)

If2c = Operators.InterpolateF2C()
anim = Plots.@animate for u in sol.u
    Plots.plot(If2c.(u.ρw) ./ u.Yc.ρ)
end
Plots.mp4(anim, joinpath(path, "vel_w.mp4"), fps = 20)

anim = Plots.@animate for u in sol.u
    Plots.plot(u.ρuₕ ./ u.Yc.ρ)
end
Plots.mp4(anim, joinpath(path, "vel_u.mp4"), fps = 20)

θs = [sum(u.Yc.ρθ) for u in sol.u]
Mass = [sum(u.Yc.ρ) for u in sol.u]

Plots.png(Plots.plot((θs .- θ_0) ./ θ_0), joinpath(path, "energy.png"))
Plots.png(Plots.plot((Mass .- mass_0) ./ mass_0), joinpath(path, "mass.png"))
