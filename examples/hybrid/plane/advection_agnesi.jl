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
    Hypsography,
    Operators
using ClimaCore.Geometry

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using DiffEqCallbacks

const MSLP = 1e5 # mean sea level pressure
const grav = 9.8 # gravitational constant
const R_d = 287.058 # R dry (gas constant / mol mass dry air)
const γ = 1.4 # heat capacity ratio
const C_p = R_d * γ / (γ - 1) # heat capacity at constant pressure
const C_v = R_d / (γ - 1) # heat capacity at constant volume
const R_m = R_d # moist R, assumed to be dry
const FT = Float64
const uᵣ = FT(0)

function warp_agnesi(
    coord;
    h₀ = 250,
    λ = 4000,
    a_c = 5000,
)       
    x = coord.x
    return h₀*exp(-(x/a_c)^2)*(cos(π*x/λ))^2
end

function hvspace_2D(
    xlim = (0, π),
    zlim = (0, 1),
    helem = 2,
    velem = 10,
    npoly = 5;
    stretch = Meshes.Uniform(),
    warp_fn = warp_agnesi,
)
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_tags = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, stretch, nelems = velem)
    vert_face_space = Spaces.FaceFiniteDifferenceSpace(vertmesh)

    # Generate Horizontal Space
    horzdomain = Domains.IntervalDomain(
        Geometry.XPoint{FT}(xlim[1]),
        Geometry.XPoint{FT}(xlim[2]);
        periodic = true,
    )
    horzmesh = Meshes.IntervalMesh(horzdomain; nelems = helem)
    horztopology = Topologies.IntervalTopology(horzmesh)
    quad = Spaces.Quadratures.GLL{npoly + 1}()
    hspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    # Extrusion
    z_surface = warp_fn.(Fields.coordinate_field(hspace))
    f_space = Spaces.ExtrudedFiniteDifferenceSpace(
        hspace,
        vert_face_space,
        Hypsography.LinearAdaption(),
        z_surface,
    )
    c_space = Spaces.CenterExtrudedFiniteDifferenceSpace(f_space)

    return (c_space, f_space)
end

function hvspace_2D_nowarp(
    xlim = (0, π),
    zlim = (0, 1),
    helem = 2,
    velem = 10,
    npoly = 5;
    stretch = Meshes.Uniform(),
)
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_tags = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, stretch, nelems = velem)
    vert_face_space = Spaces.FaceFiniteDifferenceSpace(vertmesh)

    # Generate Horizontal Space
    horzdomain = Domains.IntervalDomain(
        Geometry.XPoint{FT}(xlim[1]),
        Geometry.XPoint{FT}(xlim[2]);
        periodic = true,
    )
    horzmesh = Meshes.IntervalMesh(horzdomain; nelems = helem)
    horztopology = Topologies.IntervalTopology(horzmesh)
    quad = Spaces.Quadratures.GLL{npoly + 1}()
    hspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    # Extrusion
    f_space = Spaces.ExtrudedFiniteDifferenceSpace(
        hspace,
        vert_face_space,
    )
    c_space = Spaces.CenterExtrudedFiniteDifferenceSpace(f_space)

    return (c_space, f_space)
end

using ClimaCorePlots, Plots


# set up function space
#(hv_center_space, hv_face_space) = hvspace_2D((-30000, 30000), (0, 25000), 30, 40, 4;
#                                            stretch = Meshes.Uniform(), warp_fn=warp_agnesi);
(hv_center_space, hv_face_space) = hvspace_2D_nowarp((-30000, 30000), (0, 10000), 30, 40, 4;
                                            stretch = Meshes.Uniform());
function pressure(ρθ)
    if ρθ >= 0
        return MSLP * (R_d * ρθ / MSLP)^γ
    else
        return NaN
    end
end

Φ(z) = grav * z
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

# Reference: https://journals.ametsoc.org/view/journals/mwre/140/4/mwr-d-10-05073.1.xml, Section
function init_agnesi_2d(x, z)
    θ₀ = 280.0
    cp_d = C_p
    cv_d = C_v
    p₀ = MSLP
    g = grav
    γ = cp_d / cv_d
    
    𝒩 = 0.01
    π_exner = @. exp(-g * z / (cp_d * θ₀))
    θ = @. θ₀ * exp(𝒩 ^2 * z / g)
    ρ = @. p₀ / (R_d * θ) * (π_exner)^(cp_d/R_d)
    ρθ  = @. ρ * θ
    ρuₕ = @. ρ * Geometry.UVector(uᵣ)

    return (ρ = ρ,
            ρθ = ρθ,
            ρuₕ = ρuₕ)
end

# initial conditions
coords = Fields.coordinate_field(hv_center_space)
face_coords = Fields.coordinate_field(hv_face_space)

Yc = map(coords) do coord
    agnesi = init_agnesi_2d(coord.x, coord.z)
    agnesi
end

ρw = map(face_coords) do coord
    Geometry.WVector(0.0)
end;

Spaces.weighted_dss!(Yc)
Spaces.weighted_dss!(ρw)

Y = Fields.FieldVector(Yc = Yc, ρw = ρw)

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
    combine_momentum.(Y.Yc.ρuₕ, If2c.(Y.ρw))
end
function total_energy(Y)
    ρ = Y.Yc.ρ
    ρu = center_momentum(Y)
    ρθ = Y.Yc.ρθ
    z = Fields.coordinate_field(axes(ρ)).z
    sum(energy.(Yc, ρu, z))
end

energy_0 = total_energy(Y)
mass_0 = sum(Yc.ρ) # Computes ∫ρ∂Ω such that quadrature weighting is accounted for.

function rhs!(dY, Y, params, t)
    
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
    BW = Operators.SetBoundaryOperator(
        bottom = Operators.SetValue(Geometry.WVector(0.0)),
        top = Operators.SetValue(Geometry.WVector(0.0)),
    )
    BU = Operators.SetBoundaryOperator(
        bottom = Operators.SetValue(Geometry.UVector(0.0)),
        top = Operators.SetValue(Geometry.UVector(0.0)),
    )


    fcc = Operators.FluxCorrectionC2C(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    fcf = Operators.FluxCorrectionF2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )

    uₕ = @. Yc.ρuₕ / Yc.ρ
    w = @. ρw / If(Yc.ρ)
    wc = @. Ic(ρw) / Yc.ρ
    p = @. pressure(Yc.ρθ)
    θ = @. Yc.ρθ / Yc.ρ
    Yfρ = @. If(Yc.ρ)

    ### HYPERVISCOSITY
    # 1) compute hyperviscosity coefficients
    @. dYc.ρθ = hdiv(hgrad(θ))
    @. dYc.ρuₕ = hdiv(hgrad(uₕ))
    @. dρw = hdiv(hgrad(w))
    Spaces.weighted_dss!(dYc)
    Spaces.weighted_dss!(dρw)

    κ₄ = 1e5 # m^4/s
    @. dYc.ρθ = -κ₄ * hdiv(Yc.ρ * hgrad(dYc.ρθ))
    @. dYc.ρuₕ = -κ₄ * hdiv(Yc.ρ * hgrad(dYc.ρuₕ))
    @. dρw = -κ₄ * hdiv(Yfρ * hgrad(dρw))

    # density
    @. dYc.ρ = -∂(ρw)
    @. dYc.ρ -= hdiv(Yc.ρuₕ)

    # potential temperature
    @. dYc.ρθ += -(∂(ρw * If(Yc.ρθ / Yc.ρ)))
    @. dYc.ρθ -= hdiv(uₕ * Yc.ρθ)

    # horizontal momentum
    Ih = Ref(
        Geometry.Axis2Tensor(
            (Geometry.UAxis(), Geometry.UAxis()),
            @SMatrix [1.0]
        ),
    )
    @. dYc.ρuₕ -= uvdivf2c(ρw ⊗ If(uₕ))
    @. dYc.ρuₕ -= hdiv(Yc.ρuₕ ⊗ uₕ + p * Ih)

    # vertical momentum

    # vertical component of vertical momentum
    @. dρw +=
        BW(
            Geometry.project( # project
                Geometry.WAxis(),
                -(∂f(p)) - If(Yc.ρ) * ∂f(Φ(coords.z)),
            ) - vvdivc2f(Ic(ρw ⊗ w)),
        )
    # horizontal component of vertical momentum
    @. dYc.ρuₕ += @. Ic(BU(
            Geometry.project( # project
                Geometry.UAxis(),
                -(∂f(p)) - If(Yc.ρ) * ∂f(Φ(coords.z)),
            ),
        ))

    # vertical component of horizontal momentum
    uₕf = @. If(Yc.ρuₕ / Yc.ρ) # requires boundary conditions
    @. dρw -= hdiv(uₕf ⊗ ρw)

    ### DIFFUSION
    κ₂ = 0.0 # m^2/s
    #  1a) horizontal div of horizontal grad of horiz momentun
    @. dYc.ρuₕ += hwdiv(κ₂ * (Yc.ρ * hgrad(Yc.ρuₕ / Yc.ρ)))
    #  1b) vertical div of vertical grad of horiz momentun
    @. dYc.ρuₕ += uvdivf2c(κ₂ * (Yfρ * ∂f(Yc.ρuₕ / Yc.ρ)))

    #  1c) horizontal div of horizontal grad of vert momentum
    @. dρw += hwdiv(κ₂ * (Yfρ * hgrad(ρw / Yfρ)))
    #  1d) vertical div of vertical grad of vert momentun
    @. dρw += vvdivc2f(κ₂ * (Yc.ρ * ∂c(ρw / Yfρ)))

    #  2a) horizontal div of horizontal grad of potential temperature
    @. dYc.ρθ += hwdiv(κ₂ * (Yc.ρ * hgrad(Yc.ρθ / Yc.ρ)))
    #  2b) vertical div of vertial grad of potential temperature
    @. dYc.ρθ += ∂(κ₂ * (Yfρ * ∂f(Yc.ρθ / Yc.ρ)))

    # sponge
    β = @. rayleigh_sponge(coords.z)
    βx = @. rayleigh_sponge_x(coords.x)
    ρuᵣ = @. Yc.ρ * Geometry.UVector(uᵣ)
    @. dYc.ρuₕ -= β * (Yc.ρuₕ - ρuᵣ)
    @. dYc.ρuₕ -= βx * (Yc.ρuₕ - ρuᵣ)
    @. dρw -= If(β) * ρw
    @. dρw -= If(βx) * ρw

    Spaces.weighted_dss!(dYc)
    Spaces.weighted_dss!(dρw)
    return dY
end

dYdt = similar(Y);

dt = 0.5
energy_0 = total_energy(Y)
mass_0 = sum(Yc.ρ) # Computes ∫ρ∂Ω such that quadrature weighting is accounted for.
cb_dss = PeriodicCallback(
        int -> map(f -> Spaces.weighted_dss!(f), Fields._values(int.u)),
        dt;
        initial_affect = true,
       )
# run!
using OrdinaryDiffEq
timeend = 3600*1
prob = ODEProblem(rhs!, Y, (0.0, timeend))
sol = solve(
    prob,
    SSPRK33(),
    dt = dt,
    saveat = 3600*0.5,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
    callback = cb_dss
);

ENV["GKSwstype"] = "nul"
import Plots
Plots.GRBackend()
dirname = "agnesi_0flow"
path = joinpath(@__DIR__, "output", dirname)
mkpath(path)
# post-processing
If2c = Operators.InterpolateF2C()
u = sol.u
p1 = Plots.plot(If2c.(u[end].ρw) ./ u[end].Yc.ρ, xlim = (-30000, 30000), ylim = (0, 25000))
p2 = Plots.plot(u[end].Yc.ρuₕ ./ u[end].Yc.ρ, xlim = (-30000, 30000), ylim = (0, 25000))
p3 = Plots.plot(u[end].Yc.ρθ ./ u[end].Yc.ρ, xlim = (-30000, 30000), ylim = (0, 25000))
Plots.savefig(p1, path*"/vel_w.png")
Plots.savefig(p2, path*"/vel_u.png")
Plots.savefig(p3, path*"/theta.png")
Es = [total_energy(u) for u in sol.u]
Mass = [sum(u.Yc.ρ) for u in sol.u]
Plots.png(
    Plots.plot((Es .- energy_0) ./ energy_0),
    joinpath(path, "energy.png"),
)
Plots.png(Plots.plot((Mass .- mass_0) ./ mass_0), joinpath(path, "mass.png"))
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
    Plots.plot(u.Yc.ρuₕ ./ u.Yc.ρ)
end
Plots.mp4(anim, joinpath(path, "vel_u.mp4"), fps = 20)
