push!(LOAD_PATH, joinpath(@__DIR__, "..", ".."))

using Test
using StaticArrays, IntervalSets, LinearAlgebra, UnPack
using NCDatasets
using Dierckx

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
    Operators,
    Hypsography
using ClimaCore.Geometry

using DiffEqCallbacks

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

# Warping function as prescribed in GMD-8-317-2015

## Topography Function Generation
include("/Users/asridhar/Research/Codes/Topography/read_topo.jl")

dlonx = range(0,3.8e6,length=128*4)
dlon = range(-80,-40,length=128*4)
elev_dlon = earth_spline
relev = elev_dlon.(collect(dlon),-30)
const disc_elev = Spline1D(dlonx, relev)

function lat2meters(lat1, lon1, lon2)
    lat2 = lat1
    R = 6378.137;
    dLat = lat2 * π / 180 - lat1 * π / 180;
    dLon = lon2 * π / 180 - lon1 * π / 180;
    a = sin(dLat/2) * sin(dLat/2) + cos(lat1 * π / 180) * cos(lat2 * π / 180) * sin(dLon/2) * sin(dLon/2);
    c = 2 * atan(sqrt(a), sqrt(1-a));
    d = R * c;
    return d * 1000
end

function warp_surface(coord)
  x = Geometry.component(coord,1)
  return disc_elev(x)
end
function hvspace_2D(
    xlim = (-π, π),
    zlim = (0, 4π),
    xelem = 128,
    zelem = 64,
    npoly = 4,
    warp_fn = warp_surface,
)
    FT = Float64
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_names = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, nelems = zelem)
    vert_face_space = Spaces.FaceFiniteDifferenceSpace(vertmesh)

    horzdomain = Domains.IntervalDomain(
        Geometry.XPoint{FT}(xlim[1]),
        Geometry.XPoint{FT}(xlim[2]);
        periodic = true,
    )
    horzmesh = Meshes.IntervalMesh(horzdomain, nelems = xelem)
    horztopology = Topologies.IntervalTopology(horzmesh)
    quad = Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = Spaces.SpectralElementSpace1D(horztopology, quad)
  
    z_surface = warp_fn.(Fields.coordinate_field(horzspace))
    hv_face_space = Spaces.ExtrudedFiniteDifferenceSpace(
                    horzspace,
                    vert_face_space,
                    Hypsography.LinearAdaption(), 
                    z_surface
              )
    hv_center_space =
        Spaces.CenterExtrudedFiniteDifferenceSpace(hv_face_space)
    return (hv_center_space, hv_face_space)
end

# set up 2D domain - doubly periodic box
hv_center_space, hv_face_space = hvspace_2D((0, 3.8e6), (0, 10000))

const MSLP = 1e5 # mean sea level pressure
const grav = 9.8 # gravitational constant
const R_d = 287.058 # R dry (gas constant / mol mass dry air)
const γ = 1.4 # heat capacity ratio
const C_p = R_d * γ / (γ - 1) # heat capacity at constant pressure
const C_v = R_d / (γ - 1) # heat capacity at constant volume
const T_0 = 273.16 # triple point temperature

Φ(z) = grav * z

# Reference: https://journals.ametsoc.org/view/journals/mwre/140/4/mwr-d-10-05073.1.xml, Section 5a
# Prognostic thermodynamic variable: Total Energy 
function init_dry_density_current_2d(x, z)
    x_c = 1e6
    z_c = 6500.0
    r_c = 1.0
    x_r = 100000.0
    z_r = 750.0
    θ_b = 300.0
    θ_c = -15.0
    cp_d = C_p
    cv_d = C_v
    p_0 = MSLP
    g = grav

    # auxiliary quantities
    r = sqrt((x - x_c)^2 / x_r^2 + (z - z_c)^2 / z_r^2)
    θ_p = r < r_c ? 0.5 * θ_c * (1.0 + cospi(r / r_c)) : 0.0 # potential temperature perturbation

    θ = θ_b + θ_p # potential temperature
    π_exn = 1.0 - Φ(z) / cp_d / θ # exner function
    T = π_exn * θ # temperature
    p = p_0 * π_exn^(cp_d / R_d) # pressure
    ρ = p / R_d / T # density
    e = cv_d * (T - T_0) + Φ(z)
    ρe = ρ * e # total energy

    return (ρ = ρ, ρe = ρe)
end

# initial conditions
coords = Fields.coordinate_field(hv_center_space)
face_coords = Fields.coordinate_field(hv_face_space)

Yc = map(coord -> init_dry_density_current_2d(coord.x, coord.z), coords)
uₕ = map(_ -> Geometry.Covariant1Vector(0.0), coords)
w = map(_ -> Geometry.Covariant3Vector(0.0), face_coords)
Y = Fields.FieldVector(Yc = Yc, uₕ = uₕ, w = w)

energy_0 = sum(Y.Yc.ρe)
mass_0 = sum(Y.Yc.ρ)

function rhs_invariant!(dY, Y, _, t)

    cρ = Y.Yc.ρ # scalar on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant1Vector on centers
    cρe = Y.Yc.ρe

    dρ = dY.Yc.ρ
    dw = dY.w
    duₕ = dY.uₕ
    dρe = dY.Yc.ρe
    z = coords.z

    # 0) update w at the bottom
    # fw = -g^31 cuₕ/ g^33

    hdiv = Operators.Divergence()
    hwdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    hwgrad = Operators.WeakGradient()
    hcurl = Operators.Curl()

    If2c = Operators.InterpolateF2C()
    Ic2f = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )

    dρ .= 0 .* cρ

    cw = If2c.(fw)
    fuₕ = Ic2f.(cuₕ)
    cuw = Geometry.Covariant13Vector.(cuₕ) .+ Geometry.Covariant13Vector.(cw)

    ce = @. cρe / cρ
    cI = @. ce - Φ(z) - (norm(cuw)^2) / 2
    cT = @. cI / C_v + T_0
    cp = @. cρ * R_d * cT

    h_tot = @. ce + cp / cρ # Total enthalpy at cell centers

    ### HYPERVISCOSITY
    # 1) compute hyperviscosity coefficients
    χe = @. dρe = hwdiv(hgrad(h_tot)) # we store χe in dρe
    χuₕ = @. duₕ = hwgrad(hdiv(cuₕ))

    Spaces.weighted_dss!(dρe)
    Spaces.weighted_dss!(duₕ)

    κ₄ = 0.0 # m^4/s
    @. dρe = -κ₄ * hwdiv(cρ * hgrad(χe))
    @. duₕ = -κ₄ * (hwgrad(hdiv(χuₕ)))

    # 1) Mass conservation
    dw .= fw .* 0

    # 1.a) horizontal divergence
    dρ .-= hdiv.(cρ .* (cuw))

    # 1.b) vertical divergence
    vdivf2c = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    vdivc2f = Operators.DivergenceC2F(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    # we want the total u³ at the boundary to be zero: we can either constrain
    # both to be zero, or allow one to be non-zero and set the other to be its
    # negation

    # explicit part
    dρ .-= vdivf2c.(Ic2f.(cρ .* cuₕ))
    # implicit part
    dρ .-= vdivf2c.(Ic2f.(cρ) .* fw)

    # 2) Momentum equation

    # curl term
    hcurl = Operators.Curl()
    # effectively a homogeneous Dirichlet condition on u₁ at the boundary
    vcurlc2f = Operators.CurlC2F(
        bottom = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
        top = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
    )

    fω¹ = hcurl.(fw)
    fω¹ .+= vcurlc2f.(cuₕ)

    # cross product
    # convert to contravariant
    # these will need to be modified with topography
    fu = Geometry.Contravariant13Vector.(Ic2f.(cuₕ)) .+ Geometry.Contravariant13Vector.(fw)
    fu¹ = Geometry.project.(Ref(Geometry.Contravariant1Axis()), fu)
    fu³ = Geometry.project.(Ref(Geometry.Contravariant3Axis()), fu)
    @. dw -= fω¹ × fu¹ # Covariant3Vector on faces
    @. duₕ -= If2c(fω¹ × fu³)


    @. duₕ -= hgrad(cp) / cρ
    vgradc2f = Operators.GradientC2F(
        bottom = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
        top = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
    )
    @. dw -= vgradc2f(cp) / Ic2f(cρ)

    cE = @. (norm(cuw)^2) / 2 + Φ(z)
    @. duₕ -= hgrad(cE)
    @. dw -= vgradc2f(cE)

    # 3) total energy

    @. dρe -= hdiv(cuw * (cρe + cp))
    @. dρe -= vdivf2c(fw * Ic2f(cρe + cp))
    @. dρe -= vdivf2c(Ic2f(cuₕ * (cρe + cp)))

    # Uniform 2nd order diffusion
    ∂c = Operators.GradientF2C()
    fρ = @. Ic2f(cρ)
    κ₂ = 75.0 # m^2/s

    ᶠ∇ᵥuₕ = @. vgradc2f(cuₕ.components.data.:1)
    ᶜ∇ᵥw = @. ∂c(fw.components.data.:1)
    ᶠ∇ᵥh_tot = @. vgradc2f(h_tot)

    ᶜ∇ₕuₕ = @. hgrad(cuₕ.components.data.:1)
    ᶠ∇ₕw = @. hgrad(fw.components.data.:1)
    ᶜ∇ₕh_tot = @. hgrad(h_tot)

    hκ₂∇²uₕ = @. hwdiv(κ₂ * ᶜ∇ₕuₕ)
    vκ₂∇²uₕ = @. vdivf2c(κ₂ * ᶠ∇ᵥuₕ)
    hκ₂∇²w = @. hwdiv(κ₂ * ᶠ∇ₕw)
    vκ₂∇²w = @. vdivc2f(κ₂ * ᶜ∇ᵥw)
    hκ₂∇²h_tot = @. hwdiv(cρ * κ₂ * ᶜ∇ₕh_tot)
    vκ₂∇²h_tot = @. vdivf2c(fρ * κ₂ * ᶠ∇ᵥh_tot)

    dfw = dY.w.components.data.:1
    dcu = dY.uₕ.components.data.:1

    # Laplacian Diffusion (Uniform)
    @. dcu += hκ₂∇²uₕ
    @. dcu += vκ₂∇²uₕ
    @. dfw += hκ₂∇²w
    @. dfw += vκ₂∇²w
    @. dρe += hκ₂∇²h_tot
    @. dρe += vκ₂∇²h_tot

    Spaces.weighted_dss!(dY.Yc)
    Spaces.weighted_dss!(dY.uₕ)
    Spaces.weighted_dss!(dY.w)

    return dY
end

dYdt = similar(Y);
rhs_invariant!(dYdt, Y, nothing, 0.0);

# run!
using OrdinaryDiffEq
timeend = 5000.0
Δt = 0.2
function make_dss_func()
  _dss!(x::Fields.Field)=Spaces.weighted_dss!(x)
  _dss!(::Any)=nothing
  dss_func(Y,t,integrator) = foreach(_dss!,Fields._values(Y))
  return dss_func
end
dss_func = make_dss_func()
dss_callback = FunctionCallingCallback(dss_func, func_start=true)
prob = ODEProblem(rhs_invariant!, Y, (0.0, timeend))
integrator = OrdinaryDiffEq.init(
    prob,
    SSPRK33(),
    dt = Δt,
    saveat = 10.0,
    progress = true,
    progress_message = (dt, u, p, t) -> t,
    callback = dss_callback
);

if haskey(ENV, "CI_PERF_SKIP_RUN") # for performance analysis
    throw(:exit_profile)
end

sol = @timev OrdinaryDiffEq.solve!(integrator)

ENV["GKSwstype"] = "nul"
import Plots, ClimaCorePlots
Plots.GRBackend()

dir = "dc_invariant_etot_topo"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

anim = Plots.@animate for u in sol.u
    Plots.plot(u.Yc.ρe ./ u.Yc.ρ)
end
Plots.mp4(anim, joinpath(path, "total_energy.mp4"), fps = 20)

If2c = Operators.InterpolateF2C()
anim = Plots.@animate for u in sol.u
    ᶜuw = @. Geometry.Covariant13Vector(u.uₕ) +
       Geometry.Covariant13Vector(If2c(u.w))
    w = @. Geometry.project(Geometry.WAxis(), ᶜuw)
    Plots.plot(w)
end
Plots.mp4(anim, joinpath(path, "vel_w.mp4"), fps = 20)

anim = Plots.@animate for u in sol.u
    ᶜuw = @. Geometry.Covariant13Vector(u.uₕ) +
       Geometry.Covariant13Vector(If2c(u.w))
    u = @. Geometry.project(Geometry.UAxis(), ᶜuw)
    Plots.plot(u)
end
Plots.mp4(anim, joinpath(path, "vel_u.mp4"), fps = 20)


anim = Plots.@animate for u in sol.u
           e = u.Yc.ρe ./ u.Yc.ρ
           uₕ = Geometry.UWVector.(Geometry.Covariant13Vector.(u.uₕ))
           w = Geometry.UWVector.(Geometry.Covariant13Vector.(If2c.(u.w)))
           uw = uₕ .+ w
           z = Fields.coordinate_field(uₕ).z
           I = @. e - Φ(z) - (norm(uw)^2) / 2
           T = @. I / C_v + T_0
           p = @. u.Yc.ρ * R_d * T
           θ = @. (MSLP/p)^(R_d/C_p)*(T)
           Plots.plot(θ)
end
Plots.mp4(anim, joinpath(path, "theta.mp4"), fps = 20)

# post-processing
Es = [sum(u.Yc.ρe) for u in sol.u]
Mass = [sum(u.Yc.ρ) for u in sol.u]

@show maximum(parent(Geometry.UVector.(sol[end].uₕ)))
@show minimum(parent(Geometry.UVector.(sol[end].uₕ)))
@show maximum(parent(Geometry.WVector.(sol[end].w)))
@show minimum(parent(Geometry.WVector.(sol[end].w)))

Plots.png(
    Plots.plot((Es .- energy_0) ./ energy_0),
    joinpath(path, "energy_cons.png"),
)
Plots.png(
    Plots.plot((Mass .- mass_0) ./ mass_0),
    joinpath(path, "mass_cons.png"),
)

function linkfig(figpath, alt = "")
    # buildkite-agent upload figpath
    # link figure in logs if we are running on CI
    if get(ENV, "BUILDKITE", "") == "true"
        artifact_url = "artifact://$figpath"
        print("\033]1338;url='$(artifact_url)';alt='$(alt)'\a\n")
    end
end

linkfig(
    relpath(joinpath(path, "energy_cons.png"), joinpath(@__DIR__, "../..")),
    "Total Energy",
)
linkfig(
    relpath(joinpath(path, "mass_cons.png"), joinpath(@__DIR__, "../..")),
    "Mass",
)
