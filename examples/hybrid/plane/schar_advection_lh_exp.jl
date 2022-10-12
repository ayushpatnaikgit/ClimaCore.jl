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
    Operators,
    Hypsography
using ClimaCore.Geometry
using ClimaCore.Utilities: half
using ClimaCore.Meshes: GeneralizedExponentialStretching

using DiffEqCallbacks

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())


### This file follows the test problem described in 
# https://doi.org/10.1175/1520-0493(2002)130<2459:ANTFVC>2.0.CO;2
# Section 3(b)

const MSLP = 1e5 # mean sea level pressure
const grav = 9.8 # gravitational constant
const R_d = 287.058 # R dry (gas constant / mol mass dry air)
const γ = 1.4 # heat capacity ratio
const C_p = R_d * γ / (γ - 1) # heat capacity at constant pressure
const C_v = R_d / (γ - 1) # heat capacity at constant volume
const T_0 = 273.16 # triple point temperature
const kinematic_viscosity = 0.0 #m²/s
const hyperdiffusivity = 2e7 #m²/s

const nx = 40
const nz = 30
const np = 4
const Lx = 120000
const Lz = 25000
 
function warp_schar(coord)   
  x = Geometry.component(coord,1)
  FT = eltype(x)
  a = 5000
  λ = 4000
  h₀ = 250.0
  if abs(x) <= a
    h = h₀ * exp(-(x/a)^2) * (cos(π*x/λ))^2
  else
    h = FT(0)
  end
end

function warp_agnesi(
    coord
)       
    x = coord.x
    h₀ = 1
    a = 10_000
    return h₀ * a^2 / (x^2 + a^2)
end

function warp_trig(
    coord
)       
    x = coord.x
    h₀ = 2000
    a = 10_000
    return h₀ * cospi.(x / Lx).^2
end


function hvspace_2D(
    xlim = (-π, π),
    zlim = (0, 4π),
    xelem = nx,
    zelem = nz,
    npoly = np,
    warp_fn = warp_schar,
)
    FT = Float64
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_names = (:bottom, :top),
    )
    #vertmesh = Meshes.IntervalMesh(vertdomain, GeneralizedExponentialStretching(500.0, 5000.0), nelems = zelem)
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
hv_center_space, hv_face_space = hvspace_2D((-Lx/2, Lx/2), (0, Lz))

Φ(z) = grav * z

# Reference: https://journals.ametsoc.org/view/journals/mwre/140/4/mwr-d-10-05073.1.xml, Section 5a
# Prognostic thermodynamic variable: Total Energy 
function init_advection_over_mountain(x, z)
    θ₀ = 250.0
    cp_d = C_p
    cv_d = C_v
    p₀ = MSLP
    g = grav
    
    𝒩 = g / sqrt(cp_d * θ₀)
    π_exner = @. exp(-g * z / (cp_d * θ₀))
    θ = @. θ₀ * exp(𝒩 ^2 * z / g)
    T = @. π_exner * θ # temperature
    ρ = @. p₀ / (R_d * θ) * (π_exner)^(cp_d/R_d)
    e = @. cv_d * (T - T_0) + Φ(z) + 200.0
    ρe = @. ρ * e
    ρq = @. ρ * 0.0
    return (ρ = ρ,
            ρe = ρe,
            ρq = ρq)
end

function initial_velocity(x, z)
  FT = eltype(x)
  return @. Geometry.UWVector(FT(20), FT(0))
end

# initial conditions
coords = Fields.coordinate_field(hv_center_space)
face_coords = Fields.coordinate_field(hv_face_space)

    Ic2f = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
# Assign initial conditions to cell center, cell face variables
# Group scalars (ρ, ρe) in Yc 
# Retain uₕ and w as separate components of velocity vector (primitive variables)
Yc = map(coord -> init_advection_over_mountain(coord.x, coord.z), coords)
uₕ_local = map(coord -> initial_velocity(coord.x, coord.z), coords)
w = map(_ -> Geometry.Covariant3Vector(0.0), face_coords)
uₕ = Geometry.Covariant1Vector.(uₕ_local)
    
# ==========
u₁_bc = Fields.level(Ic2f.(uₕ), ClimaCore.Utilities.half)
gⁱʲ = Fields.level(Fields.local_geometry_field(hv_face_space), ClimaCore.Utilities.half).gⁱʲ
g13 = gⁱʲ.components.data.:3
g11 = gⁱʲ.components.data.:1
g33 = gⁱʲ.components.data.:4
u₃_bc = Geometry.Covariant3Vector.(-1 .* g13 .* u₁_bc.components.data.:1 ./ g33)
apply_boundary_w = Operators.SetBoundaryOperator(bottom = Operators.SetValue(u₃_bc))
@. w = apply_boundary_w(w)
# ==========


const u_init = uₕ

ᶜlg = Fields.local_geometry_field(hv_center_space)
ᶠlg = Fields.local_geometry_field(hv_face_space)

Y = Fields.FieldVector(Yc = Yc, uₕ = uₕ, w = w)

Spaces.weighted_dss!(Y.Yc)
Spaces.weighted_dss!(Y.uₕ)
Spaces.weighted_dss!(Y.w)
Spaces.weighted_dss!(u_init)

energy_0 = sum(Y.Yc.ρe)
mass_0 = sum(Y.Yc.ρ)

function rayleigh_sponge(z;
                         z_sponge=10000.0,
                         z_max=25000.0,
                         α = 1,  # Relaxation timescale
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
                         x_sponge=20000,
                         x_max=30000,
                         α = 0.5,  # Relaxation timescale
                         τ = 0.5,
                         γ = 2.0)
    if abs(x) >= x_sponge
        r = (abs(x) - x_sponge) / (x_max - x_sponge)
        β_sponge = α * sinpi(τ * r)^γ
        return β_sponge
    else
        return eltype(x)(0)
    end
end

function rhs_invariant!(dY, Y, _, t)

    cρ = Y.Yc.ρ # scalar on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant1Vector on centers
    cρe = Y.Yc.ρe
    cρq = Y.Yc.ρq

    dρ = dY.Yc.ρ
    dw = dY.w
    duₕ = dY.uₕ
    dρe = dY.Yc.ρe
    dρq = dY.Yc.ρq
    z = coords.z
    fz = face_coords.z
    fx = face_coords.x

    # 0) update w at the bottom

    hdiv = Operators.Divergence()
    hwdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    hwgrad = Operators.WeakGradient()
    hcurl = Operators.Curl()

    # get u_cov at first interior cell center
    # constant extrapolation to bottom face 
    # apply as boundary condition on w for interpolation operator 

    If2c = Operators.InterpolateF2C()
    Ic2f = Operators.InterpolateC2F(
        bottom = Operators.Extrapolate(),
        top = Operators.Extrapolate(),
    )
    
  f_upwind_product1 = Operators.UpwindBiasedProductC2F()
  f_upwind_product3 = Operators.Upwind3rdOrderBiasedProductC2F(
      bottom = Operators.FirstOrderOneSided(),
      top = Operators.FirstOrderOneSided(),
  )

    dρ .= 0 .* cρ

    cw = If2c.(fw)
    fuₕ = Ic2f.(cuₕ)

    # ==========
    u₁_bc = Fields.level(fuₕ,ClimaCore.Utilities.half)
    gⁱʲ = Fields.level(Fields.local_geometry_field(hv_face_space), ClimaCore.Utilities.half).gⁱʲ
    g13 = gⁱʲ.components.data.:3
    g11 = gⁱʲ.components.data.:1
    g33 = gⁱʲ.components.data.:4
    u₃_bc = Geometry.Covariant3Vector.(-1 .* g13 .* u₁_bc.components.data.:1 ./ g33)
    apply_boundary_w = Operators.SetBoundaryOperator(bottom = Operators.SetValue(u₃_bc))
    @. fw = apply_boundary_w(fw)
    # ==========

    cw = @. If2c(fw)
    cuw = @. Geometry.Covariant13Vector(cuₕ) + Geometry.Covariant13Vector(cw)
    fuw = @. Ic2f(cuw)
    

    ce = @. cρe / cρ
    cq = @. cρq / cρ
    cI = @. ce - Φ(z) - (norm(cuw)^2) / 2
    cT = @. cI / C_v + T_0
    cp = @. cρ * R_d * cT

    h_tot = @. ce + cp / cρ # Total enthalpy at cell centers

    ### HYPERVISCOSITY
    # 1) compute hyperviscosity coefficients
    χe = @. dρe = hwdiv(hgrad(h_tot)) # we store χe in dρe
    χq = @. dρq = hwdiv(hgrad(cq)) # we store χe in dρe
    χuₕ = @. duₕ = hwgrad(hdiv(cuₕ))

    Spaces.weighted_dss!(dρe)
    Spaces.weighted_dss!(duₕ)
    Spaces.weighted_dss!(dρq)

    κ₄_dynamic = hyperdiffusivity # m^4/s
    κ₄_tracer = hyperdiffusivity * 0 
    @. dρe = -κ₄_dynamic * hwdiv(cρ * hgrad(χe))
    @. dρq = -κ₄_tracer * hwdiv(cρ * hgrad(χq))
    @. duₕ = -κ₄_dynamic * (hwgrad(hdiv(χuₕ)))

    # 1) Mass conservation
    dw .= fw .* 0

    # 1.a) horizontal divergence
    dρ .-= hdiv.(cρ .* (cuw))

    # 1.b) vertical divergence

    vdivf2c = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.Contravariant3Vector.(0.0)),
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
    # effectively a homogeneous Neumann condition on u₁ at the boundary
    fω¹ = hcurl.(fw)
    fω¹_bc = Fields.level(fω¹,ClimaCore.Utilities.half)
    vcurlc2f = Operators.CurlC2F(
        bottom = Operators.SetCurl(fω¹_bc),
        top = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
    )
    fω¹ .+= vcurlc2f.(cuₕ)


#    ᶠρ_b = Fields.level(Ic2f.(cρ), ClimaCore.Utilities.half)
#    ᶠρ_t = Fields.level(Ic2f.(cρ), 30)
    
    # cross product
    # convert to contravariant
    fu¹ = @. Geometry.project(Geometry.Contravariant1Axis(), fuw) 
    fu³ = @. Geometry.project(Geometry.Contravariant3Axis(), fuw) 

    @. dw -= fω¹ × fu¹ # Covariant3Vector on faces
    @. duₕ -= If2c(fω¹ × fu³)
    @. duₕ -= hgrad(cp) / cρ

    vgradc2fP = Operators.GradientC2F(
      bottom = Operators.SetGradient(Geometry.Covariant3Vector(0.0)),
      top = Operators.SetGradient(Geometry.Covariant3Vector(0.0)),
    )
    
    fuw_bc = Fields.level(fuw, ClimaCore.Utilities.half)
    
    vgradc2fE = Operators.GradientC2F(
     #   bottom = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
     #   top = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
         bottom = Operators.SetValue(norm(fuw_bc)),
         top = Operators.SetGradient(Geometry.Covariant3Vector(0.0)),
    )
    vgradc2fΦ = Operators.GradientC2F(
        bottom = Operators.SetGradient(Geometry.Covariant3Vector(0.0)),
        top = Operators.SetGradient(Geometry.Covariant3Vector(0.0)),
    )
    
    cE = @. (norm(cuw)^2) / 2
    cΦ = @. Φ(z)
    @. duₕ -= (hgrad(cE) + hgrad(cΦ))
    @. dw -= (vgradc2fE(cE) + vgradc2fΦ(cΦ))
    
    @. dw -= vgradc2fP(cp) / Ic2f(cρ)

    # 3) total energy

    @. dρe -= hdiv(cuw * (cρe + cp))
    
    @. dρe -= vdivf2c(fw * Ic2f(cρe + cp)) 
    #@. dρe -= vdivf2c(Ic2f(cρ) * f_upwind_product1(fw, (cρe + cp)/cρ)) # Upwind Approximation - First Order
    #@. dρe -= vdivf2c(Ic2f(cρ) * f_upwind_product3(fw, (cρe + cp)/cρ)) # Upwind Approximation - Third Order
    
    @. dρe -= vdivf2c(Ic2f(cuₕ * (cρe + cp)))
    
    # 4) tracer tendencies  
    # In extruded grids
    @. dρq -= hdiv(cuw * (cρq))
    @. dρq -= vdivf2c(fw * Ic2f(cρq))
    @. dρq -= vdivf2c(Ic2f(cuₕ * (cρq)))

    # Uniform 2nd order diffusion
    ∂c = Operators.GradientF2C()
    fρ = @. Ic2f(cρ)
    κ₂ = kinematic_viscosity # m^2/s

    # Sponge tendency
    β = @. rayleigh_sponge(z)
    βf = @. rayleigh_sponge(fz)
    #@. duₕ -= β * (uₕ - u_init)
    @. dw -= βf * fw
    
    Spaces.weighted_dss!(dY.Yc)
    Spaces.weighted_dss!(dY.uₕ)
    Spaces.weighted_dss!(dY.w)

    return dY
end
dYdt = similar(Y);
rhs_invariant!(dYdt, Y, nothing, 0.0);

# run!
using OrdinaryDiffEq
Δt = min(Lx / nx / np / 300, Lz / nz / 300) * 0.50
@show Δt

timeend = 3600.0 * 20.0
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
    saveat = 500.0,
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

dir = "ScharMountain_LinearHydrostatic_HR"
path = joinpath(@__DIR__, "output", dir)
mkpath(path)

anim = Plots.@animate for u in sol.u
    Plots.plot(u.Yc.ρe ./ u.Yc.ρ)
end
Plots.mp4(anim, joinpath(path, "total_energy.mp4"), fps = 20)

anim = Plots.@animate for u in sol.u
    Plots.plot(u.Yc.ρq ./ u.Yc.ρ)
end
Plots.mp4(anim, joinpath(path, "tracer.mp4"), fps = 20)

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
    ᶜu = @. Geometry.Covariant13Vector(u.uₕ)
    ᶜw = @. Geometry.Covariant13Vector(If2c(u.w))
    w = @. Geometry.project(Geometry.Contravariant1Axis(), ᶜu) +  Geometry.project(Geometry.Contravariant1Axis(), ᶜw) 
    Plots.plot(w)
end
Plots.mp4(anim, joinpath(path, "ucontravariant1.mp4"), fps = 20)

anim = Plots.@animate for u in sol.u
    ᶜu = @. Geometry.Covariant13Vector(u.uₕ)
    ᶜw = @. Geometry.Covariant13Vector(If2c(u.w))
    w = @. Geometry.project(Geometry.Contravariant3Axis(), ᶜu) +  Geometry.project(Geometry.Contravariant3Axis(), ᶜw) 
    Plots.plot(w)
end
Plots.mp4(anim, joinpath(path, "contravariant3.mp4"), fps = 20)

# post-processing
Es = [sum(u.Yc.ρe) for u in sol.u]
Mass = [sum(u.Yc.ρ) for u in sol.u]

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
