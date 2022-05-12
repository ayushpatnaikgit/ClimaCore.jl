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
const hyperdiffusivity = 1e8*1.0 #m²/s
 
function warp_surface(coord)   
  x = Geometry.component(coord,1)
  FT = eltype(x)
  a = 25000
  λ = 8000
  h₀ = 3000
  if abs(x) <= a
    h = h₀ * (cos(π*x/2/a))^2 * (cos(π*x/λ))^2
  else
    h = FT(0)
  end
end

function hvspace_2D(
    xlim = (-π, π),
    zlim = (0, 4π),
    xelem = 75,
    zelem = 50,
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
hv_center_space, hv_face_space = hvspace_2D((-100000, 100000), (0, 25000))

Φ(z) = grav * z

# Reference: https://journals.ametsoc.org/view/journals/mwre/140/4/mwr-d-10-05073.1.xml, Section 5a
# Prognostic thermodynamic variable: Total Energy 
function init_advection_over_mountain(x, z)
    FT = eltype(x)
    θ_b = 300.0
    cp_d = C_p
    cv_d = C_v
    p_0 = MSLP
    g = grav

    π_exn = 1.0 - g * z / cp_d / θ_b # exner function
    T = π_exn * θ_b # temperature
    p = p_0 * π_exn^(cp_d / R_d) # pressure
    ρ = p / R_d / T # density
    uᵢ = 10.0
    z₁ = 4000.0
    z₂ = 5000.0
    u = FT(0)
    if z <= FT(z₁)
      u = FT(0)
    elseif z >= FT(z₂)
      u = FT(uᵢ)
    else
      u = FT(uᵢ) * sinpi(0.5 * (z - z₁) / (z₂ - z₁))
    end
    K = FT(0.5) * (u^2)
    e = cv_d * (T - T_0) + g * z  + K 
    ρe = ρ * e # total energy

    #x₀ = -50000.0
    x₀ = -11000.0
    z₀ = 9000.0
    A_x = 25000.0
    A_z = 3000.0
    r = @. sqrt((x-x₀)^2/A_x^2 + (z-z₀)^2/A_z^2)
    q₀ = 1.0
  
    if r <= 1
      q = q₀ * (cos(π*r/2))^2 
    else
      q = eltype(x)(q₀) * 0
    end
    ρq = @. ρ * q
    return (ρ = ρ,
            ρe = ρe, 
            ρq = ρq)
end

function initial_velocity(x, z)
  FT = eltype(x)
  uᵢ = 10.0
  z₁ = 4000.0
  z₂ = 5000.0
  if z<=z₁
    u = @. FT(0)
  elseif z >= z₂
    u = @. uᵢ
  else
    u = @. uᵢ * sin(π/2 * (z-z₁)/(z₂-z₁))
  end
  return @. Geometry.UWVector(u*0, FT(0))
end

# initial conditions
coords = Fields.coordinate_field(hv_center_space)
face_coords = Fields.coordinate_field(hv_face_space)

# Assign initial conditions to cell center, cell face variables
# Group scalars (ρ, ρe) in Yc 
# Retain uₕ and w as separate components of velocity vector (primitive variables)
Yc = map(coord -> init_advection_over_mountain(coord.x, coord.z), coords)
w = map(_ -> Geometry.Covariant3Vector(0.0), face_coords)
uₕ_local = map(coord -> initial_velocity(coord.x, coord.z), coords)
uₕ = Geometry.Covariant1Vector.(uₕ_local)

const u₀ = uₕ

ᶜlg = Fields.local_geometry_field(hv_center_space)
ᶠlg = Fields.local_geometry_field(hv_face_space)

Y = Fields.FieldVector(Yc = Yc, uₕ = uₕ, w = w)

energy_0 = sum(Y.Yc.ρe)
mass_0 = sum(Y.Yc.ρ)

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

    # 1.b) vertical divergence
    vdivf2c = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    vdivc2f = Operators.DivergenceC2F(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    vgradc2f = Operators.GradientC2F(
        bottom = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
        top = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
    )
    # curl term
    hcurl = Operators.Curl()
    # effectively a homogeneous Neumann condition on u₁ at the boundary
    vcurlc2f = Operators.CurlC2F(
        bottom = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
        top = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
    )

    dρ .= 0 .* cρ

    cw = If2c.(fw)
    fuₕ = Ic2f.(cuₕ)
    cuw = Geometry.Covariant13Vector.(cuₕ) .+ Geometry.Covariant13Vector.(cw)

    ce = @. cρe / cρ
    cq = @. cρq / cρ
    cI = @. ce - Φ(z) - (norm(cuw)^2) / 2
    cT = @. cI / C_v + T_0
    cpressure = @. cρ * R_d * cT
    h_tot = @. ce + cpressure / cρ # Total enthalpy at cell centers

    ### Horizontal Hyperviscosity
    # 1) compute hyperviscosity coefficients
    χe = @. dρe = hwdiv(hgrad(h_tot)) 
    χuₕ = @. duₕ = hwgrad(hdiv(cuw))
    χq = @. dρq = hwdiv(hgrad(cq)) 
    Spaces.weighted_dss!(dρe)
    Spaces.weighted_dss!(duₕ)
    Spaces.weighted_dss!(dρq)
    κ₄_dynamic = hyperdiffusivity 
    κ₄_tracer = hyperdiffusivity * 0
    @. dρe = -κ₄_dynamic * hwdiv(cρ * hgrad(χe))
    @. dρq = -κ₄_tracer * hwdiv(cρ * hgrad(χq))
    @. duₕ = -κ₄_dynamic * (hwgrad(hdiv(χuₕ)))

    # 1) MASS conservation
    dw .= fw .* 0
    # 1.a) horizontal divergence
    dρ .-= hdiv.(cρ .* (cuw))
    dρ .-= vdivf2c.(Ic2f.(cρ .* cuₕ))
    dρ .-= vdivf2c.(Ic2f.(cρ) .* fw)

##    # 2) Momentum equation
##    # MOMENTUM terms (See equations list in ClimaAtmos documentation)
    fuw = @. Geometry.Covariant13Vector(fuₕ) + Geometry.Covariant13Vector(fw)
    fω² = hcurl.(fw)
    fω² .+= vcurlc2f.(cuₕ)
##    # cross product
##    # convert to contravariant
    fu¹ = Geometry.project.(Ref(Geometry.Contravariant1Axis()), fuw)
    fu³ = Geometry.project.(Ref(Geometry.Contravariant3Axis()), fuw)
    @. duₕ -= If2c(fω² × fu³)
    @. dw -= fω² × fu¹ # Covariant3Vector on faces
#
#    # Pressure gradient in momentum equation
    @. duₕ -= hgrad(cpressure) / cρ
    @. dw -= vgradc2f(cpressure) / Ic2f(cρ)
    cE1 = @. (norm(cuw)^2) / 2 
    cE2 = @.  Φ(z)
    @. duₕ -= hgrad(cE1)
    @. duₕ -= hgrad(cE2)
    @. dw -= vgradc2f(cE1)
    @. dw -= vgradc2f(cE2)
#
#    # 3) total energy - 6 terms
    @. dρe -= hdiv(cuw * cρe)
    @. dρe -= hdiv(cuw * cpressure)
    @. dρe -= vdivf2c(fw * Ic2f(cρe))
    @. dρe -= vdivf2c(fw * Ic2f(cpressure))
    @. dρe -= vdivf2c(Ic2f(cuₕ * cρe))
    @. dρe -= vdivf2c(Ic2f(cuₕ * cpressure))
#    
#    # 4) tracer tendencies  - 3 terms
#    # In extruded grids
    @. dρq -= hdiv(cuw * (cρq))
    @. dρq -= vdivf2c(fw * Ic2f(cρq))
    @. dρq -= vdivf2c(Ic2f(cuₕ * (cρq)))
#
#    # Sponge tendency
    β = @. rayleigh_sponge(z)
    @. duₕ -= β * (uₕ - u₀)
    @. dw -= Ic2f(β) * fw
#
#    # Sets hyperdiffusion tendency on horizontal velocity to zero. 
#    # Restricts flow to that prescribed in initial conditions. 
#    # Tracer advection is only possible by flow specified at tstep =0.0
#    # Pressure gradients are assumed to be held constant by the initial 
#    # velocity / energy profiles
    #  @. duₕ *= 0.0  #(1)
    #  @. dw *= 0.0   #(2)
    #  @. dρe *= 0.0   #(3)
    #  @. dρ *= 0.0   #(4)
    #  @. dρq *= 0.0   #(5)

    Spaces.weighted_dss!(dY.Yc)
    Spaces.weighted_dss!(dY.uₕ)
    Spaces.weighted_dss!(dY.w)

    return dY
end

dYdt = similar(Y);
rhs_invariant!(dYdt, Y, nothing, 0.0);

# run!
using OrdinaryDiffEq
Δt = 0.75
timeend = 3Δt
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
    saveat = 500Δt,
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

dir = "tracer_decoupled_debug"
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
Plots.mp4(anim, joinpath(path, "contravariant1.mp4"), fps = 20)

anim = Plots.@animate for u in sol.u
    ᶜu = @. Geometry.Covariant13Vector(u.uₕ)
    ᶜw = @. Geometry.Covariant13Vector(If2c(u.w))
    w = @. Geometry.project(Geometry.Contravariant3Axis(), ᶜu) +  Geometry.project(Geometry.Contravariant3Axis(), ᶜw) 
    Plots.plot(w)
end
Plots.mp4(anim, joinpath(path, "contravariant3.mp4"), fps = 20)


anim = Plots.@animate for u in sol.u
    Y = u;
    cρ = Y.Yc.ρ # scalar on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant1Vector on centers
    cρe = Y.Yc.ρe
    cρq = Y.Yc.ρq
    z = coords.z

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

    # 1.b) vertical divergence
    vdivf2c = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    vdivc2f = Operators.DivergenceC2F(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    vgradc2f = Operators.GradientC2F(
        bottom = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
        top = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
    )
    # curl term
    hcurl = Operators.Curl()
    # effectively a homogeneous Neumann condition on u₁ at the boundary
    vcurlc2f = Operators.CurlC2F(
        bottom = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
        top = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
    )

    cw = @. If2c(fw)
    fuₕ = @. Ic2f(cuₕ)
    cuw = Geometry.Covariant13Vector.(cuₕ) .+ Geometry.Covariant13Vector.(cw)
    ce = @. cρe / cρ
    cq = @. cρq / cρ
    cI = @. ce - Φ(z) - (norm(cuw)^2) / 2
    cT = @. cI / C_v + T_0
    cpressure = @. cρ * R_d * cT
    h_tot = @. ce + cpressure / cρ # Total enthalpy at cell centers
    #######################################  

    χe = @. hwdiv(hgrad(h_tot)) 
    χuₕ = @. hwgrad(hdiv(cuw))
    χq = @. hwdiv(hgrad(cq)) 
    Spaces.weighted_dss!(χe)
    Spaces.weighted_dss!(χuₕ)
    Spaces.weighted_dss!(χq)
    κ₄_dynamic = hyperdiffusivity 
    κ₄_tracer = hyperdiffusivity * 0
    dρeh = @. -κ₄_dynamic * hwdiv(cρ * hgrad(χe))
    dρqh = @. -κ₄_tracer * hwdiv(cρ * hgrad(χq))
    duₕh = @. -κ₄_dynamic * (hwgrad(hdiv(χuₕ)))


    fuw = @. Geometry.Covariant13Vector(fuₕ) + Geometry.Covariant13Vector(fw)
    fω² = hcurl.(fw)
    fω² .+= vcurlc2f.(cuₕ)
    
    fu¹ = Geometry.project.(Ref(Geometry.Contravariant1Axis()), fuw)
    fu³ = Geometry.project.(Ref(Geometry.Contravariant3Axis()), fuw)

    cE1 = @. (norm(cuw)^2) / 2 
    cE2 = @.  Φ(z)
    
    dρ1 = @. hdiv.(cρ .* (cuw))
    dρ2 = @. vdivf2c.(Ic2f.(cρ .* cuₕ))
    dρ3 = @. vdivf2c.(Ic2f.(cρ) .* fw)

    duₕ1 = @. -If2c(fω² × fu³)
    duₕ2 = @. -hgrad(cpressure) / cρ
    duₕ3 = @. -hgrad(cE1)
    duₕ4 = @. -hgrad(cE2)
    
    dw1 = @. -fω² × fu¹ # Covariant3Vector on faces
    dw2 = @. -vgradc2f(cpressure) / Ic2f(cρ)
    dw3 = @. -vgradc2f(cE1)
    dw4 = @. -vgradc2f(cE2)

    dρe1 = @. -hdiv(cuw * cρe)
    dρe2 = @. -hdiv(cuw * cpressure)
    dρe3 = @. -vdivf2c(fw * Ic2f(cρe))
    dρe4 = @. -vdivf2c(fw * Ic2f(cpressure))
    dρe5 = @. -vdivf2c(Ic2f(cuₕ * cρe))
    dρe6 = @. -vdivf2c(Ic2f(cuₕ * cpressure))
    
    dρq1 = @. -hdiv(cuw * (cρq))
    dρq2 = @. -vdivf2c(fw * Ic2f(cρq))
    dρq3 = @. -vdivf2c(Ic2f(cuₕ * (cρq)))
    
    p1 = Plots.plot(duₕ1) # Make a line plot
    p2 = Plots.plot(duₕ2) # Make a line plot
    p3 = Plots.plot(duₕ3) # Make a line plot
    p4 = Plots.plot(duₕ4) # Make a line plot
    Plots.plot(p1, p2, p3, p4, layout = (2, 2), legend = false, size = (600,600))
end
Plots.mp4(anim, joinpath(path, "tendency_uh.mp4"), fps = 20)

anim = Plots.@animate for u in sol.u
    Y = u;
    cρ = Y.Yc.ρ # scalar on centers
    fw = Y.w # Covariant3Vector on faces
    cuₕ = Y.uₕ # Covariant1Vector on centers
    cρe = Y.Yc.ρe
    cρq = Y.Yc.ρq
    z = coords.z

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

    # 1.b) vertical divergence
    vdivf2c = Operators.DivergenceF2C(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    vdivc2f = Operators.DivergenceC2F(
        top = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
        bottom = Operators.SetValue(Geometry.Contravariant3Vector(0.0)),
    )
    vgradc2f = Operators.GradientC2F(
        bottom = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
        top = Operators.SetGradient(Geometry.Contravariant3Vector(0.0)),
    )
    # curl term
    hcurl = Operators.Curl()
    # effectively a homogeneous Neumann condition on u₁ at the boundary
    vcurlc2f = Operators.CurlC2F(
        bottom = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
        top = Operators.SetCurl(Geometry.Contravariant2Vector(0.0)),
    )

    cw = @. If2c(fw)
    fuₕ = @. Ic2f(cuₕ)
    cuw = Geometry.Covariant13Vector.(cuₕ) .+ Geometry.Covariant13Vector.(cw)
    ce = @. cρe / cρ
    cq = @. cρq / cρ
    cI = @. ce - Φ(z) - (norm(cuw)^2) / 2
    cT = @. cI / C_v + T_0
    cpressure = @. cρ * R_d * cT
    h_tot = @. ce + cpressure / cρ # Total enthalpy at cell centers
    #######################################  

    χe = @. hwdiv(hgrad(h_tot)) 
    χuₕ = @. hwgrad(hdiv(cuw))
    χq = @. hwdiv(hgrad(cq)) 
    Spaces.weighted_dss!(χe)
    Spaces.weighted_dss!(χuₕ)
    Spaces.weighted_dss!(χq)
    κ₄_dynamic = hyperdiffusivity 
    κ₄_tracer = hyperdiffusivity * 0
    dρeh = @. -κ₄_dynamic * hwdiv(cρ * hgrad(χe))
    dρqh = @. -κ₄_tracer * hwdiv(cρ * hgrad(χq))
    duₕh = @. -κ₄_dynamic * (hwgrad(hdiv(χuₕ)))


    fuw = @. Geometry.Covariant13Vector(fuₕ) + Geometry.Covariant13Vector(fw)
    fω² = hcurl.(fw)
    fω² .+= vcurlc2f.(cuₕ)
    
    fu¹ = Geometry.project.(Ref(Geometry.Contravariant1Axis()), fuw)
    fu³ = Geometry.project.(Ref(Geometry.Contravariant3Axis()), fuw)

    cE1 = @. (norm(cuw)^2) / 2 
    cE2 = @.  Φ(z)
    
    dρ1 = @. hdiv.(cρ .* (cuw))
    dρ2 = @. vdivf2c.(Ic2f.(cρ .* cuₕ))
    dρ3 = @. vdivf2c.(Ic2f.(cρ) .* fw)

    duₕ1 = @. -If2c(fω² × fu³)
    duₕ2 = @. -hgrad(cpressure) / cρ
    duₕ3 = @. -hgrad(cE1)
    duₕ4 = @. -hgrad(cE2)
    
    dw1 = @. -fω² × fu¹ # Covariant3Vector on faces
    dw2 = @. -vgradc2f(cpressure) / Ic2f(cρ)
    dw3 = @. -vgradc2f(cE1)
    dw4 = @. -vgradc2f(cE2)

    dρe1 = @. -hdiv(cuw * cρe)
    dρe2 = @. -hdiv(cuw * cpressure)
    dρe3 = @. -vdivf2c(fw * Ic2f(cρe))
    dρe4 = @. -vdivf2c(fw * Ic2f(cpressure))
    dρe5 = @. -vdivf2c(Ic2f(cuₕ * cρe))
    dρe6 = @. -vdivf2c(Ic2f(cuₕ * cpressure))
    
    dρq1 = @. -hdiv(cuw * (cρq))
    dρq2 = @. -vdivf2c(fw * Ic2f(cρq))
    dρq3 = @. -vdivf2c(Ic2f(cuₕ * (cρq)))
    
    p1 = Plots.plot(dw1) # Make a line plot
    p2 = Plots.plot(dw2) # Make a line plot
    p3 = Plots.plot(dw3) # Make a line plot
    p4 = Plots.plot(dw4) # Make a line plot
    Plots.plot(p1, p2, p3, p4, layout = (2, 2), legend = false, size = (600,600))
end
Plots.mp4(anim, joinpath(path, "tendency_w.mp4"), fps = 20)

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
