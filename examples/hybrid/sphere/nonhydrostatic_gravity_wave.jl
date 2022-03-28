using ClimaCorePlots, Plots, ClimaCoreVTK
using ClimaCore.DataLayouts

include("baroclinic_wave_utilities.jl")

const sponge = false

using OrdinaryDiffEq: ODEProblem, solve, SSPRK33

import Logging
import TerminalLoggers

Logging.global_logger(TerminalLoggers.TerminalLogger())

# Nonhydrostatic gravity wave
# Reference: https://climate.ucdavis.edu/pubs/UJ2012JCP.pdf Section 5.4

const N = 0.01 # Brunt-Vaisala frequency
const S = grav^2 / cp_d / N^2
const T_0_nhw = 300 # isothermal atmospheric temperature
const Δθ = 10.0 # maximum potential temperature perturbation
const R_t = R / 3 # width of the perturbation
const L_z = 20.0e3 # vertial wave length of the perturbation
const p_0 = 1.0e5 # reference pressure
const λ_c_nhw = 180.0 # center longitude of the cosine bell
const ϕ_c_nhw = 0.0 # center latitude of the cosine bell

r(λ, ϕ) = R * acos(sind(ϕ_c_nhw) * sind(ϕ) + cosd(ϕ_c_nhw) * cosd(ϕ) * cosd(λ - λ_c_nhw))

# Variables required for driver.jl (modify as needed)
helems, zelems, npoly = 4, 10, 4
number_of_days = 5.0
t_end = FT(60 * 60 * 24 * number_of_days)
dt = FT(400)
dt_save_to_sol = FT(60 * 60 * 1/4)
dt_save_to_disk = FT(0) # 0 means don't save to disk
ode_algorithm = OrdinaryDiffEq.Rosenbrock23
jacobian_flags = (; ∂ᶜ𝔼ₜ∂ᶠ𝕄_mode = :no_∂ᶜp∂ᶜK, ∂ᶠ𝕄ₜ∂ᶜρ_mode = :exact)

horzdomain = Domains.SphereDomain(R)
vertdomain = Domains.IntervalDomain(
    Geometry.ZPoint{FT}(FT(0)),
    Geometry.ZPoint{FT}(FT(12e3));
    boundary_tags = (:bottom, :top),
)
horzmesh = Meshes.EquiangularCubedSphere(horzdomain, helems)
vertmesh = Meshes.IntervalMesh(vertdomain, nelems = zelems)
quad = Spaces.Quadratures.GLL{npoly + 1}()

Nv = Meshes.nelements(vertmesh)
Nf_center, Nf_face = 4, 1
vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)

if usempi
    horztopology = Topologies.DistributedTopology2D(horzmesh, Context)
    comms_ctx =
        Spaces.setup_comms(Context, horztopology, quad, Nv + 1, Nf_center)
    global_topology = Topologies.Topology2D(horzmesh)
    global_horz_space = Spaces.SpectralElementSpace2D(global_topology, quad)
    global_center_space = Spaces.ExtrudedFiniteDifferenceSpace(
        global_horz_space,
        vert_center_space,
    )
    global_face_space =
        Spaces.FaceExtrudedFiniteDifferenceSpace(global_center_space)

else
    horztopology = Topologies.Topology2D(horzmesh)
    comms_ctx = nothing
end

horzspace = Spaces.SpectralElementSpace2D(horztopology, quad, comms_ctx)

hv_center_space =
    Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
hv_face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)


function initial_condition(ϕ, λ, z)
    if rd < R_t
        s = 0.5 * (1 + cos(pi * rd / R_t))
    else
        s = 0.0
    end
    p = p_0 * (1 - S / T_0 + S / T_0 * exp(-N^2 * z / grav))^(cp_d / R_d)
    θ = T_0 * exp(N^2 * z / grav) + Δθ * s * sin(2 * pi * z / L_z)
    T = θ * (p / p_0)^κ
    ρ = p / R_d / T
    e = cv_d * (T - T_tri) + grav * z
    ρe = ρ * e

    return (ρ = ρ, ρe = ρe)
end

additional_cache(ᶜlocal_geometry, ᶠlocal_geometry, dt) = merge(
    hyperdiffusion_cache(ᶜlocal_geometry, ᶠlocal_geometry; κ₄ = FT(2e17)),
    sponge ? rayleigh_sponge_cache(ᶜlocal_geometry, ᶠlocal_geometry, dt) : (;),
    held_suarez_cache(ᶜlocal_geometry),
)
function additional_tendency!(Yₜ, Y, p, t, comms_ctx = nothing)
    hyperdiffusion_tendency!(Yₜ, Y, p, t, comms_ctx)
    sponge && rayleigh_sponge_tendency!(Yₜ, Y, p, t)
    held_suarez_tendency!(Yₜ, Y, p, t)
end

center_initial_condition(local_geometry) = center_initial_condition(local_geometry, Val(:ρe), GravityWave())

function postprocessing(sol, p, output_dir, usempi = false)
    sol_global = []
    if usempi
        for sol_step in sol.u
            sol_step_values_center_global =
                DataLayouts.gather(comms_ctx, Fields.field_values(sol_step.c))
            sol_step_values_face_global =
                DataLayouts.gather(comms_ctx, Fields.field_values(sol_step.f))
            if ClimaComms.iamroot(Context)
                sol_step_global = Fields.FieldVector(
                    c = Fields.Field(
                        sol_step_values_center_global,
                        global_center_space,
                    ),
                    f = Fields.Field(
                        sol_step_values_face_global,
                        global_face_space,
                    ),
                )
                push!(sol_global, sol_step_global)
            end
        end
        if ClimaComms.iamroot(Context)
        end
    else
        sol_global = sol.u
    end

    if !usempi || (usempi && ClimaComms.iamroot(Context))
        @info "L₂ norm of ρe at t = $(sol.t[1]): $(norm(sol_global[1].c.ρe))"
        @info "L₂ norm of ρe at t = $(sol.t[end]): $(norm(sol_global[end].c.ρe))"

        anim = Plots.@animate for Y in sol_global
            ᶜv = Geometry.UVVector.(Y.c.uₕ).components.data.:2
            Plots.plot(ᶜv, level = 1, clim = (-10,10))
        end
        Plots.mp4(anim, joinpath(output_dir, "v.mp4"), fps = 5)
        
        anim = Plots.@animate for Y in sol_global
            ᶜu = Geometry.UVVector.(Y.c.uₕ).components.data.:1
            Plots.plot(ᶜu, level = 1, clim = (-10,10))
        end
        Plots.mp4(anim, joinpath(output_dir, "u.mp4"), fps = 5)
        
        anim = Plots.@animate for Y in sol_global
            ᶠw = Geometry.WVector.(Y.f.w).components.data.:1
            ᶜw = @. ᶜinterp(ᶠw)
            Plots.plot(ᶜw, level = 1, clim = (-0.005,0.005))
        end
        Plots.mp4(anim, joinpath(output_dir, "w.mp4"), fps = 5)
        
        anim = Plots.@animate for Y in sol_global
            ᶜρ = Y.c.ρ
            ᶜu = Geometry.UVVector.(Y.c.uₕ).components.data.:1
            ᶜv = Geometry.UVVector.(Y.c.uₕ).components.data.:2
            ᶠw = Geometry.WVector.(Y.f.w).components.data.:1
            ᶜw = ᶜinterp.(ᶠw)
            ᶜuvw = @. Geometry.UVWVector(Y.c.uₕ) + Geometry.UVWVector(ᶜinterp(Y.f.w))
            ᶜz = ᶜlocal_geometry.coordinates.z
            eint = @. Y.c.ρe / ᶜρ - (grav * ᶜz) - 1/2 * norm_sqr(ᶜuvw)
            T = @. eint / cv_d + T_tri 
            Plots.plot(T .- 300, level = 1)
        end
        Plots.mp4(anim, joinpath(output_dir, "T.mp4"), fps = 5)
        
        
        vtk_counter = 0
        for Y in sol_global
            vtk_counter += 1
            ᶜρ = Y.c.ρ
            ᶜu = Geometry.UVVector.(Y.c.uₕ).components.data.:1
            ᶜv = Geometry.UVVector.(Y.c.uₕ).components.data.:2
            ᶠw = Geometry.WVector.(Y.f.w).components.data.:1
            ᶜw = ᶜinterp.(ᶠw)
            ᶜuvw = @. Geometry.UVWVector(Y.c.uₕ) + Geometry.UVWVector(ᶜinterp(Y.f.w))
            ᶜz = ᶜlocal_geometry.coordinates.z
            eint = @. Y.c.ρe / ᶜρ - (grav * ᶜz) - 1/2 * norm_sqr(ᶜuvw)
            T = @. eint / cv_d + T_tri 
            ClimaCoreVTK.writevtk(joinpath(output_dir,"nhw_$(vtk_counter)"), (Tprime=T, uh1 = ᶜu, uh2 = ᶜv, w = ᶜw))
        end
    end
end
