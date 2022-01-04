push!(LOAD_PATH, joinpath(@__DIR__, "..", ".."))

using Test, Plots
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
    Topographies
using ClimaCore.Geometry

# Given a set of quadrature points (x_in, y_in, z_in) 
# Compute the warped topography coordinates, and output
# the new (x,y,z) coordinates. In the examples shown below, 
# the topography warp linearly relaxes to the top of the domain,
# other choices for relaxation functions may be applied. 

function warp_agnesi_peak(
    coord;
    a = 500,
)
    return 8 * a^3 / (coord.x^2 + 4 * a^2)
end

function warp_schar(
    coord;
    a = 250, ## Half-width parameter [m]
    h₀ = 500, ## Peak height [m]
    λ = 100, ## Wavelength
)
    x = coord.x
    h_star = abs(x) <= a ? h₀ * (cospi((x) / 2a))^2 : 0
    h = h_star * (cospi(x / λ))^2
    return h
end

# set up function space
function hvspace_2D(
    xlim = (-π, π),
    zlim = (0, 4π),
    helem = 10,
    velem = 50,
    npoly = 4;
    stretch = Meshes.Uniform(),
    warp_fn = warp_agnesi_peak,
)
    # build vertical mesh information with stretching in [0, H]
    FT = Float64
    vertdomain = Domains.IntervalDomain(
        Geometry.ZPoint{FT}(zlim[1]),
        Geometry.ZPoint{FT}(zlim[2]);
        boundary_tags = (:bottom, :top),
    )

    vertmesh = Meshes.IntervalMesh(vertdomain, stretch, nelems = velem)
    vert_face_space = Spaces.FaceFiniteDifferenceSpace(vertmesh)
    # build horizontal mesh information
    horzdomain = Domains.IntervalDomain(
        Geometry.XPoint{FT}(xlim[1])..Geometry.XPoint{FT}(xlim[2]),
        periodic = true,
    )

    # Construct Horizontal Mesh + Space
    horzmesh = Meshes.IntervalMesh(horzdomain; nelems = helem)
    horztopology = Topologies.IntervalTopology(horzmesh)
    quad = Spaces.Quadratures.GLL{npoly + 1}()
    hspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    # Apply warp
    z_surface = warp_fn.(Fields.coordinate_field(hspace))
    f_space = Spaces.ExtrudedFiniteDifferenceSpace(
        hspace,
        vert_face_space,
        z_surface,
        Topographies.LinearAdaption(),
    )
    c_space = Spaces.CenterExtrudedFiniteDifferenceSpace(
                                                        f_space
                                                       )
    return (c_space,f_space)
end

ENV["GKSwstype"] = "nul"
import Plots
Plots.GRBackend()
dirname = "agnesi_2d_3h"
path = joinpath(@__DIR__, "output", dirname)
mkpath(path)
# post-processing
If2c = Operators.InterpolateF2C()
u = sol.u
p1 = Plots.plot(If2c.(u[end].ρw) ./ u[end].Yc.ρ, xlim = (-25000, 25000), ylim = (0, 12000))
p2 = Plots.plot(u[end].Yc.ρuₕ ./ u[end].Yc.ρ, xlim = (-25000, 25000), ylim = (0, 12000))
p3 = Plots.plot(u[end].Yc.ρθ ./ u[end].Yc.ρ, xlim = (-25000, 25000), ylim = (0, 12000))
Plots.savefig(p1, path*"/vel_w.png")
Plots.savefig(p2, path*"/vel_u.png")
Plots.savefig(p3, path*"/theta.png")
