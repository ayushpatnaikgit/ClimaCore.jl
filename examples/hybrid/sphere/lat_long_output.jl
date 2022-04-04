import ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Operators
using NCDatasets
using ClimaCoreTempestRemap

using JLD2

# Required for vertical velocity interpolation to cell-centers
const ᶜinterp = Operators.InterpolateF2C()

datain = jldopen("nonhydrostatic_gravity_wave_day0.jld2")

Y = datain["Y"]; # integrator solution
t = datain["t"]; # integrator solution
Nq = 5
hspace = ClimaCore.Spaces.SpectralElementSpace2D(ClimaCore.Fields.axes(Y.c.ρ).horizontal_space.topology, ClimaCore.Spaces.Quadratures.GLL{Nq}())
vspace = ClimaCore.Spaces.CenterFiniteDifferenceSpace(ClimaCore.Fields.axes(Y.c.ρ).vertical_topology)
cspace = ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace(hspace, vspace)
fspace = ClimaCore.Spaces.FaceExtrudedFiniteDifferenceSpace(cspace)

dirname = mkpath("ncoutput")
datafile_cc = joinpath(dirname, "sphere_inertial_gravity_wave.nc")
nc = NCDataset(datafile_cc, "c")
# defines the appropriate dimensions and variables for a space coordinate
def_space_coord(nc, cspace, type="cgll")
# defines the appropriate dimensions and variables for a time coordinate (by default, unlimited size)
nc_time = def_time_coord(nc)

nc_rho = defVar(nc, "rho", Float64, cspace, ("time",))
nc_e = defVar(nc, "e", Float64, cspace, ("time",))
nc_u = defVar(nc, "u", Float64, cspace, ("time",))
nc_v = defVar(nc, "v", Float64, cspace, ("time",))
nc_w = defVar(nc, "w", Float64, cspace, ("time",))

for i = 1:length(Y.c)
    nc_time[i] = t[i]
    Yc = Y.c
    Yf = Y.f
    w = ᶜinterp.(Yf.w)
    nc_rho[:,i] = Yc.ρ
    nc_e[:,i] = Yc.ρe ./ Yc.ρ
    nc_u[:,i] = Yc.uₕ.components.data.:1
    nc_v[:,i] = Yc.uₕ.components.data.:2
    nc_w[:,i] = w.components.data.:1
end

close(nc)

# write out our cubed sphere mesh
meshfile_cc = joinpath(dirname, "mesh_cubedsphere.g") 
write_exodus(meshfile_cc, hspace.topology)

nlat = 90
nlon = 180
meshfile_rll = joinpath(dirname, "mesh_rll.g")
rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

meshfile_overlap = joinpath(dirname,"mesh_overlap.g")
overlap_mesh(meshfile_overlap, meshfile_cc, meshfile_rll)

weightfile = joinpath(dirname, "remap_weights.nc")
remap_weights(
    weightfile,
    meshfile_cc,
    meshfile_rll,
    meshfile_overlap;
    in_type = "cgll",
    in_np = Nq,
)

datafile_rll = joinpath(dirname,"data_rll.nc")
apply_remap(datafile_rll, datafile_cc, weightfile, ["rho", "e", "u", "v", "w"], verbose =true)
