using Test
using StaticArrays, IntervalSets
import ClimaCore.DataLayouts: IJFH
import ClimaCore:
    Fields, slab, Domains, Topologies, Meshes, Operators, Spaces, Geometry

using LinearAlgebra: norm
using ForwardDiff

function spectral_space_2D(; n1 = 1, n2 = 1, Nij = 4)
    domain = Domains.RectangleDomain(
        Geometry.XPoint(-3.0) .. Geometry.XPoint(5.0),
        Geometry.YPoint(-2.0) .. Geometry.YPoint(8.0),
        x1periodic = false,
        x2periodic = false,
        x1boundary = (:east, :west),
        x2boundary = (:south, :north),
    )
    mesh = Meshes.EquispacedRectangleMesh(domain, n1, n2)
    grid_topology = Topologies.GridTopology(mesh)

    quad = Spaces.Quadratures.GLL{Nij}()
    space = Spaces.SpectralElementSpace2D(grid_topology, quad)
    return space
end

@testset "1×1 2D domain space" begin
    Nij = 4
    n1 = n2 = 1
    space = spectral_space_2D(n1 = n1, n2 = n2, Nij = Nij)

    field =
        Fields.Field(IJFH{ComplexF64, Nij}(ones(Nij, Nij, 2, n1 * n2)), space)

    @test sum(field) ≈ Complex(1.0, 1.0) * 8.0 * 10.0 rtol = 10eps()
    @test sum(x -> 3.0, field) ≈ 3 * 8.0 * 10.0 rtol = 10eps()
    @test norm(field) ≈ sqrt(2.0 * 8.0 * 10.0) rtol = 10eps()
    @test norm(field, 1) ≈ norm(Complex(1.0, 1.0)) * 8.0 * 10.0 rtol = 10eps()
    @test norm(field, Inf) ≈ norm(Complex(1.0, 1.0)) rtol = 10eps()

    @test extrema(real, field) == (1.0, 1.0)

    @test Operators.matrix_interpolate(field, 4) ≈
          [Complex(1.0, 1.0) for i in 1:(4 * n1), j in 1:(4 * n2)]


    field_sin = map(x -> sin((x.x) / 2), Fields.coordinate_field(space))
    M = Operators.matrix_interpolate(field_sin, 20)
    @test size(M) == (20, 20)  # 20 x 20 for a 1 element field

    real_field = field.re

    # test broadcasting
    res = field .+ 1
    @test parent(Fields.field_values(res)) == Float64[
        f == 1 ? 2 : 1 for i in 1:Nij, j in 1:Nij, f in 1:2, h in 1:(n1 * n2)
    ]

    res = field.re .+ 1
    @test parent(Fields.field_values(res)) ==
          Float64[2 for i in 1:Nij, j in 1:Nij, f in 1:1, h in 1:(n1 * n2)]

end

@testset "Broadcasting interception for tuple-valued fields" begin
    n1 = n2 = 1
    Nij = 4
    space = spectral_space_2D(n1 = n1, n2 = n2, Nij = Nij)

    nt_field = Fields.Field(
        IJFH{NamedTuple{(:a, :b), Tuple{Float64, Float64}}, Nij}(
            ones(Nij, Nij, 2, n1 * n2),
        ),
        space,
    )
    nt_sum = sum(nt_field)
    @test nt_sum isa NamedTuple{(:a, :b), Tuple{Float64, Float64}}
    @test nt_sum.a ≈ 8.0 * 10.0 rtol = 10eps()
    @test nt_sum.b ≈ 8.0 * 10.0 rtol = 10eps()
    @test norm(nt_field) ≈ sqrt(2.0 * 8.0 * 10.0) rtol = 10eps()

    # test scalar asignment
    nt_field.a .= 0.0
    @test sum(nt_field.a) == 0.0
end

@testset "Special case handling for broadcased norm to pass through space local geometry" begin
    space = spectral_space_2D()
    u = Geometry.Covariant12Vector.(ones(space), ones(space))
    @test norm.(u) ≈ hypot(4 / 8 / 2, 4 / 10 / 2) .* ones(space)
end

@testset "FieldVector" begin
    space = spectral_space_2D()
    u = Geometry.Covariant12Vector.(ones(space), ones(space))
    x = Fields.coordinate_field(space)
    y = [1.0, 2.0, 3.0]
    z = 1.0
    Y = Fields.FieldVector(u = u, k = (x = x, y = y, z = z))

    @test propertynames(Y) == (:u, :k)
    @test propertynames(Y.k) == (:x, :y, :z)
    @test Y.u === u
    @test Y.k.x === x
    @test Y.k.y === y
    @test Y.k.z === z

    Y1 = 2 .* Y
    @test parent(Y1.u) == 2 .* parent(u)
    @test parent(Y1.k.x) == 2 .* parent(x)
    @test Y1.k.y == 2 .* y
    @test Y1.k.z === 2 * z

    Y1 .= Y1 .+ 2 .* Y
    @test parent(Y1.u) == 4 .* parent(u)
    @test parent(Y1.k.x) == 4 .* parent(x)
    @test Y1.k.y == 4 .* y
    @test Y1.k.z === 4 * z

    Y.k.z = 3.0
    @test Y.k.z === 3.0
end

@testset "FieldVector basetype replacement and deepcopy" begin
    domain_z = Domains.IntervalDomain(
        Geometry.ZPoint(-1.0) .. Geometry.ZPoint(1.0),
        periodic = true,
    )
    mesh_z = Meshes.IntervalMesh(domain_z; nelems = 10)
    topology_z = Topologies.IntervalTopology(mesh_z)

    domain_x = Domains.IntervalDomain(
        Geometry.XPoint(-1.0) .. Geometry.XPoint(1.0),
        periodic = true,
    )
    mesh_x = Meshes.IntervalMesh(domain_x; nelems = 10)
    topology_x = Topologies.IntervalTopology(mesh_x)

    domain_xy = Domains.RectangleDomain(
        Geometry.XPoint(-1.0) .. Geometry.XPoint(1.0),
        Geometry.YPoint(-1.0) .. Geometry.YPoint(1.0),
        x1periodic = true,
        x2periodic = true,
    )
    mesh_xy = Meshes.EquispacedRectangleMesh(domain_xy, 10, 10)
    topology_xy = Topologies.GridTopology(mesh_xy)

    quad = Spaces.Quadratures.GLL{4}()

    space_vf = Spaces.CenterFiniteDifferenceSpace(topology_z)
    space_ifh = Spaces.SpectralElementSpace1D(topology_x, quad)
    space_ijfh = Spaces.SpectralElementSpace2D(topology_xy, quad)
    space_vifh = Spaces.ExtrudedFiniteDifferenceSpace(space_ifh, space_vf)
    space_vijfh = Spaces.ExtrudedFiniteDifferenceSpace(space_ijfh, space_vf)

    space2field(space) = map(
        coord -> (coord, Geometry.Covariant12Vector(1.0, 2.0)),
        Fields.coordinate_field(space),
    )

    Y = Fields.FieldVector(
        field_vf = space2field(space_vf),
        field_if = slab(space2field(space_ifh), 1),
        field_ifh = space2field(space_ifh),
        field_ijf = slab(space2field(space_ijfh), 1, 1),
        field_ijfh = space2field(space_ijfh),
        field_vifh = space2field(space_vifh),
        field_vijfh = space2field(space_vijfh),
        array = [1.0, 2.0, 3.0],
        scalar = 1.0,
    )

    Yf = ForwardDiff.Dual{Nothing}.(Y, 1.0)
    Yf .= Yf .^ 2 .+ Y
    @test all(ForwardDiff.value.(Yf) .== Y .^ 2 .+ Y)
    @test all(ForwardDiff.partials.(Yf, 1) .== 2 .* Y)

    dual_field = Yf.field_vf
    dual_field_original_basetype = similar(Y.field_vf, eltype(dual_field))
    @test eltype(dual_field_original_basetype) === eltype(dual_field)
    @test eltype(parent(dual_field_original_basetype)) === Float64
    @test eltype(parent(dual_field)) === ForwardDiff.Dual{Nothing, Float64, 1}

    object_that_contains_Yf = (; Yf)
    @test axes(deepcopy(Yf).field_vf) === space_vf
    @test axes(deepcopy(object_that_contains_Yf).Yf.field_vf) === space_vf
end
