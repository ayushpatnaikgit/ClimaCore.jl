# These tests require running with `--check-bounds=[auto|no]`
using Test
using StaticArrays, IntervalSets
import ClimaCore
import ClimaCore.Utilities: PlusHalf
import ClimaCore.DataLayouts: IJFH
import ClimaCore:
    Fields, slab, Domains, Topologies, Meshes, Operators, Spaces, Geometry

using FastBroadcast
using LinearAlgebra: norm
using Statistics: mean
using ForwardDiff

function FieldFromNamedTuple(space, nt::NamedTuple)
    cmv(z) = nt
    return cmv.(Fields.coordinate_field(space))
end

include(joinpath(@__DIR__, "util_spaces.jl"))

# https://github.com/CliMA/ClimaCore.jl/issues/946
@testset "Allocations with broadcasting Refs" begin
    FT = Float64
    function foo!(Yx::Fields.Field)
        Yx .= Ref(1) .+ Yx
        return nothing
    end
    function foocolumn!(Yx::Fields.Field)
        Fields.bycolumn(axes(Yx)) do colidx
            Yx[colidx] .= Ref(1) .+ Yx[colidx]
            nothing
        end
        return nothing
    end
    for space in all_spaces(FT)
        (
            space isa Spaces.ExtrudedFiniteDifferenceSpace ||
            space isa Spaces.SpectralElementSpace1D ||
            space isa Spaces.SpectralElementSpace2D
        ) || continue
        Y = FieldFromNamedTuple(space, (; x = FT(2)))

        # Plain broadcast
        Yx = Y.x
        foo!(Yx) # compile first
        p = @allocated foo!(Yx)
        @test p == 0

        # bycolumn
        foocolumn!(Yx) # compile first
        p = @allocated foocolumn!(Yx)
        @test p == 0
    end
end

# https://github.com/CliMA/ClimaCore.jl/issues/949
@testset "Allocations with getproperty on FieldVectors" begin
    FT = Float64
    function allocs_test!(Y)
        x = Y.x
        fill!(x, 2.0)
        nothing
    end
    function callfill!(Y)
        fill!(Y, Ref((; x = 2.0)))
        nothing
    end
    for space in all_spaces(FT)
        Y = FieldFromNamedTuple(space, (; x = FT(2)))
        allocs_test!(Y)
        p = @allocated allocs_test!(Y)
        @test p == 0

        callfill!(Y)
        p = @allocated callfill!(Y)
        @test p == 0
    end
end

# https://github.com/CliMA/ClimaCore.jl/issues/963
@testset "Allocations StencilCoefs broadcasting" begin
    sc(::Type{FT}) where {FT} =
        Operators.StencilCoefs{-1, 1}((zero(FT), one(FT), zero(FT)))
    function allocs_test1!(Y)
        x = Y.x
        FT = Spaces.undertype(axes(x))
        I = sc(FT)
        x .= x .+ Ref(I)
        nothing
    end
    function allocs_test2!(Y)
        x = Y.x
        FT = Spaces.undertype(axes(x))
        I = sc(FT)
        IR = Ref(sc(FT))
        @. x += IR
        nothing
    end
    FT = Float64
    for space in all_spaces(FT)
        Y = FieldFromNamedTuple(space, (; x = sc(FT)))
        allocs_test1!(Y)
        p = @allocated allocs_test1!(Y)
        @test_broken p == 0
        allocs_test2!(Y)
        p = @allocated allocs_test2!(Y)
        @test_broken p == 0
    end
end

@testset "Allocations in @.. broadcasting" begin
    function allocs_test_field!(Y1, dt, Y2)
        x = Y1.x
        @.. x += 2.0
        nothing
    end
    function allocs_test_field_vector!(Y1, dt, Y2, Y3)
        @.. Y1 = Y2 + dt * Y3
        nothing
    end
    FT = Float32
    for space in all_spaces(FT)
        Y1 = FieldFromNamedTuple(space, (; x = FT(2.0), y = FT(2.0), z = FT(2.0)))
        Y2 = FieldFromNamedTuple(space, (; x = FT(2.0), y = FT(2.0), z = FT(2.0)))
        Y3 = FieldFromNamedTuple(space, (; x = FT(2.0), y = FT(2.0), z = FT(2.0)))
        dt = FT(2.0)
        allocs_test_field!(Y1, dt, Y2)
        p = @allocated allocs_test_field!(Y1, dt, Y2)
        @test p == 0
        allocs_test_field_vector!(Y1, dt, Y2, Y3)
        p = @allocated allocs_test_field_vector!(Y1, dt, Y2, Y3)
        @test p == 0
    end
end
