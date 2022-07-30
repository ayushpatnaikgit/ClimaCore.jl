# This was collected by running ClimaCore's test suite
# and printing all input argument types for `transform`
#! format: off
function used_transform_arg_types(::Type{FT}) where {FT}
    return [
        (CovariantAxis{(1, 2)}, Covariant1Vector{FT}),
        (CovariantAxis{(1, 2)}, Covariant13Vector{FT}),
        (CovariantAxis{(1, 2)}, AxisTensor{FT, 2, Tuple{CovariantAxis{(1,)}, CartesianAxis{(1,)}}, StaticArraysCore.SMatrix{1, 1, FT, 1},},),
        (CovariantAxis{(1, 2)}, AxisTensor{FT, 2, Tuple{CovariantAxis{(1, 3)}, CartesianAxis{(1,)}}, StaticArraysCore.SMatrix{2, 1, FT, 2},},),
        (LocalAxis{(1, 2)}, UVVector{FT}),
        (CovariantAxis{(1, 3)}, Covariant1Vector{FT}),
        (LocalAxis{(1, 3)}, UWVector{FT}),
        (CovariantAxis{(1,)}, Covariant13Vector{FT}),
        (LocalAxis{(1, 3)}, UVector{FT}),
        (CovariantAxis{(1, 3)}, Covariant13Vector{FT}),
        (LocalAxis{(1, 3)}, WVector{FT}),
        (CovariantAxis{(3,)}, Covariant13Vector{FT}),
        (CovariantAxis{(1, 2, 3)}, Covariant12Vector{FT}),
        (LocalAxis{(1, 2, 3)}, UVWVector{FT}),
        (CovariantAxis{(1, 2)}, Covariant123Vector{FT}),
        (CovariantAxis{(1, 2)}, Covariant12Vector{FT}),
        (LocalAxis{(1,)}, UVVector{FT}),
        (LocalAxis{(1,)}, UVector{FT}),
        (CovariantAxis{(1,)}, Covariant1Vector{FT}),
        (LocalAxis{(3,)}, UVVector{FT}),
        (LocalAxis{(3,)}, WVector{FT}),
        (CovariantAxis{(3,)}, Covariant3Vector{FT}),
        (CovariantAxis{(1, 2)}, Covariant3Vector{FT}),
        (LocalAxis{(1, 2, 3)}, UVVector{FT}),
        (CovariantAxis{(1, 2, 3)}, Covariant123Vector{FT}),
        (CovariantAxis{(1, 2)}, AxisTensor{FT, 2, Tuple{CovariantAxis{(1, 2)}, LocalAxis{(1, 2)}}, StaticArraysCore.SMatrix{2, 2, FT, 4},},),
        (LocalAxis{(1, 2)}, AxisTensor{FT, 2, Tuple{LocalAxis{(1, 2)}, LocalAxis{(1, 2)}}, StaticArraysCore.SMatrix{2, 2, FT, 4},},),
        (ContravariantAxis{(1, 2)}, Contravariant12Vector{FT}),
        (CovariantAxis{(1, 2, 3)}, Covariant3Vector{FT}),
        (ContravariantAxis{(1, 2, 3)}, Contravariant3Vector{FT}),
        (ContravariantAxis{(1, 2, 3)}, Contravariant123Vector{FT}),
        (LocalAxis{(1, 2, 3)}, WVector{FT}),
        (CovariantAxis{(3,)}, Covariant123Vector{FT}),
        (CovariantAxis{(1, 3)}, Covariant3Vector{FT}),
        (LocalAxis{(3,)}, UWVector{FT}),
        (LocalAxis{(3,)}, UVWVector{FT}),
        (ContravariantAxis{(1, 3)}, Contravariant13Vector{FT}),
        (ContravariantAxis{(1,)}, Contravariant13Vector{FT}),
        (ContravariantAxis{(3,)}, Contravariant13Vector{FT}),
        (LocalAxis{(1,)}, UWVector{FT}),
        (ContravariantAxis{(1, 2, 3)}, Contravariant12Vector{FT}),
        (ContravariantAxis{(1, 2)}, Contravariant123Vector{FT}),
        (ContravariantAxis{(3,)}, Contravariant123Vector{FT}),
        (LocalAxis{(1, 2)}, UVWVector{FT}),
        (LocalAxis{(1, 3)}, UVVector{FT}),
        (ContravariantAxis{(1, 2)}, Contravariant2Vector{FT}),
        (CovariantAxis{(1, 3)}, Covariant12Vector{FT}),
        (LocalAxis{(1, 2)}, UWVector{FT}),
    ]
end
#! format: on
