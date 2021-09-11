using Test
using StaticArrays, IntervalSets, LinearAlgebra

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
    Operators
import ClimaCore.Domains.Geometry: Cartesian2DPoint, ⊗

function hvspace_2D(
    xlim = (-π, π),
    zlim = (0, 4π),
    helem = 10,
    velem = 64,
    npoly = 7,
)
    FT = Float64
    vertdomain = Domains.IntervalDomain(
        FT(zlim[1]),
        FT(zlim[2]);
        x3boundary = (:bottom, :top),
    )
    vertmesh = Meshes.IntervalMesh(vertdomain, nelems = velem)
    vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)

    horzdomain = Domains.RectangleDomain(
        xlim[1]..xlim[2],
        -0..0,
        x1periodic = true,
        x2boundary = (:a, :b),
    )
    horzmesh = Meshes.EquispacedRectangleMesh(horzdomain, helem, 1)
    horztopology = Topologies.GridTopology(horzmesh)

    quad = Spaces.Quadratures.GLL{npoly + 1}()
    horzspace = Spaces.SpectralElementSpace1D(horztopology, quad)

    hv_center_space =
        Spaces.ExtrudedFiniteDifferenceSpace(horzspace, vert_center_space)
    hv_face_space = Spaces.FaceExtrudedFiniteDifferenceSpace(hv_center_space)
    return (hv_center_space, hv_face_space)
end

# V_face = Covariant3Vector
# V_center = Covariant12Vector
# V = V_center + V_face

# divergence(V) = horz_div(V) + vert_div(V)
# divergence(V) = horz_div(V_center) + horz_div(V_face) + vert_div(V_center) + vert_div(V_face)


# 1) horz_div(V_center): project to Contravariant1 + Contravariant2, take spectral derivative
# 2) horz_div(V_face): project to Contravariant1 + Contravariant2, interpolate to center, take spectral derivative
#   - will be zero if orthogional geom
# 3) vert_div(V_center): project to Contravariant3, interpolate to face, take FD deriv
#   - will be zero if orthogional geom
# 4) vert_div(V_face): project to Contravariant3, take FD deriv





@testset "1D SE, 1D FV Extruded Domain ∇ ODE Solve vertical" begin

    hv_center_space, hv_face_space = hvspace_2D()
    V = Geometry.Cartesian3Vector.(ones(Float64, hv_face_space),)

    function rhs!(dudt, u, _, t)
        A = Operators.AdvectionC2C(
            bottom = Operators.SetValue(sin(-t)),
            top = Operators.Extrapolate(),
        )
        return @. dudt = -A(V, u)
    end

    U = sin.(Fields.coordinate_field(hv_center_space).z)
    dudt = zeros(eltype(U), hv_center_space)
    rhs!(dudt, U, nothing, 0.0)

    using OrdinaryDiffEq
    Δt = 0.01
    prob = ODEProblem(rhs!, U, (0.0, 2π))
    sol = solve(prob, SSPRK33(), dt = Δt)

    htopo = Spaces.topology(hv_center_space)
    for h in 1:Topologies.nlocalelems(htopo)
        sol_column_field = ClimaCore.column(sol.u[end], 1, 1, h)
        ref_column_field = ClimaCore.column(U, 1, 1, h)
        @test norm(sol_column_field .- ref_column_field) ≤ 5e-2
    end
end

@testset "1D SE, 1D FV Extruded Domain ∇ ODE Solve horzizontal" begin

    # Advection Equation
    # ∂_t f + c ∂_x f  = 0
    # the solution translates to the right at speed c,
    # so if you you have a periodic domain of size [-π, π]
    # at time t, the solution is f(x - c * t, y)
    # here c == 1, integrate t == 2π or one full period

    function rhs!(dudt, u, _, t)
        # horizontal divergence operator applied to all levels
        hdiv = Operators.Divergence()
        @. dudt = -hdiv(u * Geometry.Cartesian1Vector(1.0))
        Spaces.weighted_dss!(dudt)
        return dudt
    end

    hv_center_space, _ = hvspace_2D()
    U = sin.(Fields.coordinate_field(hv_center_space).x)
    dudt = zeros(eltype(U), hv_center_space)
    rhs!(dudt, U, nothing, 0.0)

    using OrdinaryDiffEq
    Δt = 0.01
    prob = ODEProblem(rhs!, U, (0.0, 2π))
    sol = solve(prob, SSPRK33(), dt = Δt)

    @test norm(U .- sol.u[end]) ≤ 5e-5
end

@testset "1D SE, 1D FV Extruded Domain ∇ ODE Solve diagonally" begin

    # Advection equation in Cartesian domain with
    # uₕ = (cₕ, 0), uᵥ = (0, cᵥ)
    # ∂ₜf + ∇ₕ⋅(uₕ * f) + ∇ᵥ⋅(uᵥ * f)  = 0
    # the solution translates diagonally to the top right and
    # at time t, the solution is f(x - cₕ * t, y - cᵥ * t)
    # here cₕ == cᵥ == 1, integrate t == 2π or one full period.
    # This is only correct if the solution is localized in the vertical
    # as we don't use periodic boundary conditions in the vertical.
    #
    # NOTE: the equation setup is only correct for Cartesian domains!

    hv_center_space, hv_face_space = hvspace_2D((-1, 1), (-1, 1))

    function rhs!(dudt, u, _, t)
        h = u.h
        dh = dudt.h

        # vertical advection no inflow at bottom
        # and outflow at top
        Ic2f = Operators.InterpolateC2F(top = Operators.Extrapolate())
        divf2c = Operators.DivergenceF2C(
            bottom = Operators.SetValue(Geometry.Cartesian3Vector(0.0)),
        )
        # only upward advection
        @. dh = -divf2c(Ic2f(h) * Geometry.Cartesian3Vector(1.0))

        # only horizontal advection
        hdiv = Operators.Divergence()
        @. dh += -hdiv(h * Geometry.Cartesian1Vector(1.0))
        Spaces.weighted_dss!(dh)

        return dudt
    end

    # initial conditions
    function h_init(x_init, z_init, A = 1.0, σ = 0.2)
        coords = Fields.coordinate_field(hv_center_space)
        h = map(coords) do coord
            A * exp(-((coord.x + x_init)^2 + (coord.z + z_init)^2) / (2 * σ^2))
        end

        return h
    end

    using OrdinaryDiffEq
    U = Fields.FieldVector(h = h_init(0.5, 0.5))
    Δt = 0.01
    t_end = 1.0
    prob = ODEProblem(rhs!, U, (0.0, t_end))
    sol = solve(prob, SSPRK33(), dt = Δt)

    h_end = h_init(0.5 - t_end, 0.5 - t_end)
    @test norm(h_end .- sol.u[end].h) ≤ 5e-2
end


@testset "1D SE, 1D FV Extruded Domain ∇ ODE Solve vertical" begin

    hv_center_space, hv_face_space = hvspace_2D()
    V = Geometry.Cartesian3Vector.(ones(Float64, hv_face_space),)

    function rhs_advect!(dudt, u, _, t)
        A = Operators.AdvectionC2C(
            bottom = Operators.SetValue(Geometry.Cartesian1Vector(sin(-t))),
            top = Operators.Extrapolate(),
        )
        @. dudt = -A(V, u)
    end

    U =
        Geometry.Cartesian1Vector.(
            sin.(Fields.coordinate_field(hv_center_space).z),
        )
    dudt = similar(U)
    rhs_advect!(dudt, U, nothing, 0.0)


    function rhs_div!(dUdt, U, _, t)
        Ic2f = Operators.InterpolateC2F(top = Operators.Extrapolate())
        divf2c = Operators.DivergenceF2C(
            bottom = Operators.SetValue(
                Geometry.Cartesian3Vector(1.0) ⊗
                Geometry.Cartesian1Vector(sin(-t)),
            ),
        )
        # only upward advection
        @. dUdt = -divf2c(V ⊗ Ic2f(U))

        hdiv = Operators.Divergence()
        @. dUdt += -hdiv(Geometry.Cartesian1Vector(1.0) ⊗ U)
        Spaces.weighted_dss!(dUdt)
        return dUdt
    end

    U =
        Geometry.Cartesian1Vector.(
            sin.(Fields.coordinate_field(hv_center_space).z),
        )
    dUdt = similar(U)
    rhs_div!(dUdt, U, nothing, 0.0)

end

@testset "1D SE, 1D FV Extruded Domain scalar diffusion" begin

    hv_center_space, hv_face_space = hvspace_2D((-5, 5), (-5, 5))

    K = 1.0
    function rhs_diff!(dhdt, h, K, t)

        gradc2f = Operators.GradientC2F(
            top = Operators.SetValue(0.0),
            bottom = Operators.SetValue(0.0),
        )
        divf2c = Operators.DivergenceF2C()
        @. dhdt = divf2c(K * gradc2f(h))

        hgrad = Operators.Gradient()
        hdiv = Operators.Divergence()

        @. dhdt += hdiv(K * hgrad(h))
        Spaces.weighted_dss!(dhdt)
    end

    h = map(Fields.coordinate_field(hv_center_space)) do coord
        exp(-(coord.x^2 + coord.z^2) / 2)
    end
    dhdt = similar(h)
    rhs_diff!(dhdt, h, K, 0.0)

end

@testset "1D SE, 1D FV Extruded Domain vector diffusion" begin

    hv_center_space, hv_face_space = hvspace_2D((-5, 5), (-5, 5))

    K = 1.0
    function rhs_diff!(dhdt, h, K, t)

        gradc2f = Operators.GradientC2F(
            top = Operators.SetValue(Geometry.Cartesian1Vector(0.0)),
            bottom = Operators.SetValue(Geometry.Cartesian1Vector(0.0)),
        )
        divf2c = Operators.DivergenceF2C()
        @. dhdt = divf2c(K * gradc2f(h))

        hgrad = Operators.Gradient()
        hdiv = Operators.Divergence()

        @. dhdt += hdiv(K * hgrad(h))
        Spaces.weighted_dss!(dhdt)
    end

    h = map(Fields.coordinate_field(hv_center_space)) do coord
        Geometry.Cartesian1Vector(exp(-(coord.x^2 + coord.z^2) / 2))
    end
    dhdt = similar(h)
    rhs_diff!(dhdt, h, K, 0.0)

end
