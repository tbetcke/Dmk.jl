using Test
using Dmk.Octree.Morton
using Dmk.Octree.Constants
using TestItemRunner

@testitem "test chebychev polynomial 1d" begin
    import Dmk.Chebychev
    using LinearAlgebra
    using SIMD: Vec

    n = 30
    cheb_points = Chebychev.cheb_points(n)
    cheb_weights = Chebychev.cheb_weights(n)
    values = exp.(-cheb_points .^ 2)
    eval_points = collect(range(-1, 1, length=1000))
    expected_values = exp.(-eval_points .^ 2)
    actual_values = Chebychev.evaluate1d(eval_points, values)
    @test maximum(abs.((expected_values - actual_values) ./ (expected_values))) < 1E-13

    eval_points = copy(cheb_points)
    actual_values = Chebychev.evaluate1d(eval_points, values)
    @test actual_values == values

    # Test SIMD

    values1 = exp.(-cheb_points .^ 2)
    values2 = exp.(-2.0 .* cheb_points .^ 2)
    simd_values = [Vec{2,Float64}((val1, val2)) for (val1, val2) in zip(values1, values2)]
    actual_values = Chebychev.evaluate1d(eval_points, simd_values)

    actual_values1 = [getindex(actual_value, 1) for actual_value in actual_values]
    actual_values2 = [getindex(actual_value, 2) for actual_value in actual_values]

    expected_values1 = exp.(-eval_points .^ 2)
    expected_values2 = exp.(-2.0 .* eval_points .^ 2)

    @test maximum(abs.((expected_values1 - actual_values1) ./ (expected_values1))) < 1E-13
    @test maximum(abs.((expected_values2 - actual_values2) ./ (expected_values2))) < 1E-13


end

@testitem "test chebychev polynomial 2d" begin
    import Dmk.Chebychev
    using LinearAlgebra
    import SIMD


    m = 40
    n = 50

    npy = 200
    npx = 100
    eval_x = range(-1, 1, npx)
    eval_y = range(-1, 1, npy)

    cheb_points_y = Chebychev.cheb_points(n)
    cheb_weights_y = Chebychev.cheb_weights(n)
    cheb_points_x = Chebychev.cheb_points(m)
    cheb_weights_x = Chebychev.cheb_weights(m)

    values = zeros(m, n)

    for i in axes(values, 1)
        for j in axes(values, 2)
            values[i, j] = exp.(-(3 * cheb_points_x[i]^2 + 5 * cheb_points_y[j]^2))

        end
    end

    expected_values = zeros(npx, npy)

    for i in axes(expected_values, 1)
        for j in axes(expected_values, 2)
            expected_values[i, j] = exp.(-(3 * eval_x[i]^2 + 5 * eval_y[j]^2))

        end
    end

    actual_values = Chebychev.evaluate2d_tensor(eval_x, eval_y, values)


    @test maximum(abs.((expected_values - actual_values) ./ (expected_values))) < 1E-14

    # Test SIMD

    values1 = zeros(Float64, m, n)
    values2 = zeros(Float64, m, n)
    for i in axes(values, 1)
        for j in axes(values, 2)
            values1[i, j] = exp.(-(3 * cheb_points_x[i]^2 + 5 * cheb_points_y[j]^2))
            values2[i, j] = exp.(-(4 * cheb_points_x[i]^2 + 2 * cheb_points_y[j]^2))
        end
    end


    simd_values = zeros(SIMD.Vec{2,Float64}, m, n)
    for i in axes(values, 1)
        for j in axes(values, 2)
            simd_values[i, j] = SIMD.Vec{2,Float64}((values1[i, j], values2[i, j]))
        end
    end

    actual_values = Chebychev.evaluate2d_tensor(eval_x, eval_y, simd_values)
    actual_values1 = zeros(npx, npy)
    actual_values2 = zeros(npx, npy)
    for i in axes(actual_values, 1)
        for j in axes(actual_values, 2)
            actual_values1[i, j] = getindex(actual_values[i, j], 1)
            actual_values2[i, j] = getindex(actual_values[i, j], 2)
        end
    end

    expected_values1 = zeros(npx, npy)
    expected_values2 = zeros(npx, npy)
    for i in axes(expected_values1, 1)
        for j in axes(expected_values1, 2)
            expected_values1[i, j] = exp.(-(3 * eval_x[i]^2 + 5 * eval_y[j]^2))
            expected_values2[i, j] = exp.(-(4 * eval_x[i]^2 + 2 * eval_y[j]^2))
        end
    end

    @test maximum(abs.((expected_values1 - actual_values1) ./ (expected_values1))) < 1E-14
    @test maximum(abs.((expected_values2 - actual_values2) ./ (expected_values2))) < 1E-14

end

@testitem "test chebychev polynomial 3d" begin
    import Dmk.Chebychev
    using LinearAlgebra
    import SIMD: Vec
    import Random

    m = 70
    n = 80
    p = 90


    npx = 100
    npy = 200
    npz = 80

    eval_x = range(-1, 1, npx)
    eval_y = range(-1, 1, npy)
    eval_z = range(-1, 1, npz)

    cheb_points_y = Chebychev.cheb_points(n)
    cheb_weights_y = Chebychev.cheb_weights(n)
    cheb_points_x = Chebychev.cheb_points(m)
    cheb_weights_x = Chebychev.cheb_weights(m)
    cheb_points_z = Chebychev.cheb_points(p)
    cheb_weights_z = Chebychev.cheb_weights(p)
    values = zeros(m, n, p)
    for i in axes(values, 1)
        for j in axes(values, 2)
            for k in axes(values, 3)
                values[i, j, k] = exp.(-(3 * cheb_points_x[i]^2 + 5 * cheb_points_y[j]^2 + 6 * cheb_points_z[k]^2))
            end
        end
    end

    expected_values = zeros(npx, npy, npz)
    for i in axes(expected_values, 1)
        for j in axes(expected_values, 2)
            for k in axes(expected_values, 3)
                expected_values[i, j, k] = exp.(-(3 * eval_x[i]^2 + 5 * eval_y[j]^2 + 6 * eval_z[k]^2))
            end
        end
    end
    actual_values = Chebychev.evaluate3d_tensor_1x(eval_x, eval_y, eval_z, values)
    @test maximum(abs.((expected_values - actual_values) ./ (expected_values))) < 1E-14





    # Test SIMD

    values1 = zeros(Float64, m, n, p)
    values2 = zeros(Float64, m, n, p)
    for i in axes(values, 1)
        for j in axes(values, 2)
            for k in axes(values, 3)
                values1[i, j, k] = exp.(-(3 * cheb_points_x[i]^2 + 5 * cheb_points_y[j]^2 + 1.5 * cheb_points_z[k]^2))
                values2[i, j, k] = exp.(-(4 * cheb_points_x[i]^2 + 2 * cheb_points_y[j]^2 + 2.5 * cheb_points_z[k]^2))
            end
        end
    end

    expected_values1 = zeros(npx, npy, npz)
    expected_values2 = zeros(npx, npy, npz)
    for i in axes(expected_values1, 1)
        for j in axes(expected_values1, 2)
            for k in axes(expected_values1, 3)
                expected_values1[i, j, k] = exp.(-(3 * eval_x[i]^2 + 5 * eval_y[j]^2 + 1.5 * eval_z[k]^2))
                expected_values2[i, j, k] = exp.(-(4 * eval_x[i]^2 + 2 * eval_y[j]^2 + 2.5 * eval_z[k]^2))
            end
        end
    end

    (actual_values1, actual_values2) = Chebychev.evaluate3d_tensor_2x(eval_x, eval_y, eval_z, values1, values2)

    @test maximum(abs.((expected_values1 - actual_values1) ./ (expected_values1))) < 1E-14
    @test maximum(abs.((expected_values2 - actual_values2) ./ (expected_values2))) < 1E-14

    # Test the non-tensor version of the evaluation

    Random.seed!()

    points = 2.0 .* rand(3, 100) .- 1.0

    actual_values = Chebychev.evaluate3d(points, values)
    expected_values = [
        exp.(-(3 * points[1, i]^2 + 5 * points[2, i]^2 + 6 * points[3, i]^2)) for i in axes(points, 2)
    ]

    @test maximum(abs.((expected_values - actual_values) ./ (expected_values))) < 1E-14

    # Check that we get the right non-tensorial result for evaluation on chebychev points

    m = 10
    n = 10
    p = 10

    cheb_points_y = Chebychev.cheb_points(n)
    cheb_weights_y = Chebychev.cheb_weights(n)
    cheb_points_x = Chebychev.cheb_points(m)
    cheb_weights_x = Chebychev.cheb_weights(m)
    cheb_points_z = Chebychev.cheb_points(p)
    cheb_weights_z = Chebychev.cheb_weights(p)
    values = zeros(m, n, p)
    for i in axes(values, 1)
        for j in axes(values, 2)
            for k in axes(values, 3)
                values[i, j, k] = exp.(-(3 * cheb_points_x[i]^2 + 5 * cheb_points_y[j]^2 + 6 * cheb_points_z[k]^2))
            end
        end
    end


    points = hcat(cheb_points_x, cheb_points_y, cheb_points_z)'

    actual_values = Chebychev.evaluate3d(points, values)
    expected_values = [
        exp.(-(3 * points[1, i]^2 + 5 * points[2, i]^2 + 6 * points[3, i]^2)) for i in axes(points, 2)
    ]
    @test maximum(abs.((expected_values - actual_values) ./ (expected_values))) == 0




end