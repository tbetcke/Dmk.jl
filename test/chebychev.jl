using Test
using Dmk.Octree.Morton
using Dmk.Octree.Constants
using TestItemRunner

@testitem "test chebychev polynomial 1d" begin
    import Dmk.Chebychev
    using LinearAlgebra

    n = 21
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

end

@testitem "test chebychev polynomial 2d" begin
    import Dmk.Chebychev
    using LinearAlgebra


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

    println("max error: ", maximum(abs.((expected_values - actual_values) ./ (expected_values))))

    @test maximum(abs.((expected_values - actual_values) ./ (expected_values))) < 1E-14

end

@testitem "test chebychev polynomial 3d" begin
    import Dmk.Chebychev
    using LinearAlgebra

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
    actual_values = Chebychev.evaluate3d_tensor(eval_x, eval_y, eval_z, values)
    @test maximum(abs.((expected_values - actual_values) ./ (expected_values))) < 1E-14
end