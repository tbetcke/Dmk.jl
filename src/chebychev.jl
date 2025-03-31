

module Chebychev

using LinearAlgebra

"""Return the n Chebychev points on the interval [-1, 1]."""
function cheb_points(n::Int64)
    return [cos((j * pi / (n - 1))) for j in 0:n-1]
end

""" Return the n Chebychev weights for barycentric interpolation."""
function cheb_weights(n::Int64)
    weights = zeros(n)
    weights[1] = weights[n] = 0.5
    weights[2:n-1] .= 1.0
    return weights .* (-1.0) .^ (1:n)
end


function evaluate1d_impl(eval_points::T, values::V, cheb_points::Vector{S}, weights::Vector{S})::Vector{S} where {S<:Real,T<:AbstractVector{S},V<:AbstractVector{S}}

    npoints = length(eval_points)
    n = length(values)

    numerator = zeros(S, npoints)
    denominator = zeros(S, npoints)
    require_exact = zeros(Int64, npoints)

    for (index, cheb_point) in enumerate(cheb_points)
        for i in 1:npoints
            diff = eval_points[i] - cheb_point
            diff_inv = 1.0 / diff
            if iszero(diff)
                require_exact[i] = index
            end
            numerator[i] += values[index] * weights[index] * diff_inv
            denominator[i] += weights[index] * diff_inv
        end
    end

    result = numerator ./ denominator

    for index in 1:npoints
        if require_exact[index] > 0
            result[index] = values[require_exact[index]]
        end
    end

    result

end


"""Evaluate a 1d Chebychev interpolant."""
function evaluate1d(eval_points::T, values::V)::Vector{S} where {S<:Real,T<:AbstractVector{S},V<:AbstractVector{S}}
    n = length(values)
    cheb_points = Chebychev.cheb_points(n)
    weights = Chebychev.cheb_weights(n)
    return evaluate1d_impl(eval_points, values, cheb_points, weights)

end


"""
   evaluate2d(eval_x, eval_y, values)

Evaluate aq 2d chebychev interpolant on a 2d tensor grid.

# Arguments
- `eval_x::AbstractVector{S}`: The x-coordinates where the interpolant is to be evaluated.
- `eval_y::AbstractVector{S}`: The y-coordinates where the interpolant is to be evaluated.
- `values::AbstractMatrix{S}`: The values of the interpolant at the Chebychev points. The first axis is the y-dimension.
   The second axis is the x dimension.
"""
function evaluate2d(eval_x::TX, eval_y::TY, values::V)::Matrix{S} where {S<:Real,TX<:AbstractVector{S},TY<:AbstractVector{S},V<:AbstractMatrix{S}}

    # The first axis is the y-dimension
    # So m is the number of Chebychev points in y direction
    m = size(values, 1)
    # n is the number of Chebychev points in the x direction.
    n = size(values, 2)

    npx = length(eval_x)
    npy = length(eval_y)


    cheb_points_m = Chebychev.cheb_points(m)
    cheb_points_n = Chebychev.cheb_points(n)
    weights_m = Chebychev.cheb_weights(m)
    weights_n = Chebychev.cheb_weights(n)

    res1 = zeros(S, m, npx)

    # For each y evaluate across the x direction. 
    for j in 1:m
        res1[j, :] = Chebychev.evaluate1d_impl(eval_x, values[j, :], cheb_points_n, weights_n)
    end

    res2 = zeros(S, npy, npx)

    # For each x point now evaluate across the y direction.

    for j in 1:npx
        res2[:, j] = Chebychev.evaluate1d_impl(eval_y, res1[:, j], cheb_points_m, weights_m)
    end

    return res2
end

end