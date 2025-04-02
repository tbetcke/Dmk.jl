

module Chebychev

using LinearAlgebra

import SIMD: Vec

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


function evaluate1d_impl(eval_points, values, cheb_points, weights)

    npoints = length(eval_points)
    n = length(values)

    @assert typeof(eval_points) <: AbstractVector
    @assert typeof(values) <: AbstractVector
    @assert typeof(cheb_points) <: AbstractVector
    @assert typeof(weights) <: AbstractVector
    @assert length(cheb_points) == n
    @assert length(weights) == n

    S = eltype(values)
    R = eltype(values)

    if S <: Vec
        R = eltype(S)
    end

    # If the elements of `values` contain a SIMD vector, we want the actual
    # underlying real type.

    @assert R <: Real
    @assert eltype(cheb_points) == R
    @assert eltype(weights) == R
    @assert eltype(eval_points) == R



    numerator = zeros(S, npoints)
    denominator = zeros(S, npoints)
    require_exact = zeros(Int64, npoints)

    for (index, cheb_point) in enumerate(cheb_points)
        @inbounds for i in 1:npoints
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

    @inbounds for index in 1:npoints
        if require_exact[index] > 0
            result[index] = values[require_exact[index]]
        end
    end

    result

end


"""Evaluate a 1d Chebychev interpolant."""
function evaluate1d(eval_points, values)

    @assert typeof(eval_points) <: AbstractVector
    @assert typeof(values) <: AbstractVector
    @assert length(eval_points) > 0
    @assert length(values) > 0

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
- `values::AbstractMatrix{S}`: The values of the interpolant at the Chebychev points. The first axis is the x dimension.
   The second axis is the y dimension.
"""
function evaluate2d_tensor(eval_x, eval_y, values)

    # The first axis is the x-dimension
    # So m is the number of Chebychev points in x direction
    m = size(values, 1)
    # n is the number of Chebychev points in the y direction.
    n = size(values, 2)

    npx = length(eval_x)
    npy = length(eval_y)

    S = eltype(values)

    cheb_points_m = Chebychev.cheb_points(m)
    cheb_points_n = Chebychev.cheb_points(n)
    weights_m = Chebychev.cheb_weights(m)
    weights_n = Chebychev.cheb_weights(n)

    res1 = zeros(S, m, npy)

    # For each x evaluate across the y direction. 
    for j in 1:m
        res1[j, :] = Chebychev.evaluate1d_impl(eval_y, values[j, :], cheb_points_n, weights_n)
    end

    res2 = zeros(S, npx, npy)

    # For each y point now evaluate across the x direction.

    for j in 1:npy
        res2[:, j] = Chebychev.evaluate1d_impl(eval_x, res1[:, j], cheb_points_m, weights_m)
    end

    return res2
end

"""
   evaluate3d(eval_x, eval_y, values)

# Evaluate a 3d chebychev interpolant on a 3d tensor grid.

# Arguments
# - `eval_x::AbstractVector{S}`: The x-coordinates where the interpolant is to be evaluated.
# - `eval_y::AbstractVector{S}`: The y-coordinates where the interpolant is to be evaluated.
# - `eval_z::AbstractVector{S}`: The z-coordinates where the interpolant is to be evaluated.
# - `values::AbstractArray{S, 3}`: The values of the interpolant at the Chebychev points. The first axis is the x dimension.
#    The second axis is the y dimension. The third is the z dimension.
# """
function evaluate3d_tensor(eval_x, eval_y, eval_z, values)

    # The first axis is the z-dimension
    # So m is the number of Chebychev points in x direction
    m = size(values, 1)
    # n is the number of Chebychev points in the y direction.
    n = size(values, 2)
    # p is the number of Chebychev points in the z direction.
    p = size(values, 3)

    npx = length(eval_x)
    npy = length(eval_y)
    npz = length(eval_z)

    S = eltype(values)


    cheb_points_m = Chebychev.cheb_points(m)
    cheb_points_n = Chebychev.cheb_points(n)
    weights_m = Chebychev.cheb_weights(m)
    weights_n = Chebychev.cheb_weights(n)
    cheb_points_p = Chebychev.cheb_points(p)
    weights_p = Chebychev.cheb_weights(p)

    res1 = zeros(S, m, n, npz)

    # For each x and y evaluate across the z direction. 
    for i in 1:m
        for j in 1:n
            res1[i, j, :] = Chebychev.evaluate1d_impl(eval_z, values[i, j, :], cheb_points_p, weights_p)
        end
    end

    res2 = zeros(S, m, npy, npz)

    # Now evaluate across y

    for i in 1:m
        for j in 1:npz
            res2[i, :, j] = Chebychev.evaluate1d_impl(eval_y, res1[i, :, j], cheb_points_n, weights_n)
        end
    end

    # We have now done the y, z plane for each x point. Now evaluate for each of those points along
    # the x direction.

    res3 = zeros(S, npx, npy, npz)

    for i in 1:npy
        for j in 1:npz
            res3[:, i, j] = Chebychev.evaluate1d_impl(eval_x, res2[:, i, j], cheb_points_m, weights_m)
        end
    end

    return res3
end


# function evaluate3d_single_point_impl(x::S, y::S, z::S, values::V, cheb_points::Vector{S}, weights::Vector{S})::S where {S<:Real,V<:AbstractArray{S, 3}}

#     function eval_1d_impl(x::S, values::V, cheb_points::Vector{S}, weights::Vector{S})::S where {S<:Real,V<:AbstractVector{S, 3}}::S

#         n = length(values)

#         numerator::S = 0.0
#         denominator::S = 0.0
#         exact_index::Int64 = 0

#         for (index, cheb_point) in enumerate(cheb_points)
#             diff = x - cheb_point
#             if iszero(diff)
#                 exact_index = index
#             end
#             diff_inv = 1.0 / diff
#             numerator += values[index] * weights[index] * diff_inv
#             denominator += weights[index] * diff_inv
#         end

#         result = numerator / denominator

#         if exact_index > 0
#             result = values[exact_index]
#         end

#         return result
#     end

#     # We now do the 3d evaluation.

#         # The first axis is the z-dimension
#     # So m is the number of Chebychev points in z direction
#     m = size(values, 1)
#     # n is the number of Chebychev points in the y direction.
#     n = size(values, 2)
#     # p is the number of Chebychev points in the x direction.
#     p = size(values, 3)



#     cheb_points_m = Chebychev.cheb_points(m)
#     cheb_points_n = Chebychev.cheb_points(n)
#     weights_m = Chebychev.cheb_weights(m)
#     weights_n = Chebychev.cheb_weights(n)
#     cheb_points_p = Chebychev.cheb_points(p)
#     weights_p = Chebychev.cheb_weights(p)

#     res1 = zeros(S, m, n, npx)

#     # For each y and z evaluate across the x direction. 
#     for i in 1:m
#         for j in 1:n
#             res1[i, j, :] = Chebychev.evaluate1d_impl(eval_x, values[i, j, :], cheb_points_p, weights_p)
#         end
#     end

#     res2 = zeros(S, m, npy, npx)

#     # For each z point now evaluate across all x along the y direction.

#     for i in 1:m
#         for j in 1:npx
#             res2[i, :, j] = Chebychev.evaluate1d_impl(eval_y, res1[i, :, j], cheb_points_n, weights_n)
#         end
#     end

#     # We have now done the x,y plane for each z point. Now evaluate for each of those points along
#     # the z direction.

#     res3 = zeros(S, npz, npy, npx)

#     for i in 1:npy
#         for j in 1:npx
#             res3[:, i, j] = Chebychev.evaluate1d_impl(eval_z, res2[:, i, j], cheb_points_m, weights_m)
#         end
#     end



# end

end
