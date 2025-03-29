module Tools

using LinearAlgebra

"""
    sphere_points(n::Int64) -> Array{Float64,2}
Return n points sampled from the surface of a unit sphere with radius 0.4.
"""
function sphere_points(n::Int64)
    points = rand(3, n) .- 0.5

    # Normalize the points so that they correspond to a unit sphere of radius 0.4.
    for i in axes(points, 2)
        points[:, i] ./= norm(points[:, i])
        points[:, i] .*= 0.4
    end
    return points


end

end