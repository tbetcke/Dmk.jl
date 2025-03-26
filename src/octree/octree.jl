module Octree

# Top-level code for Octree module

include("constants.jl")
include("morton.jl")

struct AdaptiveOctree
    max_depth::Int64
    number_of_points::Int64
    points::Array{Float64,2}
    point_morton_indices::Vector{Morton.MortonKey}
end

function adaptive_octree(points::Array{Float64,2}, max_depth::Int64, max_points_per_box::Int64)

    @assert size(points, 1) == 3

    number_of_points = size(points, 2)
    point_morton_indices = Vector{Morton.MortonKey}(undef, number_of_points)

    for i in 1:number_of_points
        point_morton_indices[i] = Morton.from_physical_coordinates(points[:, i]..., max_depth)
    end

end


end