module Octree

# Top-level code for Octree module

include("constants.jl")
include("morton.jl")

import .Morton: MortonKey
struct AdaptiveOctree
    max_depth::Int64
    number_of_points::Int64
    points::Array{Float64,2}
    point_morton_indices::Vector{Morton.MortonKey}
    refined_tree::Vector{Morton.MortonKey}
    leaf_to_points::Dict{Morton.MortonKey,Vector{Int64}}
    points_to_leafs::Vector{Morton.MortonKey}
end

"""
    adaptive_octree(points::Array{Float64,2}, max_depth::Int64, max_points_per_box::Int64) -> AdaptiveOctree

Constructs an adaptive octree from the given `points` with a maximum depth of `max_depth` and a maximum number of points per box of `max_points_per_box`.
"""
function adaptive_octree(points::Array{Float64,2}, max_depth::Int64, max_points_per_box::Int64)

    @assert size(points, 1) == 3

    number_of_points = size(points, 2)
    point_morton_indices = Vector{Morton.MortonKey}(undef, number_of_points)

    for i in 1:number_of_points
        point_morton_indices[i] = Morton.from_physical_coordinate(points[:, i]..., max_depth)
    end

    # Sort the points by their Morton keys.
    sorted_indices = sortperm(point_morton_indices)
    # Now generate the refined tree with the sorted order of keys.
    sorted_point_keys = point_morton_indices[sorted_indices]
    refined_tree = refine_tree(MortonKey(0), sorted_point_keys, max_depth, max_points_per_box)
    # The refined tree is also sorted. We use the refined tree keys as bins to sort the points accordingly.
    counts = sort_to_bins(sorted_point_keys, refined_tree)
    #We now know how many points are associated with each leaf key. We hence create a dictionary that maps from
    #leaf keys to the indices of the points associated with that leaf key.

    leaf_to_points = Dict{MortonKey,Vector{Int64}}()
    idx = 1
    for (count_index, count) in enumerate(counts)
        leaf_to_points[refined_tree[count_index]] = sorted_indices[idx:idx+count-1]
        idx += count
    end

    #We also create the inverted map of points to leaf keys.
    points_to_leafs = Vector{MortonKey}(undef, number_of_points)
    for (leaf_key, indices) in leaf_to_points
        for idx in indices
            points_to_leafs[idx] = leaf_key
        end
    end

    AdaptiveOctree(max_depth, number_of_points, points, point_morton_indices, refined_tree, leaf_to_points, points_to_leafs)

end

"""
    refine_tree(root::MortonKey, sorted_keys::Vector{MortonKey}, max_depth::Int64, max_points_per_box::Int64) -> Vector{MortonKey}

Refines a tree starting at `root` with the given `sorted_keys` into a tree with a maximum depth of `max_depth` and a maximum number of points per box of `max_points_per_box`.
"""
function refine_tree(root::MortonKey, sorted_keys::Vector{MortonKey}, max_depth::Int64, max_points_per_box::Int64)::Vector{MortonKey}

    if length(sorted_keys) <= max_points_per_box || Morton.level(root) >= max_depth
        return [root]
    end

    fine_keys = Vector{MortonKey}()
    children = Morton.children(root)
    counts = sort_to_bins(sorted_keys, children)

    idx = 1
    for (child_index, count) in enumerate(counts)
        append!(fine_keys, refine_tree(children[child_index], sorted_keys[idx:idx+count-1], max_depth, max_points_per_box))
        idx += count
    end

    fine_keys

end


"""
   sort_to_bins(sorted_keys::T, bins::V) -> Vector{Int64}

Sorts a list of sorted keys into n bins. We require that sorted_keys[1] >= bins[1].
The bins are assumed to be sorted in ascending order and defined by [bins[j], bins[j+1]) for j=1,.., n - 1, that is
the bins are half open intervals. The last bin is defined as the half-open interval [bins[n], ∞).

The function returns a vector of length n with the number of elements in each bin.
"""
function sort_to_bins(sorted_keys::T, bins::V)::Vector{Int64} where {S,T<:AbstractVector{S},V<:AbstractVector{S}}

    nbins = length(bins)

    @assert bins[1] <= sorted_keys[1]

    if nbins == 1
        return Vector{Int64}([length(sorted_keys)])
    end

    bin_counts = zeros(Int64, nbins)

    count = 0
    idx = 1
    break_outer = false

    for key in sorted_keys
        if break_outer
            # Break out of the outer loop.
            break
        end
        # The key falls into the current bin
        if bins[idx] <= key < bins[idx+1]
            bin_counts[idx] += 1
            count += 1
        else
            # If already in last bin simply break out
            if idx == nbins - 1
                break
            end
            # Move the bin forward until it fits.
            while idx < nbins - 1
                idx += 1
                if bins[idx] <= key < bins[idx+1]
                    # The bin fits. Proceed as normal and break out of the
                    # inner while loop.
                    bin_counts[idx] += 1
                    count += 1
                    break
                elseif idx == nbins - 1
                    # We have no more fitting bins. Set the corresponding flag.
                    # Need this so we can break out of the outer for loop
                    break_outer = true
                end
            end
        end
    end

    # Need to add the remaining elements to the last bin
    # that is half open to infinity.

    bin_counts[end] = length(sorted_keys) - count

    bin_counts
end


end