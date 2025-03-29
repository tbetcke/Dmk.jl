using Test
using Dmk.Octree.Morton
using Dmk.Octree.Constants
using TestItemRunner

@testitem "test bin sorting" begin
    using Dmk.Octree

    sorted_keys = [1, 2, 3, 4, 5]
    bins = [1, 3, 4]
    counts = Octree.sort_to_bins(sorted_keys, bins)
    @test counts == [2, 1, 2]

    sorted_keys = [1, 2, 4, 5, 6, 19]
    bins = [1, 3, 4]

    counts = Octree.sort_to_bins(sorted_keys, bins)

    @test counts == [2, 0, 4]

    sorted_keys = [1]

    bins = [3, 5]

    @test_throws AssertionError Octree.sort_to_bins(sorted_keys, bins)

    sorted_keys = [1]
    bins = [1, 5, 6]
    counts = Octree.sort_to_bins(sorted_keys, bins)
    @test counts == [1, 0, 0]
end

@testitem "test refine tree" begin
    import Dmk.Octree
    import Dmk.Octree.Morton
    import Random
    using LinearAlgebra

    Random.seed!(1234)
    points = 0.9 .* rand(3, 1000) .- 0.5

    # Normalize the points so that they correspond to the unit sphere.
    for i in axes(points, 2)
        points[:, i] ./= norm(points[:, i])
        points[:, i] .*= 0.4
    end

    keys = sort(Octree.points_to_keys(points))
    tree = Octree.refine_tree(Morton.root(), keys, 4, 30)
    @test Morton.is_linear_complete_and_balanced(tree)

end

@testitem "test adaptive octree" begin

    import Dmk.Octree
    import Dmk.Tools
    import Dmk.Octree.Morton

    points = Tools.sphere_points(10000)

    octree = Octree.adaptive_octree(points, 4, 30)

    @test Morton.is_linear_complete_and_balanced(octree.leaf_keys)

    # Now check the neighbours of the keys

    test_failed = false

    for (key, key_data) in octree.all_keys
        if key == Morton.root()
            continue
        end
        for neighbour in key_data.neighbours
            # Check that the stored neighbours are actually in the tree

            @test haskey(octree.all_keys, neighbour)

            # If neighbour is on the same level we must be in the neighbour list of the neighbour
            if Morton.level(key) == Morton.level(neighbour)
                @test key in octree.all_keys[neighbour].neighbours
            else
                # Neighbour must be one level lower. It cannot be on a higher level
                # because then the neighbour would have an interior key as neighbour
                # on the same level and we would not be in the else branch.
                parent = Morton.parent(key)
                @test parent in octree.all_keys[neighbour].neighbours
                # In that case neighbour must also be a leaf key
                @test octree.all_keys[neighbour].key_type == Octree.leaf
            end

        end

        # Check for each Morton neighbour on same level that the neighbour or its parent are
        # in the tree.

        for morton_neighbour in filter((key) -> Morton.is_valid(key), Morton.neighbours(key))
            @test morton_neighbour in key_data.neighbours || Morton.parent(morton_neighbour) in key_data.neighbours
        end

    end


end

