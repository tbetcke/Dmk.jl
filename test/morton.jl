# Unit tests for Morton encoding and decoding

using Test
using Dmk.Octree.Morton
using Dmk.Octree.Constants
using TestItemRunner

@testitem "Root key" begin

    key = Dmk.Octree.Morton.root()
    @test Dmk.Octree.Morton.is_valid(key)
    @test Dmk.Octree.Morton.level(key) == 0
    @test Dmk.Octree.Morton.decode(key) == (0, (0, 0, 0))
end

@testitem "Encoding and decoding" begin

    key = Dmk.Octree.Morton.from_index_and_level(Dmk.Octree.Constants.LEVEL_SIZE - 1,
        Dmk.Octree.Constants.LEVEL_SIZE - 1,
        Dmk.Octree.Constants.LEVEL_SIZE - 1,
        Dmk.Octree.Constants.DEEPEST_LEVEL)
    @test Dmk.Octree.Morton.is_valid(key)
    @test Dmk.Octree.Morton.level(key) == Dmk.Octree.Constants.DEEPEST_LEVEL
    @test Dmk.Octree.Morton.decode(key) == (Dmk.Octree.Constants.DEEPEST_LEVEL,
        (Dmk.Octree.Constants.LEVEL_SIZE - 1,
            Dmk.Octree.Constants.LEVEL_SIZE - 1,
            Dmk.Octree.Constants.LEVEL_SIZE - 1))

end

@testitem "Decode Tables" begin
    using Dmk.Octree.Constants

    for (index, actual) in enumerate(Constants.Z_LOOKUP_DECODE)
        index = index - 1
        expected = index & 1
        expected |= ((index >> 3) & 1) << 1
        expected |= ((index >> 6) & 1) << 2

        @test actual == expected
    end

    for (index, actual) in enumerate(Constants.Y_LOOKUP_DECODE)
        index = index - 1
        expected = (index >> 1) & 1
        expected |= ((index >> 4) & 1) << 1
        expected |= ((index >> 7) & 1) << 2

        @test actual == expected
    end

    for (index, actual) in enumerate(Constants.X_LOOKUP_DECODE)
        index = index - 1
        expected = (index >> 2) & 1
        expected |= ((index >> 5) & 1) << 1
        expected |= ((index >> 8) & 1) << 2

        @test actual == expected
    end


end

@testitem "Encode tables" begin
    using Dmk.Octree.Constants

    for (index, actual) in enumerate(Constants.Z_LOOKUP_ENCODE)
        sum = 0

        index = index - 1

        for shift in 0:7
            sum |= ((index & 1) << (3 * shift))
            index >>= 1
        end

        @test sum == actual
    end

    for (index, actual) in enumerate(Constants.Y_LOOKUP_ENCODE)
        sum = 0

        index = index - 1

        for shift in 0:7
            sum |= ((index & 1) << (3 * shift + 1))
            index >>= 1
        end

        @test sum == actual
    end

    for (index, actual) in enumerate(Constants.X_LOOKUP_ENCODE)
        sum = 0

        index = index - 1

        for shift in 0:7
            sum |= ((index & 1) << (3 * shift + 2))
            index >>= 1
        end

        @test sum == actual
    end

end

@testitem "parent" begin

    index = (15, 39, 45)
    key = Dmk.Octree.Morton.from_index_and_level(index..., 9)
    parent = Dmk.Octree.Morton.parent(key)

    expected_index = (7, 19, 22)

    (actual_level, actual_index) = Dmk.Octree.Morton.decode(parent)

    @test actual_level == 8
    @test actual_index == expected_index

end

@testitem "ancestor" begin
    import Dmk.Octree.Morton: parent

    index = (15, 39, 45)
    key = Dmk.Octree.Morton.from_index_and_level(index..., 9)
    @test Dmk.Octree.Morton.is_ancestor(key, key)

    ancestor = parent(parent(key))
    @test Dmk.Octree.Morton.is_ancestor(ancestor, key)
end

@testitem "finest common ancestor" begin

    import Dmk.Octree.Morton: parent
    import Dmk.Octree.Constants

    index = (15, 39, 45)
    key = Dmk.Octree.Morton.from_index_and_level(index..., 9)
    ancestor = parent(parent(key))

    @test Dmk.Octree.Morton.finest_common_ancestor(key, ancestor) == ancestor

    # The finest ancestor of the following keys should be the root of the tree.

    key1 = Dmk.Octree.Morton.from_index_and_level(0, 0, 0, Constants.DEEPEST_LEVEL - 1)
    key2 = Dmk.Octree.Morton.from_index_and_level(Constants.LEVEL_SIZE - 1, Constants.LEVEL_SIZE - 1, Constants.LEVEL_SIZE - 1, Constants.DEEPEST_LEVEL)

    @test Dmk.Octree.Morton.finest_common_ancestor(key1, key2) == Dmk.Octree.Morton.root()
end

@testitem "children of a key" begin

    index = (15, 39, 45)
    key = Dmk.Octree.Morton.from_index_and_level(index..., 9)
    children = Dmk.Octree.Morton.children(key)

    @assert length(children) == 8

    for child in children
        @test Dmk.Octree.Morton.parent(child) == key
    end

end

@testitem "siblings of a key" begin

    index = (15, 39, 45)
    key = Dmk.Octree.Morton.from_index_and_level(index..., 9)
    siblings = Dmk.Octree.Morton.siblings(key)

    @assert length(siblings) == 8

    for sibling in siblings
        @test Dmk.Octree.Morton.parent(sibling) == Dmk.Octree.Morton.parent(key)
    end
end

@testitem "neighbours of a key" begin

    # Check neighbours of root
    key = Dmk.Octree.Morton.root()
    neighbours = Dmk.Octree.Morton.neighbours(key)

    # Take a single key in the interior and check all its neighbours

    key = Dmk.Octree.Morton.from_index_and_level(15, 39, 45, 9)
    neighbours = Dmk.Octree.Morton.neighbours(key)

    for neighbour in neighbours
        (n_level, (n_x, n_y, n_z)) = Dmk.Octree.Morton.decode(neighbour)
        @test abs(n_x - 15) <= 1
        @test abs(n_y - 39) <= 1
        @test abs(n_z - 45) <= 1
        @test n_level == 9
        @test neighbour != key
    end

end

@testitem "child index" begin

    key = Dmk.Octree.Morton.from_index_and_level(15, 39, 45, 9)
    children = Dmk.Octree.Morton.children(key)

    for (index, child) in enumerate(children)
        @test Dmk.Octree.Morton.child_index(child) + 1 == index
    end

end

@testitem "finest outer descendent" begin
    using Dmk.Octree.Morton
    using Dmk.Octree.Constants

    key = Morton.from_index_and_level(0, 0, 0, 1)
    finest_outer_descendent = Morton.finest_outer_descendent(key)

    @test finest_outer_descendent == Morton.from_index_and_level(0, 0, 0, Constants.DEEPEST_LEVEL)

    key = Morton.from_index_and_level(Constants.LEVEL_SIZE - 1,
        Constants.LEVEL_SIZE - 1,
        Constants.LEVEL_SIZE - 1,
        Constants.DEEPEST_LEVEL)

    @test Morton.finest_outer_descendent(key) == key

    @test Morton.finest_outer_descendent(Morton.parent(Morton.parent(Morton.parent(key)))) == key


end

@testitem "next deepest nondescendent key" begin
    import Dmk.Octree.Morton
    import Dmk.Octree.Constants

    key = Morton.from_index_and_level(1, 1, 1, 1)
    @test Morton.next_deepest_nondescendent_key(key) == Morton.invalid_key()

    key = Morton.parent(Morton.parent(Morton.parent(Morton.deepest_last())))
    Morton.next_deepest_nondescendent_key(Morton.children(key)[7]) == Morton.deepest_last()

end


@testitem "fill between keys" begin
    import Dmk.Octree.Morton: MortonKey, fill_between_keys, is_complete
    import Dmk.Octree.Morton

    function sanity_check(key1::MortonKey, key2::MortonKey, keys::Vector{MortonKey})
        max_level = max(Morton.level(key1), Morton.level(key2))

        keys = [key1; keys; key2]

        for (k1, k2) in zip(keys[1:end-1], keys[2:end])
            @test k1 < k2
            @test !Morton.is_ancestor(k1, k2)
        end

        for k in keys
            @test Morton.level(k) <= max_level
        end
    end

    key1 = Morton.from_index_and_level(0, 1, 0, 4)
    key2 = Morton.from_index_and_level(8, 4, 13, 4)
    keys = fill_between_keys(key1, key2)

    @test !isempty(keys)

    sanity_check(key1, key2, keys)

    # Correct result for passing same key twice

    keys = Morton.fill_between_keys(key1, key1)
    @test isempty(keys)

    # Two consecutive keys should also be empty.

    children = Morton.children(key2)
    keys = Morton.fill_between_keys(children[1], children[2])
    @test isempty(keys)


    # The following should have deepest level same as the level of key2.

    key1 = Morton.children(Morton.root())[1]
    key2 = Morton.children(Morton.children(Morton.root())[8])[8]

    keys = Morton.fill_between_keys(key1, key2)
    max_level = maximum((key) -> Morton.level(key), keys)

    @test max_level == 2

    # Also check that the resulting keys are complete.

    @test is_complete([key1; keys; key2])


end

@testitem "test complete tree" begin
    import Dmk.Octree.Morton
    import Dmk.Octree.Morton: MortonKey

    function sanity_checks(keys::Vector{MortonKey}, complete_region::Vector{MortonKey})

        # Check that the returned region is linear
        @test Morton.is_linear(complete_region)

        # Check that the first key of the region is an ancestor of the first key of the deepest level
        # and that the last key of the region is an ancestor of the last key of the deepest level.

        deepest_first = Morton.deepest_first()
        deepest_last = Morton.deepest_last()

        @test Morton.is_ancestor(complete_region[1], deepest_first)
        @test Morton.is_ancestor(complete_region[end], deepest_last)

        # Check that the original keys are all contained in the complete region.

        for key in keys
            @test key in complete_region
        end

        # Check that the maximum level of the complete region is not higher than the maximum level of the input keys.

        max_level = maximum((key) -> Morton.level(key), keys)
        @test maximum((key) -> Morton.level(key), complete_region) <= max_level
    end

    # Create 3 Morton keys around which to complete region. The keys are not sorted.

    key1 = Morton.from_index_and_level(17, 30, 55, 10)
    key2 = Morton.from_index_and_level(17, 540, 55, 10)
    key3 = Morton.from_index_and_level(17, 30, 799, 11)

    keys = [key1, key2, key3]

    complete_region = Morton.complete_region(keys)

    sanity_checks(keys, complete_region)

    # For an empty slice the complete region method should just add the root of the tree.
    keys = Vector{MortonKey}()
    complete_region = Morton.complete_region(keys)
    @test length(complete_region) == 1
    @test complete_region[1] == Morton.root()

    # Choose a region where the first and last key are ancestors of deepest first and deepest last.

    keys = [Morton.deepest_first(), Morton.deepest_last()]

    complete_region = Morton.complete_region(keys)

    sanity_checks(keys, complete_region)


end

@testitem "test balancing" begin
    import Dmk.Octree.Morton

    balanced = Morton.balance([Morton.from_index_and_level(0, 1, 0, 2)])

    @test Morton.is_linear_complete_and_balanced(balanced)

    key1 = Morton.from_index_and_level(17, 35, 48, 9)
    key2 = Morton.from_index_and_level(355, 25, 67, 9)
    key3 = Morton.from_index_and_level(0, 0, 0, 8)

    balanced = Morton.balance([key1, key2, key3])

    @test Morton.is_linear_complete_and_balanced(balanced)

    # We start with all keys on level 1. We recurse twice down and replace
    # the first key by its descendents two levels down and linearize. The
    # resulting octree is complete and linear but not balanced.

    keys = Vector{Morton.MortonKey}(Morton.children(Morton.root()))
    descendents = Morton.children(keys[1])

    for key in descendents
        append!(keys, Morton.children(key))
    end

    keys = Morton.linearize(keys)

    # This tree is complete and linear but not balanced.

    @test Morton.is_linear(keys)
    @test Morton.is_complete(keys)
    @test !Morton.is_linear_complete_and_balanced(keys)

    # We now balance it

    keys = Morton.balance(keys)

    # Check again

    @test Morton.is_linear_complete_and_balanced(keys)

end

