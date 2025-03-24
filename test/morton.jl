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



# let index = [15, 39, 45];

#     let key = Morton::from_index_and_level(index, 9);
#     // The finest ancestor with itself is the key itself.
#     assert_eq!(key.finest_common_ancestor(key), key);
#     // Finest ancestor with ancestor two levels up is the ancestor.
#     let ancestor = key.parent().parent();
#     assert_eq!(key.finest_common_ancestor(ancestor), ancestor);

#     // Finest ancestor  of the following keys should be the root of the tree.

#     let key1 = Morton::from_index_and_level([0, 0, 0], DEEPEST_LEVEL as usize - 1);
#     let key2 = Morton::from_index_and_level(
#         [
#             LEVEL_SIZE as usize - 1,
#             LEVEL_SIZE as usize - 1,
#             LEVEL_SIZE as usize - 1,
#         ],
#         DEEPEST_LEVEL as usize,
#     );

#     assert_eq!(
#         key1.finest_common_ancestor(key2),
#         Morton::from_index_and_level([0, 0, 0], 0)
#     );

#     // The finest ancestor of these two keys should be at level 1.

#     let key1 = Morton::from_index_and_level([0, 0, 62], 6);
#     let key2 = Morton::from_index_and_level([0, 0, 63], 6);
#     let expected = Morton::from_index_and_level([0, 0, 31], 5);

#     assert_eq!(key1.finest_common_ancestor(key2), expected);


# fn test_z_encode_table() {
#     for (mut index, actual) in Z_LOOKUP_ENCODE.iter().enumerate() {
#         let mut sum: u64 = 0;

#         for shift in 0..8 {
#             sum |= ((index & 1) << (3 * shift)) as u64;
#             index >>= 1;
#         }

#         assert_eq!(sum, *actual);
#     }
# }

# fn test_z_decode_table() {
#     for (index, &actual) in Z_LOOKUP_DECODE.iter().enumerate() {
#         let mut expected: u64 = (index & 1) as u64;
#         expected |= (((index >> 3) & 1) << 1) as u64;
#         expected |= (((index >> 6) & 1) << 2) as u64;

#         assert_eq!(actual, expected);
#     }
# }

# @testset "Morton encoding and decoding" begin
#     @testset "Morton encoding" begin
#         @test "Root key" begin
#             key = root()
#             @test is_valid(key)
#             @test level(key) == 0
#             @test decode(key) == (0, (0, 0, 0))
#         end

#         @test "Invalid key" begin
#             key = invalid_key()
#             @test !is_valid(key)
#         end

#         @test "Key from index and level" begin
#             key = from_index_and_level(1, 1, 1, 0)
#             @test is_valid(key)
#             @test level(key) == 0
#             @test decode(key) == (0, (0, 0, 0))

#             key = from_index_and_level(1, 1, 1, 1)
#             @test is_valid(key)
#             @test level(key) == 1
#             @test decode(key) == (1, (0, 0, 0))

#             key = from_index_and_level(2, 2, 2, 1)
#             @test is_valid(key)
#             @test level(key) == 1
#             @test decode(key) == (73, (1, 1, 1))

#             key = from_index_and_level(2, 2, 2, 2)
#             @test is_valid(key)
#             @test level(key) == 2
#             @test decode(key) == (292, (2, 2, 2))

#             key = from_index_and_level(4, 4, 4, 3)
#             @test is_valid(key)
#             @test level(key) == 3
#             @test decode(key) == (1176, (4, 4, 4))

#             key = from_index_and_level(8, 8, 8, 4)
#             @test is_valid(key)
#             @test level(key) == 4
#             @test decode(key) == (4704, (8, 8, 8))

#             key = from_index_and_level(16, 16, 16, 5)
#             @test is_valid(key)
#             @test level(key) == 5
#             @test decode(key) == (18816, (16, 16,