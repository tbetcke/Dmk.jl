# Unit tests for Morton encoding and decoding

using Test
using Dmk.Octree.MortonKey
using Dmk.Octree.Constants
using TestItemRunner

@testitem "Root key" begin

    key = Dmk.Octree.MortonKey.root()
    @test Dmk.Octree.MortonKey.is_valid(key)
    @test Dmk.Octree.MortonKey.level(key) == 0
    @test Dmk.Octree.MortonKey.decode(key) == (0, (0, 0, 0))
end

@testitem "Encoding and decoding" begin

    key = Dmk.Octree.MortonKey.from_index_and_level(Dmk.Octree.Constants.LEVEL_SIZE - 1,
        Dmk.Octree.Constants.LEVEL_SIZE - 1,
        Dmk.Octree.Constants.LEVEL_SIZE - 1,
        Dmk.Octree.Constants.DEEPEST_LEVEL)
    @test Dmk.Octree.MortonKey.is_valid(key)
    @test Dmk.Octree.MortonKey.level(key) == Dmk.Octree.Constants.DEEPEST_LEVEL
    @test Dmk.Octree.MortonKey.decode(key) == (Dmk.Octree.Constants.DEEPEST_LEVEL,
        (Dmk.Octree.Constants.LEVEL_SIZE - 1,
            Dmk.Octree.Constants.LEVEL_SIZE - 1,
            Dmk.Octree.Constants.LEVEL_SIZE - 1))

end

@testitem "Decode Tables" begin
    using Dmk.Octree.Constants

    for (index::UInt64, actual) in enumerate(Constants.Z_LOOKUP_DECODE)
        index = index - 1
        expected = index & 1
        expected |= ((index >> 3) & 1) << 1
        expected |= ((index >> 6) & 1) << 2

        @test actual == expected
    end

    for (index::UInt64, actual) in enumerate(Constants.Y_LOOKUP_DECODE)
        index = index - 1
        expected = (index >> 1) & 1
        expected |= ((index >> 4) & 1) << 1
        expected |= ((index >> 7) & 1) << 2

        @test actual == expected
    end

    for (index::UInt64, actual) in enumerate(Constants.X_LOOKUP_DECODE)
        index = index - 1
        expected = (index >> 2) & 1
        expected |= ((index >> 5) & 1) << 1
        expected |= ((index >> 8) & 1) << 2

        @test actual == expected
    end


end

@testitem "Encode tables" begin
    using Dmk.Octree.Constants

    for (index::UInt64, actual) in enumerate(Constants.Z_LOOKUP_ENCODE)
        sum::UInt64 = 0

        index = index - 1

        for shift in 0:7
            sum |= ((index & 1) << (3 * shift))
            index >>= 1
        end

        @test sum == actual
    end

    for (index::UInt64, actual) in enumerate(Constants.Y_LOOKUP_ENCODE)
        sum::UInt64 = 0

        index = index - 1

        for shift in 0:7
            sum |= ((index & 1) << (3 * shift + 1))
            index >>= 1
        end

        @test sum == actual
    end

    for (index::UInt64, actual) in enumerate(Constants.X_LOOKUP_ENCODE)
        sum::UInt64 = 0

        index = index - 1

        for shift in 0:7
            sum |= ((index & 1) << (3 * shift + 2))
            index >>= 1
        end

        @test sum == actual
    end

end

@testitem "parent" begin

    index::Tuple{UInt64,UInt64,UInt64} = (15, 39, 45)
    key = Dmk.Octree.MortonKey.from_index_and_level(index..., UInt64(9))
    parent = Dmk.Octree.MortonKey.parent(key)

    expected_index = (7, 19, 22)

    (actual_level, actual_index) = Dmk.Octree.MortonKey.decode(parent)

    @test actual_level == 8
    @test actual_index == expected_index

end

@testitem "ancestor" begin
    import Dmk.Octree.MortonKey: parent

    index::Tuple{UInt64,UInt64,UInt64} = (15, 39, 45)
    key = Dmk.Octree.MortonKey.from_index_and_level(index..., UInt64(9))
    @test Dmk.Octree.MortonKey.is_ancestor(key, key)

    ancestor = parent(parent(key))
    @test Dmk.Octree.MortonKey.is_ancestor(ancestor, key)
end




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