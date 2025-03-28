using Test
using Dmk.Octree.Morton
using Dmk.Octree.Constants
using TestItemRunner

@testitem "test bin sorting" begin
    using Dmk.Octree

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