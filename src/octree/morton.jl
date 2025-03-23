
module Morton

import ..Octree.Constants

struct MortonKey
    value::UInt64
end

Base.:(==)(a::MortonKey, b::MortonKey)::Bool = a.value == b.value

"The largest possible key."
function upper_bound()
    typemax(UInt64)
end

"Return an invalid key."
function invalid_key()
    MortonKey(Constants.HIGHEST_BIT_MASK)
end

"Check if a key is valid."
function is_valid(key::MortonKey)::Bool
    key.value >> 63 !== 1
end

"Return the root key."
function root()
    MortonKey(0)
end

"""
    from_index_and_level(x_index::UInt64, y_index::UInt64, z_index::UInt64, level::UInt64)

Create a key from an index and a level.
"""
function from_index_and_level(x_index::UInt64, y_index::UInt64, z_index::UInt64, level::UInt64)::MortonKey

    level_diff = Constants.DEEPEST_LEVEL - level

    #If we are not on the deepest level we need to shift the box.
    #The box with x-index one on DEEPEST_LEVEL-1 has index two on
    #DEEPEST_LEVEL.

    x = x_index << level_diff
    y = y_index << level_diff
    z = z_index << level_diff

    key = (Constants.X_LOOKUP_ENCODE[1+((x>>Constants.BYTE_DISPLACEMENT)&Constants.BYTE_MASK)]
           | Constants.Y_LOOKUP_ENCODE[1+((y>>Constants.BYTE_DISPLACEMENT)&Constants.BYTE_MASK)]
           | Constants.Z_LOOKUP_ENCODE[1+((z>>Constants.BYTE_DISPLACEMENT)&Constants.BYTE_MASK)])

    key = ((key << 24)
           | Constants.X_LOOKUP_ENCODE[1+(x&Constants.BYTE_MASK)]
           | Constants.Y_LOOKUP_ENCODE[1+(y&Constants.BYTE_MASK)]
           | Constants.Z_LOOKUP_ENCODE[1+(z&Constants.BYTE_MASK)])

    key = key << Constants.LEVEL_DISPLACEMENT
    MortonKey(key | level)


end

"Return the level of a key."
function level(key::MortonKey)::UInt64
    key.value & Constants.LEVEL_MASK

end

"Decode a key into `(level, (x_index, y_index, z_index))`."
function decode(key::MortonKey)::Tuple{UInt64,Tuple{UInt64,UInt64,UInt64}}
    function decode_key_helper(key::UInt64, lookup_table::Vector{UInt64})::UInt64
        N_LOOPS::UInt64 = 6
        coord::UInt64 = 0

        for index in 0:N_LOOPS-1
            coord |= lookup_table[((key>>(index*9))&Constants.NINE_BIT_MASK)+1] << (3 * index)
        end

        coord
    end

    level = Morton.level(key)

    level_diff = Constants.DEEPEST_LEVEL - level
    key = key.value >> Constants.LEVEL_DISPLACEMENT

    x = decode_key_helper(key, Constants.X_LOOKUP_DECODE)
    y = decode_key_helper(key, Constants.Y_LOOKUP_DECODE)
    z = decode_key_helper(key, Constants.Z_LOOKUP_DECODE)

    y = y >> level_diff
    z = z >> level_diff
    x = x >> level_diff

    (level, (x, y, z))

end

"Return the parent of a key."
function parent(key::MortonKey)::MortonKey
    level = Morton.level(key)

    @assert level > 0

    bit_displacement = Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - level)
    mask = ~(7 << bit_displacement)

    MortonKey((key.value & mask) - 1)
end

"Return true if `key` is an ancestor of or identical to `other_key`."
function is_ancestor(key::MortonKey, other_key::MortonKey)::Bool
    my_level = Morton.level(key)
    other_level = Morton.level(other_key)
    if ~Morton.is_valid(key) || ~Morton.is_valid(other_key)
        return false
    end

    if key == other_key
        return true
    elseif my_level > other_level
        return false
    else
        # We shift both keys out to 3 * DEEPEST_LEVEL - my_level
        # This gives identical bit sequences if my_key is an ancestor of other_key
        my_key = key.value >> (Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - my_level))
        other_key = other_key.value >> (Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - my_level))

        return my_key == other_key
    end

end

function finest_common_ancestor(key::MortonKey, other_key::MortonKey)::MortonKey

    if key == other_key
        return key
    end

    my_level = Morton.level(key)
    other_level = Morton.level(other_key)

    # We bring both keys to the minimum of the two levels.

    level = min(my_level, other_level)

    # Remove the level information and bring second key to the same level as first key
    # After the following operation the least significant bits are associated with `first_level`.

    first_key = key.value >> (Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - level))
    second_key = other_key.value >> (Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - level))

    # Now move both keys up until they are identical.
    # At the same time we reduce the first level.

    count = 0
    while first_key != second_key
        count += 1
        first_key >>= 3
        second_key >>= 3
    end

    # We now return the ancestor at the given level.

    new_level = level - count

    first_key <<= 3 * (Constants.DEEPEST_LEVEL - new_level) + Constants.LEVEL_DISPLACEMENT

    MortonKey(first_key | new_level)
end

"Return true if the key is root."
function is_root(key::MortonKey)::Bool
    key == 0
end

"Return the children of a key."
function children(key::MortonKey)::Tuple{MortonKey,MortonKey,MortonKey,MortonKey,MortonKey,MortonKey,MortonKey,MortonKey}

    level = Morton.level(key)

    @assert level < Constants.DEEPEST_LEVEL

    child_level = level + 1

    shift = Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - child_level)

    key = key.value

    (MortonKey(1 + (key | 0 << shift)),
        MortonKey(1 + (key | 1 << shift)),
        MortonKey(1 + (key | 2 << shift)),
        MortonKey(1 + (key | 3 << shift)),
        MortonKey(1 + (key | 4 << shift)),
        MortonKey(1 + (key | 5 << shift)),
        MortonKey(1 + (key | 6 << shift)),
        MortonKey(1 + (key | 7 << shift)))


end

end