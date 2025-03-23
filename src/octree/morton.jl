
module MortonKey

import ..Octree.Constants

struct MKey
    value::UInt64
end

function ==(key1::MKey, key2::MKey)::Bool
    key1.value == key2.value
end

"The largest possible key."
function upper_bound()
    typemax(UInt64)
end

"Return an invalid key."
function invalid_key()
    MKey(Constants.HIGHEST_BIT_MASK)
end

"Check if a key is valid."
function is_valid(key::MKey)::Bool
    key.value >> 63 !== 1
end

"Return the root key."
function root()
    MKey(0)
end

"""
    from_index_and_level(x_index::UInt64, y_index::UInt64, z_index::UInt64, level::UInt64)

Create a key from an index and a level.
"""
function from_index_and_level(x_index::UInt64, y_index::UInt64, z_index::UInt64, level::UInt64)::MKey

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
    MKey(key | level)


end

"Return the level of a key."
function level(key::MKey)::UInt64
    key.value & Constants.LEVEL_MASK

end

"Decode a key into `(level, (x_index, y_index, z_index))`."
function decode(key::MKey)::Tuple{UInt64,Tuple{UInt64,UInt64,UInt64}}
    function decode_key_helper(key::UInt64, lookup_table::Vector{UInt64})::UInt64
        N_LOOPS::UInt64 = 6
        coord::UInt64 = 0

        for index in 0:N_LOOPS-1
            coord |= lookup_table[((key>>(index*9))&Constants.NINE_BIT_MASK)+1] << (3 * index)
        end

        coord
    end

    level = MortonKey.level(key)

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
function parent(key::MKey)::MKey
    level = MortonKey.level(key)

    @assert level > 0

    bit_displacement = Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - level)
    mask = ~(7 << bit_displacement)

    MKey((key.value & mask) - 1)
end

"Return true if `key` is an ancestor of or identical to `other_key`."
function is_ancestor(key::MKey, other_key::MKey)::Bool
    my_level = MortonKey.level(key)
    other_level = MortonKey.level(other_key)
    if ~MortonKey.is_valid(key) || ~MortonKey.is_valid(other_key)
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

function finest_common_ancestor(key::MKey, other_key::MKey)::Mkey

    if key.value == other.value
        return key
    end

    my_level = MortonKey.level(key)
    other_level = MortonKey.level(other_key)

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

    MKey(first_key | new_level)

end



end