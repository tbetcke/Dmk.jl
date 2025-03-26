
module Morton

import ..Octree.Constants

import StaticArrays: SVector, @SVector

struct MortonKey
    value::Int64
end

Base.:(==)(a::MortonKey, b::MortonKey)::Bool = a.value == b.value
Base.:copy(a::MortonKey)::MortonKey = MortonKey(a.value)
Base.:(<)(a::MortonKey, b::MortonKey)::Bool = a.value < b.value
Base.:(<=)(a::MortonKey, b::MortonKey)::Bool = a.value <= b.value
Base.:(>)(a::MortonKey, b::MortonKey)::Bool = a.value > b.value
Base.:(>=)(a::MortonKey, b::MortonKey)::Bool = a.value >= b.value

function Base.isless(a::MortonKey, b::MortonKey)
    a.value < b.value
end

"The largest possible key."
function upper_bound()
    typemax(Int64)
end

"Return an invalid key."
function invalid_key()
    MortonKey(-1)
end

"Check if a key is valid."
function is_valid(key::MortonKey)::Bool
    key.value >= 0
end

"Return the root key."
function root()
    MortonKey(0)
end

"""Return the deepest last key."""
function deepest_last()
    Morton.from_index_and_level(Constants.LEVEL_SIZE - 1, Constants.LEVEL_SIZE - 1, Constants.LEVEL_SIZE - 1, Constants.DEEPEST_LEVEL)
end

"""Return the deepest first key."""
function deepest_first()
    MortonKey(Constants.DEEPEST_LEVEL)
end

"""
    from_index_and_level(x_index::Int64, y_index::Int64, z_index::Int64, level::Int64)

Create a key from an index and a level.
"""
function from_index_and_level(x_index::Int64, y_index::Int64, z_index::Int64, level::Int64)::MortonKey

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
function level(key::MortonKey)::Int64
    key.value & Constants.LEVEL_MASK

end

"Decode a key into `(level, (x_index, y_index, z_index))`."
function decode(key::MortonKey)::Tuple{Int64,Tuple{Int64,Int64,Int64}}
    function decode_key_helper(key::Int64, lookup_table::Vector{Int64})::Int64
        N_LOOPS::Int64 = 6
        coord::Int64 = 0

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
    if !Morton.is_valid(key) || !Morton.is_valid(other_key)
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
function children(key::MortonKey)::SVector{8,MortonKey}

    level = Morton.level(key)

    @assert level < Constants.DEEPEST_LEVEL

    child_level = level + 1

    shift = Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - child_level)

    key = key.value

    SVector(MortonKey(1 + (key | 0 << shift)),
        MortonKey(1 + (key | 1 << shift)),
        MortonKey(1 + (key | 2 << shift)),
        MortonKey(1 + (key | 3 << shift)),
        MortonKey(1 + (key | 4 << shift)),
        MortonKey(1 + (key | 5 << shift)),
        MortonKey(1 + (key | 6 << shift)),
        MortonKey(1 + (key | 7 << shift)))

end

"Return the siblings of a key"
function siblings(key::MortonKey)::SVector{8,MortonKey}

    @assert !Morton.is_root(key)
    Morton.children(Morton.parent(key))
end

"Return the neighbours of a key. If a direction has no neighbour an invalid key is returned for that direction."
function neighbours(key::MortonKey)::SVector{26,MortonKey}

    (level, (x_index, y_index, z_index)) = Morton.decode(key)

    if level == 0
        return @SVector [Morton.invalid_key() for i in 1:26]
    end

    level_size = 1 << level

    function compute_neighbour(x_index::Int64, y_index::Int64, z_index::Int64)::MortonKey

        if x_index < 0 || x_index >= level_size || y_index < 0 || y_index >= level_size || z_index < 0 || z_index >= level_size
            return Morton.invalid_key()
        end

        Morton.from_index_and_level(x_index, y_index, z_index, level)

    end

    return @SVector [compute_neighbour(x_index + d[1], y_index + d[2], z_index + d[3])
                     for d in Constants.DIRECTIONS]


end

"Return the index of a child (between 0 and 7)"
function child_index(key::MortonKey)::Int64
    if key == Morton.root()
        return 0
    end

    level = Morton.level(key)

    shift = Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - level)

    ((key.value >> shift) % 8)
end

"Return the finest descendent that is opposite the joint corner with the siblings."
function finest_outer_descendent(key::MortonKey)::MortonKey

    if Morton.is_root(key)
        return Morton.from_index_and_level(0, 0, 0, Constants.DEEPEST_LEVEL)
    end

    # Get the index of the current key as a child.

    outer_index = Morton.child_index(key)
    child_level = 1 + Morton.level(key)

    new_key = copy(key)

    while child_level <= Constants.DEEPEST_LEVEL
        shift = Constants.LEVEL_DISPLACEMENT + 3 * (Constants.DEEPEST_LEVEL - child_level)
        new_key = MortonKey(1 + (new_key.value | outer_index << shift))
        child_level += 1
    end

    new_key



end

"""Return the next key that is not a descendent of the current key."""
function next_deepest_nondescendent_key(key::MortonKey)::MortonKey

    if Morton.is_ancestor(key, Morton.deepest_last())
        return Morton.invalid_key()
    end

    level = Morton.level(key)

    level_diff = Constants.DEEPEST_LEVEL - level
    shift = Constants.LEVEL_DISPLACEMENT + 3 * level_diff

    child_index = (key.value >> shift) % 8

    if child_index < 7
        # Take the next sibling and go to the deepest level
        MortonKey(key.value + (1 << shift) + level_diff)
    else
        # If we are the last sibling we go to the parent
        # and take the next deepest nondescendent key from there
        Morton.next_deepest_nondescendent_key(Morton.parent(key))
    end

end

"Get interior keys"
function interior_keys(keys::T)::Set{MortonKey} where {T<:AbstractVector{MortonKey}}

    result = Set{MortonKey}()

    for key in keys
        if Morton.is_root(key)
            continue
        end

        p = Morton.parent(key)
        while Morton.level(p) > 0 && !(p in result)
            push!(result, p)
            p = Morton.parent(p)
        end
    end

    push!(result, Morton.root())

    result
end

"Get all keys"
function all_keys(keys::T)::Set{MortonKey} where {T<:AbstractVector{MortonKey}}

    result = Morton.interior_keys(keys)
    union!(result, keys)

    result
end



"Linearize by removing duplicates and overlaps."
function linearize(keys::T)::Vector{MortonKey} where {T<:AbstractVector{MortonKey}}

    new_keys = Vector{MortonKey}()

    if !isempty(keys)
        sorted_keys = sort(keys)
        for (indx1, indx2) in zip(1:length(sorted_keys)-1, 2:length(sorted_keys))
            key1 = sorted_keys[indx1]
            key2 = sorted_keys[indx2]

            if key1 == key2 || Morton.is_ancestor(key1, key2)
                continue
            end
            push!(new_keys, key1)
        end
        push!(new_keys, sorted_keys[end])
        new_keys
    end

    new_keys
end

"Fill the region between two keys with a minimal number of keys."
function fill_between_keys(key1::MortonKey, key2::MortonKey)::Vector{MortonKey}

    # Make sure that key1 is smaller or equal key2

    if key2 > key1
        key1, key2 = key1, key2
    end

    # If key1 is ancestor of key2 return empty list. Note that
    # is_ancestor is true if key1 is equal to key2.

    if Morton.is_ancestor(key1, key2)
        return Vector{MortonKey}()
    end

    # The finest common ancestor is always closer to the root than either key
    # if key1 is not an ancestor of key2 or vice versa.

    ancestor = Morton.finest_common_ancestor(key1, key2)
    children = Morton.children(ancestor)

    result = Vector{MortonKey}()

    work_set = Set(children)

    while !isempty(work_set)
        item = pop!(work_set)

        # If the item is either key we don't want it in the result.
        if item == key1 || item == key2
            continue
            # We want items that are strictly between the two keys and are not ancestors of either.
            # We do not check specifically if item is an ancestor of key1 as then it would be smaller than key1.
        elseif key1 < item && item < key2 && !Morton.is_ancestor(item, key2)
            push!(result, item)
        else
            # If the item is an ancestor of key1 or key2 just refine to the children and try again.
            # Note we already exclude that item is identical to key1 or key2.
            # So if item is an ancestor of either its children cannot have a level larger than key1 or key2.
            if Morton.is_ancestor(item, key1) || Morton.is_ancestor(item, key2)
                children = Morton.children(item)
                push!(work_set, children...)
            end
        end
    end

    sort!(result)

    result

end

"Complete a region ensuring that the given keys are part of the leafs."
function complete_region(keys::T)::Vector{MortonKey} where {T<:AbstractVector{MortonKey}}

    # First make sure that the input sequence is sorted.
    keys = sort(keys)

    result = Vector{MortonKey}()

    # Special case of empty keys.
    if isempty(keys)
        push!(result, Morton.from_index_and_level(0, 0, 0, 0))
        return result
    end

    # If just the root is given return that.
    if length(keys) == 1 && keys[1] == Morton.root()
        return copy(keys)
    end

    deepest_first = Morton.deepest_first()
    deepest_last = Morton.deepest_last()

    # If the first key is not an ancestor of the deepest possible first element in the
    # tree get the finest ancestor between the two and use the first child of that.

    first_key = keys[1]
    last_key = keys[end]

    if !Morton.is_ancestor(first_key, deepest_first)
        ancestor = Morton.finest_common_ancestor(deepest_first, first_key)
        pushfirst!(keys, Morton.children(ancestor)[1])
    end

    if !Morton.is_ancestor(last_key, deepest_last)
        ancestor = Morton.finest_common_ancestor(deepest_last, last_key)
        push!(keys, Morton.children(ancestor)[Constants.NSIBLINGS])
    end

    # Now just iterate over the keys by tuples of two and fill the region between two keys.

    for (key1, key2) in zip(keys[1:end-1], keys[2:end])
        push!(result, key1)
        append!(result, Morton.fill_between_keys(key1, key2))
    end

    # Push the final key
    push!(result, keys[end])
    # We do not sort the keys. They are already sorted.
    result
end

"2:1 balance a tree given by a list of leaf keys. The result is not only balanced but also complete."
function balance(keys::T)::Vector{MortonKey} where {T<:AbstractVector{MortonKey}}

    if Morton.isempty(keys)
        return Vector{MortonKey}()
    end

    deepest_level = maximum((key) -> Morton.level(key), keys)

    if deepest_level == 0
        return [Morton.root()]
    end

    # Start with keys at deepest level

    work_list = filter((key) -> Morton.level(key) == deepest_level, keys)

    result = Vector{MortonKey}()

    for level in deepest_level:-1:1
        parents = Set{MortonKey}()
        new_work_list = Vector{MortonKey}()

        for key in work_list
            parent = Morton.parent(key)
            if !(parent in parents)
                push!(parents, parent)
                append!(result, Morton.siblings(key))
                append!(new_work_list, filter((key) -> Morton.is_valid(key), Morton.neighbours(parent)))
            end
        end

        append!(new_work_list, filter((key) -> Morton.level(key) == level - 1, keys))

        work_list = new_work_list
    end

    Morton.linearize(result)
end

"Return true if a vector of keys is linear."
function is_linear(keys::T)::Bool where {T<:AbstractVector{MortonKey}}

    if isempty(keys)
        return true
    end

    for (key1, key2) in zip(keys[1:end-1], keys[2:end])
        if key1 >= key2 || Morton.is_ancestor(key1, key2)
            return false
        end
    end

    true

end

"Return true if a vector of keys is complete."
function is_complete(keys::T)::Bool where {T<:AbstractVector{MortonKey}}

    if isempty(keys)
        return false
    end

    # Get the interior keys

    interior_keys = Morton.interior_keys(keys)

    # Get all keys

    all_keys = Morton.all_keys(keys)

    # Make sure that all 8 children of each interior key are in the set of all keys.

    for key in interior_keys
        children = Morton.children(key)
        for child in children
            if !(child in all_keys)
                return false
            end
        end
    end

    true

end

"Return true if a vector of keys is linear, complete and balanced."
function is_linear_complete_and_balanced(keys::T)::Bool where {T<:AbstractVector{MortonKey}}

    if isempty(keys)
        return false
    end

    # We add for each key the neighbours of the parents and linearize. If the key was balanced
    # this will give the same result as the input.

    new_keys = Vector{MortonKey}()

    for key in keys
        push!(new_keys, key, filter((key) -> Morton.is_valid(key), Morton.neighbours(Morton.parent(key)))...)
    end

    new_keys = Morton.linearize(new_keys)

    for (key1, key2) in zip(keys, new_keys)
        if key1 != key2
            return false
        end
    end

    true


end
end