# kd tree for fast searching of a list of objects
export KDTree

"""
A K-D tree of objects in `K` dimensions.

This is a binary tree.  Each node in the tree divides space along
coordinate `ix` into objects that overlap the region `≤ x` and `≥ x`.
The leaves of the tree are just a list of objects `o`.

This implementation is a little different from a standard K-D tree,
because a given object may be in both branches of the tree.  (Normally,
K-D trees are used for nearest-neighbor searches for a list of *points*,
not objects of nonzero size.)
"""
type KDTree{K}
    o::Vector{Object{K}}
    ix::Int
    x::Float64
    left::KDTree  # objects ≤ x in coordinate ix
    right::KDTree # objects > x in coordinate ix
    KDTree(o::AbstractVector{Object{K}}) = new(o, 0)
    function KDTree(x::Real, ix::Integer, left::KDTree{K}, right::KDTree{K})
        1 ≤ ix ≤ K || throw(BoundsError())
        new(Object{K}[], ix, x, left, right)
    end
end

Base.ndims{K}(::KDTree{K}) = K

"""
    KDTree(objects::AbstractVector{Object{K}})

Construct a K-D tree (`KDTree`) representation of a list of
`objects` in order to enable rapid searching of an object list.

When searching the tree, objects that appear earlier in `objects`
take precedence over objects that appear later.
"""
function KDTree{K}(o::AbstractVector{Object{K}})
    (length(o) <= 4 || K == 0) && return KDTree{K}(o)

    # figure out the best dimension ix to divide over,
    # the dividing plane x, and the number (nl,nr) of
    # objects that fall into the left and right subtrees
    b = map(bounds, o) # cartesian bounding boxes of all objects
    ix = 0
    x = 0.0
    nl = nr = typemax(Int)
    for i = 1:K
        mx = median(map(b -> 0.5*(b[1][i] + b[2][i]), b))
        mnl = count(b -> b[1][i] ≤ mx, b) # lower bound ≤ mx
        mnr = count(b -> b[2][i] > mx, b) # upper bound > mx
        if max(mnl,mnr) < max(nl,nr)
            ix = i
            x = mx
            nl = mnl
            nr = mnr
        end
    end

    # don't bother subdividing if it doesn't reduce the # of objects much
    4*min(nl,nr) > 3*length(o) && return KDTree{K}(o)

    # create the arrays of objects in each subtree
    ol = Array{Object{K}}(nl)
    or = Array{Object{K}}(nr)
    il = ir = 0
    for k in eachindex(o)
        if b[k][1][ix] ≤ x
            ol[il += 1] = o[k]
        end
        if b[k][2][ix] > x
            or[ir += 1] = o[k]
        end
    end

    return KDTree{K}(x, ix, KDTree(ol), KDTree(or))
end

depth(kd::KDTree) = kd.ix == 0 ? 0 : max(depth(kd.left), depth(kd.right)) + 1

Base.show{K}(io::IO, kd::KDTree{K}) = print(io, "KDTree{$K} of depth ", depth(kd))

function _show(io::IO, kd::KDTree, indent)
    indentstr = " "^indent
    if kd.ix == 0
        println(io, indentstr, length(kd.o), " objects")
    else
        println(io, indentstr, "if x[", kd.ix, "] ≤ ", kd.x, ':')
        _show(io, kd.left, indent + 2)
        println(io, indentstr, "else:")
        _show(io, kd.right, indent + 2)
    end
end

function Base.show{K}(io::IO, ::MIME"text/plain", kd::KDTree{K})
    println(io, kd, ':')
    _show(io, kd, 0)
end

function Base.findin{N}(p::SVector{N}, o::Vector{Object{N}})
    for i in eachindex(o)
        if p ∈ o[i]
            return Nullable{Object{N}}(o[i])
        end
    end
    return Nullable{Object{N}}()
end

function Base.findin{N}(p::SVector{N}, kd::KDTree{N})
    if isempty(kd.o)
        if p[kd.ix] ≤ kd.x
            return findin(p, kd.left)
        else
            return findin(p, kd.right)
        end
    else
        return findin(p, kd.o)
    end
end

"""
    findin(p::AbstractVector, kd::KDTree)

Return a `Nullable` container of the first object in `kd` that contains
the point `p`, or an `isnull` container if no object was found.
"""
Base.findin{N}(p::AbstractVector, kd::KDTree{N}) = findin(SVector{N}(p), kd)
Base.findin{N}(p::AbstractVector, o::Vector{Object{N}}) = findin(SVector{N}(p), o)
