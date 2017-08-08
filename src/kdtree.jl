# kd tree for fast searching of a list of shapes
export KDTree

"""
A K-D tree of shapes in `K` dimensions.

This is a binary tree.  Each node in the tree divides space along
coordinate `ix` into shapes that overlap the region `≤ x` and `≥ x`.
The leaves of the tree are just a list of shapes `s`.

This implementation is a little different from a standard K-D tree,
because a given shape may be in both branches of the tree.  (Normally,
K-D trees are used for nearest-neighbor searches for a list of *points*,
not shapes of nonzero size.)
"""
mutable struct KDTree{K,S<:Shape{K}}
    s::Vector{S}
    ix::Int
    x::Float64
    left::KDTree{K,S}  # shapes ≤ x in coordinate ix
    right::KDTree{K,S} # shapes > x in coordinate ix
    KDTree{K,S}(s::AbstractVector{S}) where {K,S<:Shape{K}} = new(s, 0)
    function KDTree{K,S}(ix::Integer, x::Real, left::KDTree{K,S}, right::KDTree{K,S}) where {K,S<:Shape{K}}
        1 ≤ ix ≤ K || throw(BoundsError())
        new(S[], ix, x, left, right)
    end
end

Base.ndims(::KDTree{K}) where {K} = K

"""
    KDTree(s::AbstractVector{<:Shape{K}})

Construct a K-D tree (`KDTree`) representation of a list of
`shapes` in order to enable rapid searching of an shape list.

When searching the tree, shapes that appear earlier in `s`
take precedence over shapes that appear later.
"""
function KDTree(s::AbstractVector{S}) where {K,S<:Shape{K}}
    (length(s) <= 4 || K == 0) && return KDTree{K,S}(s)

    # figure out the best dimension ix to divide over,
    # the dividing plane x, and the number (nl,nr) of
    # shapes that fall into the left and right subtrees
    b = map(bounds, s) # cartesian bounding boxes of all shapes
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

    # don't bother subdividing if it doesn't reduce the # of shapes much
    4*min(nl,nr) > 3*length(s) && return KDTree{K,S}(s)

    # create the arrays of shapes in each subtree
    sl = Array{S}(nl)
    sr = Array{S}(nr)
    il = ir = 0
    for k in eachindex(s)
        if b[k][1][ix] ≤ x
            sl[il += 1] = s[k]
        end
        if b[k][2][ix] > x
            sr[ir += 1] = s[k]
        end
    end

    return KDTree{K,S}(ix, x, KDTree(sl), KDTree(sr))
end

depth(kd::KDTree) = kd.ix == 0 ? 0 : max(depth(kd.left), depth(kd.right)) + 1

Base.show(io::IO, kd::KDTree{K,S}) where {K,S} = print(io, "KDTree{$K,$S} of depth ", depth(kd))

function _show(io::IO, kd::KDTree, indent)
    indentstr = " "^indent
    if kd.ix == 0
        println(io, indentstr, length(kd.s), " shapes")
    else
        println(io, indentstr, "if x[", kd.ix, "] ≤ ", kd.x, ':')
        _show(io, kd.left, indent + 2)
        println(io, indentstr, "else:")
        _show(io, kd.right, indent + 2)
    end
end

function Base.show(io::IO, ::MIME"text/plain", kd::KDTree)
    println(io, kd, ':')
    _show(io, kd, 0)
end

function Base.findin(p::SVector{N}, s::Vector{S}) where {N,S<:Shape{N}}
    for i in eachindex(s)
        if p ∈ s[i]
            return Nullable{S}(s[i])
        end
    end
    return Nullable{S}()
end

function Base.findin(p::SVector{N}, kd::KDTree{N}) where {N}
    if isempty(kd.s)
        if p[kd.ix] ≤ kd.x
            return findin(p, kd.left)
        else
            return findin(p, kd.right)
        end
    else
        return findin(p, kd.s)
    end
end

"""
    findin(p::AbstractVector, kd::KDTree)

Return a `Nullable` container of the first shape in `kd` that contains
the point `p`, or an `isnull` container if no shape was found.
"""
Base.findin(p::AbstractVector, kd::KDTree{N}) where {N} = findin(SVector{N}(p), kd)
Base.findin(p::AbstractVector, s::Vector{<:Shape{N}}) where {N} = findin(SVector{N}(p), s)
