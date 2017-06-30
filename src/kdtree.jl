# kd tree for fast searching of a list of shapes
export KDTree

"""
A K-D tree of shapes in `K` dimensions.

This is a binary tree.  Each node in the tree divides space along
coordinate `ix` into shapes that overlap the region `≤ x` and `≥ x`.
The leaves of the tree are just a list of shapes `o`.

This implementation is a little different from a standard K-D tree,
because a given shape may be in both branches of the tree.  (Normally,
K-D trees are used for nearest-neighbor searches for a list of *points*,
not shapes of nonzero size.)
"""
type KDTree{K}
    o::Vector{Shape{K}}
    ix::Int
    x::Float64
    left::KDTree  # shapes ≤ x in coordinate ix
    right::KDTree # shapes > x in coordinate ix
    (::Type{KDTree{K}}){K}(o::AbstractVector{Shape{K}}) = new{K}(o, 0)  # inner constructor compatible with both v0.5 and v0.6
    function (::Type{KDTree{K}}){K}(x::Real, ix::Integer, left::KDTree{K}, right::KDTree{K})  # inner constructor compatible with both v0.5 and v0.6
        1 ≤ ix ≤ K || throw(BoundsError())
        new{K}(Shape{K}[], ix, x, left, right)
    end
end

Base.ndims{K}(::KDTree{K}) = K

"""
    KDTree(shapes::AbstractVector{Shape{K}})

Construct a K-D tree (`KDTree`) representation of a list of
`shapes` in order to enable rapid searching of an shape list.

When searching the tree, shapes that appear earlier in `shapes`
take precedence over shapes that appear later.
"""
function KDTree{K}(o::AbstractVector{Shape{K}})
    (length(o) <= 4 || K == 0) && return KDTree{K}(o)

    # figure out the best dimension ix to divide over,
    # the dividing plane x, and the number (nl,nr) of
    # shapes that fall into the left and right subtrees
    b = map(bounds, o) # cartesian bounding boxes of all shapes
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
    4*min(nl,nr) > 3*length(o) && return KDTree{K}(o)

    # create the arrays of shapes in each subtree
    ol = Array{Shape{K}}(nl)
    or = Array{Shape{K}}(nr)
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
        println(io, indentstr, length(kd.o), " shapes")
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

function Base.findin{N}(p::SVector{N}, o::Vector{Shape{N}})
    for i in eachindex(o)
        if p ∈ o[i]
            return Nullable{Shape{N}}(o[i])
        end
    end
    return Nullable{Shape{N}}()
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

Return a `Nullable` container of the first shape in `kd` that contains
the point `p`, or an `isnull` container if no shape was found.
"""
Base.findin{N}(p::AbstractVector, kd::KDTree{N}) = findin(SVector{N}(p), kd)
Base.findin{N}(p::AbstractVector, o::Vector{Shape{N}}) = findin(SVector{N}(p), o)
