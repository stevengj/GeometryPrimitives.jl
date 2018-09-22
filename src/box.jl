export Box

# Below, the box axes describe the directions of edges.  For example, if a₁, a₂, a₃ are the
# axes of a 3D box, then the volume inside the box is the collection of points x₁a₁ + x₂a₂
# + x₃a₃, where -rₙ ≤ xₙ ≤ rₙ.  This means aₙ and aₘ with n ≠ m span a face of the box.
#
# The box projection matrix (Box.p) is the matrix that produces the coordinates (x₁, x₂, x₃)
# of a given point.  In other words, for a point x inside the box, the nth entry of Box.p * x
# must have magnitude ≤ rₙ.  This means that Box.p * aₙ = eₙ, because aₙ has only the aₙ-component
# and the component must be 1.  This also means that each row of Box.p is orthogononal to
# two box axes, and therefore normal to the face spanned by the two box axes.  (Note that
# the rows of Box.p are not unit normals to the faces, becuase they are not unit vectors.)
mutable struct Box{N,L,D} <: Shape{N,L,D}
    c::SVector{N,Float64}  # center of box
    r::SVector{N,Float64}  # "radii" (semi-axes) in axis directions
    p::SMatrix{N,N,Float64,L}  # projection matrix to box coordinates
    data::D  # auxiliary data
    Box{N,L,D}(c,r,p,data) where {N,L,D} = new(c,r,p,data)  # suppress default outer constructor
end

Box(c::SVector{N,<:Real},
    d::SVector{N,<:Real},
    axes::SMatrix{N,N,<:Real}=SMatrix{N,N,Float64}(I),
    data::D=nothing) where {N,D} =
    Box{N,N*N,D}(c, 0.5d, inv(axes ./ sqrt.(sum(abs2,axes,dims=Val(1)))), data)

Box(c::AbstractVector{<:Real},  # center of box
    d::AbstractVector{<:Real},  # size of box in axis directions
    axes::AbstractMatrix{<:Real}=Matrix{Float64}(I,length(c),length(c)),  # columns are axes vectors (each being parallel to two sets of faces in 3D)
    data=nothing) =
    (N = length(c); Box(SVector{N}(c), SVector{N}(d), SMatrix{N,N}(axes), data))

Base.:(==)(b1::Box, b2::Box) = b1.c==b2.c && b1.r==b2.r && b1.p==b2.p && b1.data==b2.data
Base.isapprox(b1::Box, b2::Box) = b1.c≈b2.c && b1.r≈b2.r && b1.p≈b2.p && b1.data==b2.data
Base.hash(b::Box, h::UInt) = hash(b.c, hash(b.r, hash(b.p, hash(b.data, hash(:Box, h)))))

function Base.in(x::SVector{N,<:Real}, b::Box{N}) where {N}
    d = b.p * (x - b.c)
    for i = 1:N
        abs(d[i]) > b.r[i] && return false  # boundary is considered inside
    end
    return true
end

function surfpt_nearby(x::SVector{N,<:Real}, b::Box{N}) where {N}
    ax = inv(b.p)  # axes: columns are unit vectors

    # Below, the rows of n are the unit normals to the faces of the box.  Appropriate signs
    # will be multiplied later.  (Note that the signs will be such that the rows of n are
    # always outward directions, even if x is inside the box.)
    n = (b.p ./ sqrt.(sum(abs2,b.p,dims=Val(2))[:,1]))  # b.p normalized in row direction

    # Below, θ[i], the angle between ax[:,i] and n[i,:], is always acute (i.e, cosθ .≥ 0),
    # because the diagonal entries of b.p * ax are positive (= 1) and the diagonal entries
    # of n * ax are the scaled version of the diagonal entries of b.p * ax with positive
    # scale factors.
    cosθ = sum(ax.*n', dims=Val(1))[1,:]  # equivalent to diag(n*ax)
    # cosθ = diag(n*ax)  # faster than SVector(ntuple(i -> ax[:,i]⋅n[i,:], Val(N)))
    # @assert all(cosθ .≥ 0)

    d = b.p * (x - b.c)
    n = n .* copysign.(1.0,d)  # operation returns SMatrix (reason for leaving n untransposed)
    absd = abs.(d)
    onbound = abs.(b.r.-absd) .≤ Base.rtoldefault(Float64) .* b.r  # basically b.r .≈ absd but faster
    isout = (b.r.<absd) .| onbound
    ∆ = (b.r .- absd) .* cosθ  # entries can be negative
    if count(isout) == 0  # x strictly inside box; ∆ all positive
        l∆x, i = findmin(∆)  # find closest face
        nout = n[i,:]
        ∆x = l∆x * nout
    else  # x outside box or on boundary in one or multiple directions
        ∆x = n' * (∆ .* isout)  # project out .!isout directions

        # Below, (.!isout .| onbound) tests if each entry of ∆ (not ∆x) is either projected
        # out or on the boundary.  If the ith entry of ∆ is projected out, ∆x (not ∆) does
        # not have any component in the n[i,:] direction.  Even if the ith entry of ∆ is not
        # projected out, if it is too small then ∆x barely has any component in the n[i,:]
        # direction.  If all the entries of ∆ satisfy one of the two conditions, ∆x ≈ 0 and
        # we cannot use -∆x for nout.  In that case, take n'*onbound as nout.  When x is
        # outside the box only in one dimension, n'*onbound is basically the row of n along
        # that dimension, which reduces to normalize(-∆x) even if ∆x is very small.
        nout = all(.!isout .| onbound) ? n'*onbound : -∆x
        nout = normalize(nout)
    end

    return x+∆x, nout
end

translate(b::Box{N,L,D}, ∆::SVector{N,<:Real}) where {N,L,D} = Box{N,L,D}(b.c+∆, b.r, b.p, b.data)

signmatrix(b::Box{1}) = SMatrix{1,1}(1)
signmatrix(b::Box{2}) = SMatrix{2,2}(1,1, -1,1)
signmatrix(b::Box{3}) = SMatrix{3,4}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1)

function bounds(b::Box)
    A = inv(b.p) .* b.r'
    m = maximum(abs.(A * signmatrix(b)), dims=Val(2))[:,1] # extrema of all 2^N corners of the box
    return (b.c-m,b.c+m)
end
