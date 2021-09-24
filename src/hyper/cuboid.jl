export Cuboid

# Below, the cuboid axes describe the directions of edges.  For example, if a₁, a₂, a₃ are
# the axes of a 3D cuboid, then the volume inside the cuboid is the collection of points
# x₁a₁ + x₂a₂ + x₃a₃, where -rₙ ≤ xₙ ≤ rₙ.  This means aₙ and aₘ with n ≠ m span a face of
# the cuboid.
#
# The cuboid projection matrix (Cuboid.p) is the matrix that produces the coordinates
# (x₁, x₂, x₃) of a given point.  In other words, for a point x inside the cuboid, the nth
# entry of Cuboid.p * x must have magnitude ≤ rₙ.  This means that Cuboid.p * aₙ = eₙ,
# because aₙ has only the aₙ-component and the component must be 1.  This also means that
# each row of Cuboid.p is orthogononal to two cuboid axes, and therefore normal to the face
# spanned by the two cuboid axes.  (Note that the rows of Cuboid.p are not unit normals to
# the faces, becuase they are not unit vectors.)
mutable struct Cuboid{N,N²,D} <: Shape{N,N²,D}
    c::SVector{N,Float64}  # center of cuboid
    r::SVector{N,Float64}  # "radii" (semi-axes) in axis directions
    p::SMatrix{N,N,Float64,N²}  # projection matrix to cuboid coordinates
    data::D  # auxiliary data
    Cuboid{N,N²,D}(c,r,p,data) where {N,N²,D} = new(c,r,p,data)  # suppress default outer constructor
end

Cuboid(c::SVector{N,<:Real},
    d::SVector{N,<:Real},
    axes::SMatrix{N,N,<:Real}=SMatrix{N,N,Float64}(I),
    data::D=nothing) where {N,D} =
    Cuboid{N,N*N,D}(c, 0.5d, inv(axes ./ sqrt.(sum(abs2,axes,dims=Val(1)))), data)

Cuboid(c::AbstractVector{<:Real},  # center of cuboid
    d::AbstractVector{<:Real},  # size of cuboid in axis directions
    axes::AbstractMatrix{<:Real}=Matrix{Float64}(I,length(c),length(c)),  # columns are axes vectors (each being parallel to two sets of faces in 3D)
    data=nothing) =
    (N = length(c); Cuboid(SVector{N}(c), SVector{N}(d), SMatrix{N,N}(axes), data))

Base.:(==)(b1::Cuboid, b2::Cuboid) = b1.c==b2.c && b1.r==b2.r && b1.p==b2.p && b1.data==b2.data
Base.isapprox(b1::Cuboid, b2::Cuboid) = b1.c≈b2.c && b1.r≈b2.r && b1.p≈b2.p && b1.data==b2.data
Base.hash(b::Cuboid, h::UInt) = hash(b.c, hash(b.r, hash(b.p, hash(b.data, hash(:Cuboid, h)))))

function level(x::SVector{N,<:Real}, b::Cuboid{N}) where {N}
    d = b.p * (x - b.c)

    return maximum(abs.(d) ./ b.r) - 1.0
end

function surfpt_nearby(x::SVector{N,<:Real}, b::Cuboid{N}) where {N}
    ax = inv(b.p)  # axes: columns are unit vectors

    # Below, the rows of n are the unit normals to the faces of the cuboid.  Appropriate
    # signs will be multiplied later.  (Note that the signs will be such that the rows of n
    # are always outward directions, even if x is inside the cuboid.)
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
    onbnd = abs.(b.r.-absd) .≤ Base.rtoldefault(Float64) .* b.r  # basically b.r .≈ absd but faster
    isout = (b.r.<absd) .| onbnd
    ∆ = (b.r .- absd) .* cosθ  # entries can be negative
    if count(isout) == 0  # x strictly inside cuboid; ∆ all positive
        l∆x, i = findmin(∆)  # find closest face
        nout = n[i,:]
        ∆x = l∆x * nout
    else  # x outside cuboid or on boundary in one or multiple directions
        ∆x = n' * (∆ .* isout)  # project out .!isout directions

        # Below, (.!isout .| onbnd) tests if each entry of ∆ (not ∆x) is either projected
        # out or on the boundary.  If the ith entry of ∆ is projected out, ∆x (not ∆) does
        # not have any component in the n[i,:] direction.  Even if the ith entry of ∆ is not
        # projected out, if it is too small then ∆x barely has any component in the n[i,:]
        # direction.  If all the entries of ∆ satisfy one of the two conditions, ∆x ≈ 0 and
        # we cannot use -∆x for nout.  In that case, take n'*onbnd as nout.  When x is
        # outside the cuboid only in one dimension, n'*onbnd is basically the row of n along
        # that dimension, which reduces to normalize(-∆x) even if ∆x is very small.
        nout = all(.!isout .| onbnd) ? n'*onbnd : -∆x
        nout = normalize(nout)
    end

    return x+∆x, nout
end

translate(b::Cuboid{N,N²,D}, ∆::SVector{N,<:Real}) where {N,N²,D} = Cuboid{N,N²,D}(b.c+∆, b.r, b.p, b.data)

signmatrix(b::Cuboid{1}) = SMatrix{1,1}(1)
signmatrix(b::Cuboid{2}) = SMatrix{2,2}(1,1, -1,1)
signmatrix(b::Cuboid{3}) = SMatrix{3,4}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1)

function bounds(b::Cuboid)
    A = inv(b.p) .* b.r'
    m = maximum(abs.(A * signmatrix(b)), dims=Val(2))[:,1] # extrema of all 2^N corners of the cuboid
    return (b.c-m,b.c+m)
end
