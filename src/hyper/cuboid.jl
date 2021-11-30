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
mutable struct Cuboid{N,N²} <: Shape{N,N²}
    c::SFloat{N}  # center of cuboid
    r::SFloat{N}  # "radii" (semi-axes) in axis directions
    p::S²Float{N,N²}  # projection matrix to cuboid coordinates
    Cuboid{N,N²}(c,r,p) where {N,N²} = new(c,r,p)  # suppress default outer constructor
end

Cuboid(c::SReal{N},
       s::SReal{N},
       axes::S²Real{N}=S²Float{N}(I)
       ) where {N} =
    Cuboid{N,N*N}(c, 0.5s, inv(axes ./ sqrt.(sum(abs2,axes,dims=Val(1)))))

Cuboid(c::AbsVecReal,  # center of cuboid
       s::AbsVecReal,  # size of cuboid in axis directions
       axes::AbsMatReal=MatFloat(I,length(c),length(c))) =  # columns are axes vectors (each being parallel to two sets of faces in 3D)
    (N = length(c); Cuboid(SVec{N}(c), SVec{N}(s), S²Mat{N}(axes)))

Cuboid(d::Tuple2{AbsVecReal}) =  # end points of diagonal
    Cuboid((d[1]+d[2])/2, abs.(d[2]-d[1]))

Base.:(==)(s1::Cuboid, s2::Cuboid) = s1.c==s2.c && s1.r==s2.r && s1.p==s2.p
Base.isapprox(s1::Cuboid, s2::Cuboid) = s1.c≈s2.c && s1.r≈s2.r && s1.p≈s2.p
Base.hash(s::Cuboid, h::UInt) = hash(s.c, hash(s.r, hash(s.p, hash(:Cuboid, h))))

function level(x::SReal{N}, s::Cuboid{N}) where {N}
    d = s.p * (x - s.c)

    return 1.0 - maximum(abs.(d) ./ s.r)
end

function surfpt_nearby(x::SReal{N}, s::Cuboid{N}) where {N}
    ax = inv(s.p)  # axes: columns are unit vectors

    # Below, the rows of n are the unit normals to the faces of the cuboid.  Appropriate
    # signs will be multiplied later.  (Note that the signs will be such that the rows of n
    # are always outward directions, even if x is inside the cuboid.)
    n = (s.p ./ sqrt.(sum(abs2,s.p,dims=Val(2))[:,1]))  # s.p normalized in row direction

    # Below, θ[i], the angle between ax[:,i] and n[i,:], is always acute (i.e, cosθ .≥ 0),
    # because the diagonal entries of s.p * ax are positive (= 1) and the diagonal entries
    # of n * ax are the scaled version of the diagonal entries of s.p * ax with positive
    # scale factors.
    cosθ = sum(ax.*n', dims=Val(1))[1,:]  # equivalent to diag(n*ax)
    # cosθ = diag(n*ax)  # faster than SVec(ntuple(i -> ax[:,i]⋅n[i,:], Val(N)))
    # @assert all(cosθ .≥ 0)

    d = s.p * (x - s.c)
    n = n .* copysign.(1.0,d)  # operation returns SMatrix (reason for leaving n untransposed)
    absd = abs.(d)
    onbnd = abs.(s.r.-absd) .≤ Base.rtoldefault(Float) .* s.r  # basically s.r .≈ absd but faster
    isout = (s.r.<absd) .| onbnd
    ∆ = (s.r .- absd) .* cosθ  # entries can be negative
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

signmatrix(::Cuboid{1}) = S²Mat{1}(1)
signmatrix(::Cuboid{2}) = S²Mat{2}(1,1, -1,1)
signmatrix(::Cuboid{3}) = SMat{3,4}(1,1,1, -1,1,1, 1,-1,1, 1,1,-1)

function bounds(s::Cuboid)
    A = inv(s.p) .* s.r'
    m = maximum(abs.(A * signmatrix(s)), dims=Val(2))[:,1] # extrema of all 2^N corners of the cuboid
    return (s.c-m,s.c+m)
end
