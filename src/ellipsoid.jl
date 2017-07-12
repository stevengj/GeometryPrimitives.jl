export Ellipsoid

mutable struct Ellipsoid{N,D,L} <: Shape{N,D}
    c::SVector{N,Float64} # Ellipsoid center
    ri2::SVector{N,Float64} # inverse square of "radius" (semi-axis) in each direction
    p::SMatrix{N,N,Float64,L} # projection matrix to Ellipsoid coordinates
    data::D             # auxiliary data
    Ellipsoid{N,D,L}(c,ri2,p,data) where {N,D,L} = new(c,ri2,p,data)
end

Ellipsoid(c::SVector{N}, d::SVector{N},
          axes::SMatrix{N,N,<:Real,L}=@SMatrix(eye(N)),  # columns are axes unit vectors
          data::D=nothing) where {N,D,L} =
    Ellipsoid{N,D,L}(c, (0.5d) .^ -2, inv((axes' ./ sqrt.(sum(abs2,axes,Val{1}))[1,:])'), data)
# Use this after StaticArrays issue 242 is fixed:
#    Ellipsoid{N,D,L}(c, (0.5d) .^ -2, inv(axes ./ sqrt.(sum(abs2,axes,Val{1}))), data)

Ellipsoid(c::AbstractVector, d::AbstractVector, axes::AbstractMatrix=eye(length(c)), data=nothing) =
    (N = length(c); Ellipsoid(SVector{N}(c), SVector{N}(d), SMatrix{N,N}(axes), data))

Base.:(==)(b1::Ellipsoid, b2::Ellipsoid) = b1.c==b2.c && b1.ri2==b2.ri2 && b1.p==b2.p && b1.data==b2.data
Base.hash(b::Ellipsoid, h::UInt) = hash(b.c, hash(b.ri2, hash(b.p, hash(b.data, hash(:Ellipsoid, h)))))

Base.in(x::SVector{N}, b::Ellipsoid{N}) where {N} = sum((b.p * (x - b.c)).^2 .* b.ri2) ≤ 1.0

normal(x::SVector{N}, b::Ellipsoid{N}) where {N} = normalize(Ac_mul_B(b.p, b.ri2 .* (b.p * (x - b.c))))

function boundpts(b::Ellipsoid{N}) where {N}
    # Return the points tangential to the bounding box.
    # For N = 3, it returns three points at which the direction normals are +x, +y, +z
    # directions, respectively.

    r2 = 1 ./ b.ri2
    ndir = b.p  # Cartesian directions in ellipsoid coordinates: b.p * eye(3)

    # In the ellipsoid coordinates, the point on the ellipsoid where the direction normal is
    # n is (n .* r2) / sqrt(r2' * n.^2).  Below is the broadcasted version of this over a
    # matrix n, whose each column is a direction normal.  Once calculated, we need to
    # change the coordinates back to the original coordinates.
    #
    # This operation can be written M = b.p' * ((ndir .* r2) ./ sqrt.(r2' * ndir.^2)), but
    # the resulting M is not an SMatrix.  (The calculation involves broadcasted division by
    # a row vector, which leads to a non-SMatrix.)  Therefore, we first calculate M'
    # (which remains SMatrix) and recover M.
    M = (((ndir .* r2)' ./ sqrt.(ndir'.^2 * r2)) * b.p)'

    return M
end

function bounds(b::Ellipsoid{N}) where {N}
    M = boundpts(b)

    # Note that when M does not include NaN, we can simply set m = diag(M), because the
    # first (second, third) column of M is the point on the ellipsoid at which the normal
    # vector is +x (+y, +z).  Therefore, the x (y, z) point of the first (second, third)
    # column has the largest x (y, z) coordinate.
    #
    # However, if one of a, b, c is zero, the shape is a disk.  Then one column of M can be
    # completely filled with NaN.  This column must not be counted as a bounding point, so
    # we apply NaN-ignoring maximum by StaticArrays.reducedim along the row direction.
    #
    # For the efficient implementation of NaN-ignoring maximum, see
    # https://discourse.julialang.org/t/inconsistency-between-findmax-and-maximum-with-respect-to-nan/4214/8
    m = reducedim((x,y) -> x<y ? y : x, M, Val{2}, -Inf)[:,1]

    return (b.c-m,b.c+m)
end
