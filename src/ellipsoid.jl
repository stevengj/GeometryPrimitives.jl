export Ellipsoid

type Ellipsoid{N,D,L} <: Shape{N}
    c::SVector{N,Float64} # Ellipsoid center
    ri2::SVector{N,Float64} # inverse square of "radius" (semi-axis) in each direction
    p::SMatrix{N,N,Float64,L} # projection matrix to Ellipsoid coordinates
    data::D             # auxiliary data
    (::Type{Ellipsoid{N,D,L}}){N,D,L}(c,ri2,p,data) = new{N,D,L}(c,ri2,p,data)
end

Ellipsoid{N,D,L,R<:Real}(c::SVector{N}, d::SVector{N},
                         axes::SMatrix{N,N,R,L}=@SMatrix(eye(N)),  # columns are axes unit vectors
                         data::D=nothing) =
    Ellipsoid{N,D,L}(c, (d*0.5) .^ -2, inv(axes ./ sqrt.(sum(abs2,axes,1))), data)

Ellipsoid(c::AbstractVector, d::AbstractVector, axes::AbstractMatrix=eye(length(c)), data=nothing) =
    (N = length(c); Ellipsoid(SVector{N}(c), SVector{N}(d), SMatrix{N,N}(axes), data))

Base.:(==)(b1::Ellipsoid, b2::Ellipsoid) = b1.c==b2.c && b1.ri2==b2.ri2 && b1.p==b2.p && b1.data==b2.data
Base.hash(b::Ellipsoid, h::UInt) = hash(b.c, hash(b.ri2, hash(b.p, hash(b.data, hash(:Ellipsoid, h)))))

Base.in{N}(x::SVector{N}, b::Ellipsoid{N}) = sum((b.p * (x - b.c)).^2 .* b.ri2) ≤ 1.0

normal{N}(x::SVector{N}, b::Ellipsoid{N}) = normalize(Ac_mul_B(b.p, b.ri2 .* (b.p * (x - b.c))))

function boundpts{N}(b::Ellipsoid{N})
    # Return the points tangential to the bounding box.
    # For N = 3, it returns three points at which the direction normals are +x, +y, +z
    # directions, respectively.

    r2 = 1./b.ri2
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

# NaN-ignoring max.
nanmax(x,y) = isnan(x) ? y : (isnan(y) ? x : max(x,y))

function bounds{N}(b::Ellipsoid{N})
    M = boundpts(b)

    # Note that when M does not include NaN, we can simply set m = diag(M), because the
    # first (second, third) column of M is the point on the ellipsoid at which the normal
    # vector is +x (+y, +z).  Therefore, the x (y, z) point of the first (second, third)
    # column has the largest x (y, z) coordinate.
    #
    # However, if one of a, b, c is zero, the shape is a disk.  Then one column of M can be
    # completely filled with NaN.  This column must not be counted as a bounding point, so
    # we apply nanmax by StaticArrays.reducedim along the row direction.
    m = reducedim(nanmax, M, Val{2})[:,1]

    return (b.c-m,b.c+m)
end
