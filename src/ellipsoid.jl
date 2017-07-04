export Ellipsoid

type Ellipsoid{N,D} <: Shape{N}
    c::SVector{N,Float64} # Ellipsoid center
    p::SMatrix{N,N,Float64} # projection matrix to Ellipsoid coordinates
    ri2::SVector{N,Float64} # inverse square of "radius" (semi-axis) in each direction
    data::D             # auxiliary data
    (::Type{Ellipsoid{N,D}}){N,D}(c,ri2,p,data) = new{N,D}(c,ri2,p,data)  # inner constructor compatible with both v0.5 and v0.6
end

function Ellipsoid(c::AbstractVector, d::AbstractVector, axes=eye(length(c),length(c)), # columns are axes unit vectors
                   data=nothing)
    length(c) == length(d) == size(axes,1) == size(axes,2) || throw(DimensionMismatch())
    return Ellipsoid{length(c),typeof(data)}(c, inv(axes ./ sqrt.(sum(abs2,axes,1))), (d*0.5) .^ -2, data)
end

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
