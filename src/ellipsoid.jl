export Ellipsoid

immutable Ellipsoid{N,D} <: Object{N}
    c::Point{N,Float64} # Ellipsoid center
    ri2::Vec{N,Float64} # inverse square of "radius" (semi-axis) in each direction
    p::Mat{N,N,Float64} # projection matrix to Ellipsoid coordinates
    data::D             # auxiliary data
    Ellipsoid(c,ri2,p,data) = new(c,ri2,p,data)
end

function Ellipsoid(c, r, axes=eye(length(c),length(c)), # columns are axes unit vectors
                   data=nothing)
    length(c) == length(r) == size(axes,1) == size(axes,2) || throw(DimensionMismatch())
    return Ellipsoid{length(c),typeof(data)}(c, r .^ -2, inv(axes ./ sqrt(sumabs2(axes,2))),
                                             data)
end

Base.in{N}(x::Point, b::Ellipsoid{N}) = sum((b.p * (x - b.c)).^2 .* b.ri2) ≤ 1.0

normal(x::Point, b::Ellipsoid) = normalize(b.p' * (b.ri2 .* (b.p * (x - b.c))))

function bounds(b::Ellipsoid)
    # this is the bounding box for the axes-aligned box around the ellipsoid,
    # which may be bigger than necessary.  TODO: compute true bounding box
    a = abs(inv(b.p) * (b.ri2 .^ -0.5))
    return (b.c-a,b.c+a)
end
