export Ellipsoid

struct Ellipsoid{N,D} <: Object{N}
    c::SVector{N,Float64} # Ellipsoid center
    p::SMatrix{N,N,Float64} # projection matrix to Ellipsoid coordinates
    ri2::SVector{N,Float64} # inverse square of "radius" (semi-axis) in each direction
    data::D             # auxiliary data
    Ellipsoid{N,D}(c,ri2,p,data) where {N,D} = new(c,ri2,p,data)
end

function Ellipsoid(c::AbstractVector, d::AbstractVector, axes=eye(length(c),length(c)), # columns are axes unit vectors
                   data=nothing)
    length(c) == length(d) == size(axes,1) == size(axes,2) || throw(DimensionMismatch())
    return Ellipsoid{length(c),typeof(data)}(c, inv(axes ./ sqrt.(sum(abs2,axes,2))), (d*0.5) .^ -2, data)
end

Base.in(x::SVector{N}, b::Ellipsoid{N}) where N = sum((b.p * (x - b.c)).^2 .* b.ri2) ≤ 1.0

normal(x::SVector{N}, b::Ellipsoid{N}) where N = normalize(Ac_mul_B(b.p, b.ri2 .* (b.p * (x - b.c))))

function bounds(b::Ellipsoid)
    # this is the bounding box for the axes-aligned box around the ellipsoid,
    # which may be bigger than necessary.  TODO: compute true bounding box
    A = inv(b.p) .* (b.ri2 .^ -0.5)' # array of scaled axes vectors.
    m = maximum(A * signmatrix(b), 2)[:,1] # extrema of all 2^N corners of the box
    return (b.c-m,b.c+m)
end
