export Cylinder

immutable Cylinder{N,D} <: Object{N}
    c::SVector{N,Float64} # Cylinder center
    a::SVector{N,Float64}   # axis unit vector
    r::Float64          # radius
    h2::Float64         # height * 0.5
    data::D             # auxiliary data
end
Cylinder{D}(c::AbstractVector, r::Real, a::AbstractVector, h::Real=Inf, data::D=nothing) =
    Cylinder{length(c),D}(c, normalize(a), r, h * 0.5, data)

function Base.in(x::SVector, s::Cylinder)
    d = x - s.c
    p = dot(d, s.a)
    abs(p) > s.h2 && return false
    return sumabs2(d - p*s.a) ≤ s.r^2
end

function normal(x::SVector, s::Cylinder)
    d = x - s.c
    p = dot(d, s.a)
    p > s.h2 && return s.a
    p < -s.h2 && return -s.a
    return normalize(d - p*s.a)
end

const rotate2 = @SMatrix [0.0 1.0; -1.0 0.0] # 2x2 90° rotation matrix
function endcircles(s::Cylinder{2})
    b = rotate2 * s.a
    axes = @SMatrix [s.a[1] b[1]; s.a[2] b[2]]
    d = 2*s.r
    return(Ellipsoid(s.c + s.a*s.h2, SVector(0.0, d), axes),
           Ellipsoid(s.c - s.a*s.h2, SVector(0.0, d), axes))
end

function endcircles(s::Cylinder{3})
    u = abs(s.a[3]) < abs(s.a[1]) ? SVector(0,0,1) : SVector(1,0,0)
    b1 = cross(s.a, u)
    b2 = cross(b1, s.a)
    axes = [s.a[1] b1[1] b2[1]; s.a[2] b1[2] b2[2]; s.a[3] b1[3] b2[3]]
    d = 2*s.r
    return(Ellipsoid(s.c + s.a*s.h2, SVector(0.0, d, d), axes),
           Ellipsoid(s.c - s.a*s.h2, SVector(0.0, d, d), axes))
end

function bounds(s::Cylinder)
    e1, e2 = endcircles(s)
    l1, u1 = bounds(e1)
    l2, u2 = bounds(e2)
    return min(l1,l2), max(u1,u2)
end
