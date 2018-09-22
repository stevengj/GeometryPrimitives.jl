export Cylinder

mutable struct Cylinder{N,L,D} <: Shape{N,L,D}
    c::SVector{N,Float64}  # center
    r::Float64  # radius
    a::SVector{N,Float64}  # axis unit vector
    h2::Float64  # height * 0.5
    data::D  # auxiliary data
    Cylinder{N,L,D}(c,r,a,h2,data) where {N,L,D} = new(c,r,a,h2,data)  # suppress default outer constructor
end

Cylinder(c::SVector{N,<:Real}, r::Real, a::SVector{N,<:Real}, h::Real=Inf, data::D=nothing) where {N,D} =
    Cylinder{N,N*N,D}(c, r, normalize(a), 0.5h, data)

Cylinder(c::AbstractVector{<:Real}, r::Real, a::AbstractVector{<:Real}, h::Real=Inf, data=nothing) =
    (N = length(c); Cylinder(SVector{N}(c), r, SVector{N}(a), h, data))

Base.:(==)(s1::Cylinder, s2::Cylinder) = s1.c==s2.c && s1.r==s2.r && s1.a==s2.a && s1.h2==s2.h2 && s1.data==s2.data
Base.hash(s::Cylinder, h::UInt) = hash(s.c, hash(s.r, hash(s.a, hash(s.h2, hash(s.data, hash(:Cylinder, h))))))

function Base.in(x::SVector{N,<:Real}, s::Cylinder{N}) where {N}
    d = x - s.c
    p = d ⋅ s.a

    return abs(p) ≤ s.h2 && sum(abs2, d-p*s.a) ≤ s.r^2
end

function surfpt_nearby(x::SVector{N,<:Real}, s::Cylinder{N}) where {N}
    d = x - s.c
    p = d ⋅ s.a  # scalar
    q = d - p*s.a  # vector
    lp, lq = abs(p), norm(q)
    pout = copysign(1.0,p) * s.a
    qout = lq≠0 ? q/lq : @SVector(zeros(N))  # qout is used only when lq ≠ 0 below

    onbndp = abs(lp-s.h2) ≤ Base.rtoldefault(Float64) * s.h2
    onbndq = abs(lq-s.r) ≤ Base.rtoldefault(Float64) * s.r
    isoutp = (s.h2<lp) || onbndp
    isoutq = (s.r<lq) || onbndq
    ∆p, ∆q = s.h2-lp, s.r-lq
    if !isoutp && !isoutq  # x strictly inside cylinder
        (nout, ∆x) = (∆p≤∆q || lq==0) ? (pout, ∆p*pout) : (qout, ∆q*qout)  # qout is used when lq ≠ 0
    else  # x outside cylinder or on boundary
        ∆x = isoutp*∆p*pout + (isoutq && lq≠0)*∆q*qout  # qout is used when lq ≠ 0
        nout = (!isoutp || onbndp) && (!isoutq || onbndq) ? onbndp*pout + (onbndq && lq≠0)*qout : -∆x  # "if onbound in projected directions"
        nout = normalize(nout)
    end

    return x+∆x, nout
end

translate(s::Cylinder{N,L,D}, ∆::SVector{N,<:Real}) where {N,L,D} = Cylinder{N,L,D}(s.c+∆, s.r, s.a, s.h2, s.data)

const rotate2 = @SMatrix [0.0 1.0; -1.0 0.0] # 2x2 90° rotation matrix

function endcircles(s::Cylinder{2})
    b = rotate2 * s.a
    axes = [s.a b]  # SMatrix; no allocations
    return (Ellipsoid(s.c + s.a*s.h2, SVector(0.0, s.r), axes),
            Ellipsoid(s.c - s.a*s.h2, SVector(0.0, s.r), axes))
end

function endcircles(s::Cylinder{3})
    u, v = orthoaxes(s.a)
    axes = [u v s.a]  # SMatrix; no allocations
    return (Ellipsoid(s.c + s.a*s.h2, SVector(0.0, s.r, s.r), axes), # top disk
            Ellipsoid(s.c - s.a*s.h2, SVector(0.0, s.r, s.r), axes))  # bottom disk
end

function bounds(s::Cylinder)
    e1, e2 = endcircles(s)
    l1, u1 = bounds(e1)
    l2, u2 = bounds(e2)
    return min.(l1,l2), max.(u1,u2)
end
