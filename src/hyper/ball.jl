export Ball

mutable struct Ball{N,N²} <: Shape{N,N²}
    c::SVector{N,Float64}  # center of ball
    r::Float64  # radius
    Ball{N,N²}(c,r) where {N,N²} = new(c,r)  # suppress default outer constructor
end

Ball(c::SVector{N,<:Real}, r::Real) where {N} = Ball{N,N*N}(c, r)
Ball(c::AbstractVector{<:Real}, r::Real) = (N = length(c); Ball(SVector{N}(c), r))

Base.:(==)(s1::Ball, s2::Ball) = s1.c==s2.c && s1.r==s2.r
Base.isapprox(s1::Ball, s2::Ball) = s1.c≈s2.c && s1.r≈s2.r
Base.hash(s::Ball, h::UInt) = hash(s.c, hash(s.r, hash(:Ball, h)))

level(x::SVector{N,<:Real}, s::Ball{N}) where {N} = √(sum(abs2, x - s.c) / s.r^2) - 1.0

function surfpt_nearby(x::SVector{N,<:Real}, s::Ball{N}) where {N}
    nout = x==s.c ? SVector(ntuple(k -> k==1 ? 1.0 : 0.0, Val(N))) :  # nout = e₁ for x == s.c
                    normalize(x-s.c)
    return s.c+s.r*nout, nout
end

bounds(s::Ball) = (s.c.-s.r, s.c.+s.r)
