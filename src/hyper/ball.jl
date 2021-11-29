export Ball

mutable struct Ball{N,N²} <: Shape{N,N²}
    c::SFloat{N}  # center of ball
    r::Float  # radius
    Ball{N,N²}(c,r) where {N,N²} = new(c,r)  # suppress default outer constructor
end

Ball(c::SReal{N}, r::Real) where {N} = Ball{N,N*N}(c, r)
Ball(c::AbsVecReal, r::Real) = (N = length(c); Ball(SVec{N}(c), r))

Base.:(==)(s1::Ball, s2::Ball) = s1.c==s2.c && s1.r==s2.r
Base.isapprox(s1::Ball, s2::Ball) = s1.c≈s2.c && s1.r≈s2.r
Base.hash(s::Ball, h::UInt) = hash(s.c, hash(s.r, hash(:Ball, h)))

level(x::SReal{N}, s::Ball{N}) where {N} = 1.0 - √(sum(abs2, x - s.c) / s.r^2)

function surfpt_nearby(x::SReal{N}, s::Ball{N}) where {N}
    nout = x==s.c ? SVec(ntuple(k -> k==1 ? 1.0 : 0.0, Val(N))) :  # nout = e₁ for x == s.c
                    normalize(x-s.c)
    return s.c+s.r*nout, nout
end

bounds(s::Ball) = (s.c.-s.r, s.c.+s.r)
