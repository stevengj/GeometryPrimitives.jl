export Ball

mutable struct Ball{N,N²,D} <: Shape{N,N²,D}
    c::SVector{N,Float64}  # center of ball
    r::Float64  # radius
    data::D  # auxiliary data
    Ball{N,N²,D}(c,r,data) where {N,N²,D} = new(c,r,data)  # suppress default outer constructor
end

Ball(c::SVector{N,<:Real}, r::Real, data::D=nothing) where {N,D} = Ball{N,N*N,D}(c, r, data)
Ball(c::AbstractVector{<:Real}, r::Real, data=nothing) = (N = length(c); Ball(SVector{N}(c), r, data))

Base.:(==)(s1::Ball, s2::Ball) = s1.c==s2.c && s1.r==s2.r && s1.data==s2.data
Base.isapprox(s1::Ball, s2::Ball) = s1.c≈s2.c && s1.r≈s2.r && s1.data==s2.data
Base.hash(s::Ball, h::UInt) = hash(s.c, hash(s.r, hash(s.data, hash(:Ball, h))))

Base.in(x::SVector{N,<:Real}, s::Ball{N}) where {N} = sum(abs2, x - s.c) ≤ s.r^2

function surfpt_nearby(x::SVector{N,<:Real}, s::Ball{N}) where {N}
    nout = x==s.c ? SVector(ntuple(k -> k==1 ? 1.0 : 0.0, Val(N))) :  # nout = e₁ for x == s.c
                    normalize(x-s.c)
    return s.c+s.r*nout, nout
end

translate(s::Ball{N,N²,D}, ∆::SVector{N,<:Real}) where {N,N²,D} = Ball{N,N²,D}(s.c+∆, s.r, s.data)

bounds(s::Ball) = (s.c.-s.r, s.c.+s.r)
