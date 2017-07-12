export Sphere

mutable struct Sphere{N,D} <: Shape{N,D}
    c::SVector{N,Float64} # sphere center
    r::Float64          # radius
    data::D             # auxiliary data
    Sphere{N,D}(c,r,data) where {N,D} = new(c,r,data)  # suppress default outer constructor
end

Sphere(c::SVector{N}, r::Real, data::D=nothing) where {N,D} = Sphere{N,D}(c, r, data)
Sphere(c::AbstractVector, r::Real, data=nothing) = (N = length(c); Sphere(SVector{N}(c), r, data))

Base.:(==)(s1::Sphere, s2::Sphere) = s1.c==s2.c && s1.r==s2.r && s1.data==s2.data
Base.hash(s::Sphere, h::UInt) = hash(s.c, hash(s.r, hash(s.data, hash(:Sphere, h))))

Base.in(x::SVector{N}, s::Sphere{N}) where {N} = sum(abs2,x - s.c) ≤ s.r^2
normal(x::SVector{N}, s::Sphere{N}) where {N} = normalize(x - s.c)
bounds(s::Sphere) = (s.c-s.r, s.c+s.r)
