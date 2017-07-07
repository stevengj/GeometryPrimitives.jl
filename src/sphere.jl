export Sphere

type Sphere{N,D} <: Shape{N}
    c::SVector{N,Float64} # sphere center
    r::Float64          # radius
    data::D             # auxiliary data
    (::Type{Sphere{N,D}}){N,D}(c,r,data) = new{N,D}(c,r,data)
end

Sphere{N,D}(c::SVector{N}, r::Real, data::D=nothing) = Sphere{N,D}(c, r, data)
Sphere(c::AbstractVector, r::Real, data=nothing) = (N = length(c); Sphere(SVector{N}(c), r, data))

Base.:(==)(s1::Sphere, s2::Sphere) = s1.c==s2.c && s1.r==s2.r && s1.data==s2.data
Base.hash(s::Sphere, h::UInt) = hash(s.c, hash(s.r, hash(s.data, hash(:Sphere, h))))

Base.in{N}(x::SVector{N}, s::Sphere{N}) = sum(abs2,x - s.c) ≤ s.r^2
normal{N}(x::SVector{N}, s::Sphere{N}) = normalize(x - s.c)
bounds(s::Sphere) = (s.c-s.r, s.c+s.r)
