export Sphere

struct Sphere{N,D} <: Object{N}
    c::SVector{N,Float64} # sphere center
    r::Float64          # radius
    data::D             # auxiliary data
end
Sphere(c::AbstractVector, r::Real, data::D=nothing) where D = Sphere{length(c),D}(c, r, data)
Base.in(x::SVector{N}, s::Sphere{N}) where N = sum(abs2, x - s.c) ≤ s.r^2
normal(x::SVector{N}, s::Sphere{N}) where N = normalize(x - s.c)
bounds(s::Sphere) = (s.c-s.r, s.c+s.r)
