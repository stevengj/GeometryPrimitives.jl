export Sphere

immutable Sphere{N,D} <: Shape{N}
    c::SVector{N,Float64} # sphere center
    r::Float64          # radius
    data::D             # auxiliary data
end
Sphere{D}(c::AbstractVector, r::Real, data::D=nothing) = Sphere{length(c),D}(c, r, data)
Base.in{N}(x::SVector{N}, s::Sphere{N}) = sum(abs2,x - s.c) ≤ s.r^2
normal{N}(x::SVector{N}, s::Sphere{N}) = normalize(x - s.c)
bounds(s::Sphere) = (s.c-s.r, s.c+s.r)
