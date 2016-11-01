export Sphere

immutable Sphere{N,D} <: Object{N}
    c::Point{N,Float64} # sphere center
    r::Float64          # radius
    data::D             # auxiliary data
end
Sphere{D}(c, r::Real, data::D=nothing) = Sphere{length(c),D}(c, r, data)
Base.in(x::Point, s::Sphere) = sumabs2(x - s.c) ≤ s.r^2
normal(x::Point, s::Sphere) = normalize(x - s.c)
bounds(s::Sphere) = (s.c-s.r, s.c+s.r)
