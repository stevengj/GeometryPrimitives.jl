module GeometryPrimitives

using StaticArrays
using LinearAlgebra
using Statistics: mean
using Makie

export Shape, Shape1, Shape2, Shape3
export level, surfpt_nearby, normal, bounds, translate

abstract type Shape{N,N²} end # a solid geometric shape in N dimensions (N² = N*N is needed in some shapes, e.g., Cuboid)
const Shape1 = Shape{1,1}
const Shape2 = Shape{2,4}
const Shape3 = Shape{3,9}

Base.ndims(::Shape{N}) where {N} = N

# relative tolerance * x (assumed ≥ 0) for approximate comparisons,
# defined to square root of machine precision like Base.rtoldefault
rtol(x::AbstractFloat) = sqrt(eps(typeof(x))) * x

# The following functions return Any due to the limitations of Julia's dispatch system.
# Therefore, always call them with return type assertions.  See
# https://discourse.julialang.org/t/extending-base-in-type-stably/5341/12
# https://github.com/JuliaLang/julia/issues/23210
level(x::AbstractVector{<:Real}, s::Shape{N}) where {N} = level(SVector{N}(x), s)
Base.in(x::AbstractVector{<:Real}, s::Shape{N}) where {N} = level(x,s) ≥ 0
surfpt_nearby(x::AbstractVector{<:Real}, s::Shape{N}) where {N} = surfpt_nearby(SVector{N}(x), s)
normal(x::AbstractVector{<:Real}, s::Shape) = surfpt_nearby(x, s)[2]  # outward direction even for x inside s
translate(s::Shape{N}, ∆::AbstractVector{<:Real}) where {N} = translate(s, SVector{N}(∆))
translate(s::Shape{N}, ∆::SVector{N,<:Real}) where {N} = (s2 = deepcopy(s); s2.c += ∆; s2)  # default implementation

function orthoaxes(n::SVector{3,<:Real})
    u_temp = abs(n[3]) < abs(n[1]) ? SVector(0,0,1) : SVector(1,0,0)
    v = normalize(n × u_temp)
    u = v × n

    return u, v
end

function orthoaxes(n::SVector{N,<:Real}) where {N}
    u_temp = abs(n[3]) < abs(n[1]) ? SVector(0,0,1) : SVector(1,0,0)
    v = normalize(n × u_temp)
    u = v × n

    return u, v
end


include("hyper/hyper.jl")
include("planar/planar.jl")
include("prism/prism.jl")
include("util/util.jl")

end # module
