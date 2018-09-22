module GeometryPrimitives
using Compat, StaticArrays, LinearAlgebra
using Statistics: mean
export Shape, surfpt_nearby, normal, bounds, translate

abstract type Shape{N,L,D} end # a solid geometric shape in N dimensions (L = N*N is needed in some shapes, e.g., Box)

Base.ndims(o::Shape{N}) where {N} = N

# The following functions return Any due to the limitations of Julia's dispatch system.
# Therefore, always call them with type assertions.  See
# https://discourse.julialang.org/t/extending-base-in-type-stably/5341/12
# https://github.com/JuliaLang/julia/issues/23210
Base.in(x::AbstractVector{<:Real}, o::Shape{N}) where {N} = SVector{N}(x) in o
surfpt_nearby(x::AbstractVector{<:Real}, o::Shape{N}) where {N} = surfpt_nearby(SVector{N}(x), o)
normal(x::AbstractVector{<:Real}, o::Shape) = surfpt_nearby(x, o)[2]  # outward direction even for x inside o
translate(o::Shape{N}, ∆::AbstractVector{<:Real}) where {N} = translate(o, SVector{N}(∆))

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


include("box.jl")
include("ellipsoid.jl")
include("sphere.jl")
include("prism/prism.jl")
include("periodize.jl")
include("kdtree.jl")
include("vxlcut.jl")

end # module
