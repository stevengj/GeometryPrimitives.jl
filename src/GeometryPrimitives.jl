module GeometryPrimitives

using Compat, StaticArrays

@compat abstract type Object{N} end # a solid geometric object in N dimensions
Base.ndims{N}(o::Object{N}) = N

export Object, normal, bounds

Base.in{N}(x::AbstractVector, o::Object{N}) = SVector{N}(x) in o
normal{N}(x::AbstractVector, o::Object{N}) = normal(SVector{N}(x), o)

include("sphere.jl")
include("box.jl")
include("ellipsoid.jl")
include("cylinder.jl")

include("kdtree.jl")

end # module
