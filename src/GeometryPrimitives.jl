module GeometryPrimitives

using Compat, StaticArrays

abstract type Object{N} end # a solid geometric object in N dimensions
Base.ndims(o::Object{N}) where N = N

export Object, normal, bounds

Base.in(x::AbstractVector, o::Object{N}) where N = SVector{N}(x) in o
normal(x::AbstractVector, o::Object{N}) where N = normal(SVector{N}(x), o)

include("sphere.jl")
include("box.jl")
include("ellipsoid.jl")
include("cylinder.jl")

include("kdtree.jl")

end # module
