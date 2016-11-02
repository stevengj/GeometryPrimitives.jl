module GeometryPrimitives

using Compat, StaticArrays

abstract Object{N} # a solid geometric object in N dimensions

export Object, normal, bounds

Base.in{N}(x::AbstractVector, o::Object{N}) = SVector{N}(x) in o
normal{N}(x::AbstractVector, o::Object{N}) = normal(SVector{N}(x), o)

include("sphere.jl")
include("box.jl")
include("ellipsoid.jl")
include("cylinder.jl")

end # module
