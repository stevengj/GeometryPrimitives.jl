module GeometryPrimitives

using Compat, FixedSizeArrays

abstract Object{N} # a solid geometric object in N dimensions

export Object, normal, bounds

Base.in(x::AbstractVector, o::Object) = Point(x) in o
normal(x::AbstractVector, o::Object) = normal(Point(x), o)

include("sphere.jl")
include("box.jl")
include("ellipsoid.jl")
include("cylinder.jl")

end # module
