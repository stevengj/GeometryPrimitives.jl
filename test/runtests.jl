using GeometryPrimitives, StaticArrays, Base.Test
Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))

s = Sphere([3,4], 5)
@test [3,9] ∈ s
@test [3,9.1] ∉ s
@test normal([-1,2],s) == normalize([-1,2] - [3,4])
@test bounds(s) == ([-2,-1],[8,9])

b = Box([0,0], [2,4])
@test [0.3,-1.5] ∈ b
@test [0.3,-2.5] ∉ b
@test normal([1.1,0],b) == [1,0]
@test normal([-1.1,0],b) == [-1,0]
@test normal([1.1,2.01],b) == [0,1]
@test bounds(b) == ([-1,-2],[1,2])
@test bounds(Box([0,0], [2,4], [1 1; 1 -1])) ≈ ([-3*√0.5,-3*√0.5], [3*√0.5,3*√0.5])

c = Cylinder([0,0,0], 0.3, [0,0,1], 2.2)
@test [0.2,0.2,1] ∈ c
@test SVector(0.2,0.2,1.2) ∉ c
@test [0.2,0.25,1] ∉ c
@test normal([0.1,0.2,-1.3], c) == [0,0,-1]
@test normal([0.31, 0, 0.3], c) == [1,0,0]
@test bounds(c) == ([-0.3,-0.3,-1.1],[0.3,0.3,1.1])

o = Object{2}[Sphere([i,0], 1, i) for i in 0:20]
kd = KDTree(o)
@test GeometryPrimitives.depth(kd) == 4
@test get(findin([10.1,0], kd)).data == 10
@test isnull(findin([10.1,1], kd))
