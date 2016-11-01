using GeometryPrimitives, FixedSizeArrays, Base.Test

s = Sphere([3,4], 5)
@test [3,9] ∈ s
@test [3,9.1] ∉ s

b = Box([0,0], [1,2])
@test [0.3,-1.5] ∈ b
@test [0.3,-2.5] ∉ b

c = Cylinder([0,0,0], 0.3, [0,0,1], 2.2)
@test Point(0.2,0.2,1) ∈ c
@test Point(0.2,0.2,1.2) ∉ c
@test Point(0.2,0.25,1) ∉ c
