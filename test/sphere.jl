@testset "Sphere" begin
    s = Sphere([3,4], 5)
    @test s == deepcopy(s)
    @test hash(s) == hash(deepcopy(s))
    @test ndims(s) == 2
    @test [3,9] ∈ s
    @test [3,9.1] ∉ s

    @test ((x,nout) = surfpt_nearby([3,4],s); norm(x-[3,4])≈5 && norm(nout)≈1)  # handle point at center properly
    @test all([surfpt_nearby([3+ρ*5sx,4],s)[1] ≈ [3+5sx,4] for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,1)])
    @test all([surfpt_nearby([3,4+ρ*5sy],s)[1] ≈ [3,4+5sy] for ρ = (one⁻⁻,1,one⁺⁺), sy = (-1,1)])
    @test all([surfpt_nearby([3+ρ*5sx/√2,4+ρ*5sy/√2],s)[1] ≈ [3+5sx/√2,4+5sy/√2]
               for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1), sy = (-1,1)])
    @test all([(p = [3+5sx,4]; surfpt_nearby(p,s) ≈ (p,[sx,0])) for sx = (-1,1)])  # handle point on boundary properly
    @test all([(p = [3,4+5sy]; surfpt_nearby(p,s) ≈ (p,[0,sy])) for sy = (-1,1)])  # handle point on boundary properly

    @test normal([-1,2],s) == normalize([-1,2] - [3,4])
    @test bounds(s) == ([-2,-1],[8,9])
    @test checkbounds(s)
    @test checkbounds(Sphere([1,2,3], 2))

    @test (∆ = rand(2); translate(s,∆) ≈ Sphere([3,4]+∆, 5))
end
