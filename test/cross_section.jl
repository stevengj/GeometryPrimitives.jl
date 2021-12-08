@testset "CrossSection, ball" begin
    r = 2
    b = Ball([0,0,0], r)
    n = rand(3)
    s = CrossSection(b, n, 0)

    @test s == deepcopy(s)
    @test hash(s) == hash(deepcopy(s))
    @test ndims(s) == 2

    @test all(([r,0], [-r,0], [0,r], [0,-r], [0,0]) .∈ Ref(s))
    @test !any(([r*one⁺,0], [-r*one⁺,0], [0,r*one⁺], [0,-r*one⁺]) .∈ Ref(s))

    ∆ = rand(2)
    s2 = translate(s, ∆)
    @test all(([r,0], [-r,0], [0,r], [0,-r], [0,0]) .+ Ref(∆) .∈ Ref(s2))
    @test !any(([r*one⁺,0], [-r*one⁺,0], [0,r*one⁺], [0,-r*one⁺]) .+ Ref(∆) .∈ Ref(s2))

end  # @testset "CrossSection, ball"
