@testset "triangle" begin
    r = 2  # radius
    b = Polygon{3}([0,0], r)  # regular triangle inscribed in circle of radius r (triangel side = √3r/2)

    va = [0,r]  # apex
    vl = [-√3r/2,-r/2]  # left vertex
    vr = [√3r/2,-r/2]  # right vertex
    vs = (va, vl, vr)

    hb = [0,-r/2]  # midpoint of base
    hr = [√3r/4,r/4]  # midpoint of right side
    hl = [-√3r/4,r/4]  # midpoint of left side

    @test b == deepcopy(b)
    @test hash(b) == hash(deepcopy(b))
    @test [0,0] ∈ b  # center
    @test (s = one⁻; s.*va ∈ b && s.*vl ∈ b && s.*vr ∈ b)  # inside apex, left, right
    @test (s = one⁺; s.*va ∉ b && s.*vl ∉ b && s.*vr ∉ b)  # outside apex, left, right

    @test (xn = surfpt_nearby([0,0],b); xn≈(hb,normalize(hb)) || xn≈(hr,normalize(hr)) || xn≈(hl,normalize(hl)))  # handle point at center properly

    @test all([(p₀ = hb; p = s.*p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p₀))) for s = (one⁻⁻,1,one⁺⁺)])  # around base
    @test all([(p₀ = hr; p = s.*p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p₀))) for s = (one⁻⁻,1,one⁺⁺)])  # around right side
    @test all([(p₀ = hl; p = s.*p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p₀))) for s = (one⁻⁻,1,one⁺⁺)])  # around left side

    @test all([(p = one⁺⁺.*p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p₀))) for p₀ = vs])  # outside corners
    @test all([(p = one⁺⁺.*[p₀[1]*0.98, p₀[2]*1.02]; surfpt_nearby(p,b) ≈ (p₀, normalize(p-p₀))) for p₀ = vs])  # outside corners
    @test all([(p = p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p₀))) for p₀ = vs])  # at corners

    @test bounds(b) ≈ (vl, [vr[1],va[2]])
    @test checkbounds(b)

    @test (∆ = rand(2); translate(b,∆) ≈ Polygon{3}(∆,r,0.0))

    @test b ≈ Polygon([va vl vr])
    @test Polygon([va vl vr]) == Polygon([va vr vl]) == Polygon([vl va vr]) == Polygon([vl vr va]) == Polygon([vr va vl]) == Polygon([vr vl va])
    @test b ≈ Isosceles((vl,vr), 1.5r)
end  # @testset "triangle"
