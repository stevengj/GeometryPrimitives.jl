@testset "triangular prism" begin
    r = 1
    h2 = 1.1
    c = Prism([0,0,0], Polygon{3}([0,0],r), 2h2)

    va = [0,r]  # apex
    vl = [-√3r/2,-r/2]  # left vertex
    vr = [√3r/2,-r/2]  # right vertex
    vs = (va, vl, vr)

    hb = [0,-r/2]  # midpoint of base line of base triangle
    hr = [√3r/4,r/4]  # midpoint of right side of base triangle
    hl = [-√3r/4,r/4]  # midpoint of left side of base triangle
    hs = (hb, hr, hl)

    @test c == deepcopy(c)
    @test hash(c) == hash(deepcopy(c))
    @test [0,0,0] ∈ c  # center
    @test all([(s = one⁻; [s.*va;t*h2] ∈ c && [s.*vl;t*h2] ∈ c && [s.*vr;t*h2] ∈ c) for t = (-one⁻,0,one⁻)]) # inside apex, left, right
    @test all([(s = one⁺; [s.*va;t*h2] ∉ c && [s.*vl;t*h2] ∉ c && [s.*vr;t*h2] ∉ c) for t = (-one⁺,0,one⁺)]) # inside apex, left, right

    @test (xn = surfpt_nearby([0,0,0],c); xn≈([hb;0],normalize([hb;0])) || xn≈([hr;0],normalize([hr;0])) || xn≈([hl;0],normalize([hl;0])))  # handle point at center properly

    @test all([(p = one⁺⁺.*[p₀;q₀]; surfpt_nearby(p,c) ≈ ([p₀;q₀], normalize(p - [p₀;q₀]))) for p₀ = vs, q₀ = (-h2,h2)])  # outside corners

    # Below, s = one⁻ is excluded because it is an inside point, for which nout is calculated differently.
    @test all([(p = s.*[p₀;q₀]; surfpt_nearby(p,c) ≈ ([p₀;q₀], normalize([normalize(p₀);q₀/abs(q₀)]))) for p₀ = vs, q₀ = (-h2,h2), s = (1,one⁺)])  # around corners
    @test all([(p = s.*[p₀;q₀]; surfpt_nearby(p,c) ≈ ([p₀;q₀], normalize([normalize(p₀);0]))) for p₀ = vs, q₀ = (-h2/2,0,h2/2), s = (1,one⁺)])  # around vertical edges

    @test all([(p = [0.5.*p₀;s*q₀]; surfpt_nearby(p,c) ≈ ([0.5.*p₀;q₀], [0,0,sign(q₀)])) for p₀ = vs, q₀ = (-h2,h2), s = (one⁻⁻,one⁻,1,one⁺,one⁺⁺)])  # around bases
    @test all([(p = [s.*p₀;q₀]; surfpt_nearby(p,c) ≈ ([p₀;q₀], [normalize(p₀);0])) for p₀ = hs, q₀ = (-h2/2,0,h2/2), s = (one⁻⁻,one⁻,1,one⁺,one⁺⁺)])  # around sides

    @test bounds(c) ≈ ([vl;-h2],[vr[1],va[2],h2])
    @test checkbounds(c)

    @test (∆ = rand(3); translate(c,∆) ≈ Prism([0,0,0].+∆, Polygon{3}([0,0],r), 2h2))
end  # @testset "triangular prism"


@testset "triangular prism, rotated" begin
    ax1 = normalize([1,0,-1])
    ax2 = [0,1,0]
    ax3 = normalize([1,0,1])  # ax3 = ax1 × ax2
    ax = [ax1 ax2 ax3]

    r = 1
    h2 = 1.1
    c = Prism([0,0,0], Polygon{3}([0,0],r), 2h2, ax)

    va = [0,r]  # apex
    vl = [-√3r/2,-r/2]  # left vertex
    vr = [√3r/2,-r/2]  # right vertex
    vs = (va, vl, vr)

    hb = [0,-r/2]  # midpoint of base line of base triangle
    hr = [√3r/4,r/4]  # midpoint of right side of base triangle
    hl = [-√3r/4,r/4]  # midpoint of left side of base triangle
    hs = (hb, hr, hl)

    @test c == deepcopy(c)
    @test hash(c) == hash(deepcopy(c))
    @test [0,0,0] ∈ c  # center
    @test all([(s = one⁻; ax*[s.*va;t*h2] ∈ c && ax*[s.*vl;t*h2] ∈ c && ax*[s.*vr;t*h2] ∈ c) for t = (-one⁻,0,one⁻)]) # inside apex, left, right
    @test all([(s = one⁺; ax*[s.*va;t*h2] ∉ c && ax*[s.*vl;t*h2] ∉ c && ax*[s.*vr;t*h2] ∉ c) for t = (-one⁺,0,one⁺)]) # inside apex, left, right

    @test (xn = surfpt_nearby([0,0,0],c); xn≈(ax*[hb;0],ax*normalize([hb;0])) || xn≈(ax*[hr;0],ax*normalize([hr;0])) || xn≈(ax*[hl;0],ax*normalize([hl;0])))  # handle point at center properly

    @test all([(p = one⁺⁺.*[p₀;q₀]; surfpt_nearby(ax*p,c) ≈ (ax*[p₀;q₀], ax*normalize(p - [p₀;q₀]))) for p₀ = vs, q₀ = (-h2,h2)])  # outside corners

    # Below, s = one⁻ is excluded because it is an inside point, for which nout is calculated differently.
    @test all([(p = s.*[p₀;q₀]; surfpt_nearby(ax*p,c) ≈ (ax*[p₀;q₀], ax*normalize([normalize(p₀);q₀/abs(q₀)]))) for p₀ = vs, q₀ = (-h2,h2), s = (1,one⁺)])  # around corners
    @test all([(p = s.*[p₀;q₀]; surfpt_nearby(ax*p,c) ≈ (ax*[p₀;q₀], ax*normalize([normalize(p₀);0]))) for p₀ = vs, q₀ = (-h2/2,0,h2/2), s = (1,one⁺)])  # around vertical edges

    @test all([(p = [0.5.*p₀;s*q₀]; surfpt_nearby(ax*p,c) ≈ (ax*[0.5.*p₀;q₀], ax*[0,0,sign(q₀)])) for p₀ = vs, q₀ = (-h2,h2), s = (one⁻⁻,one⁻,1,one⁺,one⁺⁺)])  # around bases
    @test all([(p = [s.*p₀;q₀]; surfpt_nearby(ax*p,c) ≈ (ax*[p₀;q₀], ax*[normalize(p₀);0])) for p₀ = hs, q₀ = (-h2/2,0,h2/2), s = (one⁻⁻,one⁻,1,one⁺,one⁺⁺)])  # around sides

    @test (vtot = [[vl;h2] [vl;-h2] [vr;h2] [vr;-h2] [va;h2] [va;-h2]]; bounds(c) ≈ (minimum(ax*vtot, dims=2)[:,1], maximum(ax*vtot, dims=2)[:,1]))
    @test checkbounds(c)

    @test (∆ = rand(3); translate(c,∆) ≈ Prism([0,0,0].+∆, Polygon{3}([0,0],r), 2h2, ax))
end  # @testset "triangular prism, rotated"
