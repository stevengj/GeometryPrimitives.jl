@testset "half disk" begin
    ϕ = π/6  # start angle
    θ = ϕ-π/2  # outward normal to diameter
    c = 0.01 .* [cos(θ), sin(θ)]  # center of disk
    r = 2  # radius
    ∆ϕ = π  # half disk
    b = Sector(c, r, ϕ, ∆ϕ)  # note sector includes origin

    ∆vₛ = r * [cos(ϕ), sin(ϕ)]  # start of arc
    ∆vₑ = r * [cos(ϕ+∆ϕ), sin(ϕ+∆ϕ)]  # end of arc
    ∆vs = (∆vₛ, ∆vₑ)

    @test b == deepcopy(b)
    @test b ≈ Sector(c, r, ϕ+2π, ∆ϕ)
    @test b ≈ Sector(c, r, ϕ-2π, ∆ϕ)
    @test b ≈ Sector(c, r, ϕ+∆ϕ, -∆ϕ)
    @test b ≈ Sector(c, r, ϕ+∆ϕ-2π, -∆ϕ)
    @test hash(b) == hash(deepcopy(b))
    @test [0,0] ∈ b  # origin
    @test c ∈ b  # center
    @test (s = one⁻; s.*(∆vₛ+c) ∈ b && s.*(∆vₑ+c) ∈ b)  # slightly inside of end points of arc
    @test (s = one⁺; s.*(∆vₑ+c) ∉ b && s.*(∆vₑ+c) ∉ b)  # slightly outside of end points of arc

    @test all([(xn = surfpt_nearby(s.*c,b); xn ≈ (c, [cos(θ), sin(θ)])) for s = (one⁻⁻, 1, one⁺⁺)])  # handle point at center properly

    @test all([(p₀ = ∆vₛ/2 + c; p = s.*p₀; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1]) ≈ ϕ && n ≈ [cos(θ), sin(θ)]) for s = (one⁻⁻,1,one⁺⁺)])  # around start side
    @test all([(p₀ = ∆vₑ/2 + c; p = s.*p₀; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1]) ≈ ϕ+∆ϕ-2π && n ≈ [cos(θ), sin(θ)]) for s = (one⁻⁻,1,one⁺⁺)])  # around end side
    @test all([(α = ϕ+0.7∆ϕ; d̂ = [cos(α),sin(α)]; p = s.*r.*d̂ + c; surfpt_nearby(p,b) ≈ (r.*d̂.+c, d̂)) for s = (one⁻⁻,1,one⁺⁺)])  # around arc

    @test all([(p₀ = ∆v + c; p = one⁺⁺.*p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p-p₀))) for ∆v = ∆vs])  # outside end points of arc
    @test (p = ∆vₛ + c; surfpt_nearby(p,b) ≈ (∆vₛ+c, normalize([cos(ϕ),sin(ϕ)]+[cos(θ),sin(θ)])))   # at start point of arc
    @test (p = ∆vₑ + c; surfpt_nearby(p,b) ≈ (∆vₑ+c, normalize([cos(ϕ+∆ϕ),sin(ϕ+∆ϕ)]+[cos(θ),sin(θ)])))   # at end point of arc

    @test bounds(b) ≈ ([c[1]-r, c[2]+r*sin(ϕ+∆ϕ)], [c[1]+r*cos(ϕ), c[2]+r])
    @test checkbounds(b)

    @test (∆ = rand(2); translate(b,∆) ≈ Sector(c+∆, r, ϕ, ∆ϕ))
end  # @testset "half disk"

@testset "Sector, acute angle" begin
    ϕ = π/6  # start angle
    θ₁ = ϕ-π/2  # outward normal to start edge
    ∆ϕ = π/6  # half disk
    θ₂ = (ϕ+∆ϕ) + π/2  # outward normal to end edge

    ϕ₀ = ϕ + ∆ϕ/2
    c = -0.01 .* [cos(ϕ₀), sin(ϕ₀)]  # center
    r = 2  # radius
    b = Sector(c, r, ϕ, ∆ϕ)  # note sector includes origin

    ∆vₛ = r * [cos(ϕ), sin(ϕ)]  # start of arc
    ∆vₑ = r * [cos(ϕ+∆ϕ), sin(ϕ+∆ϕ)]  # end of arc
    ∆vs = (∆vₛ, ∆vₑ)

    @test b == deepcopy(b)
    @test b ≈ Sector(c, r, ϕ+2π, ∆ϕ)
    @test b ≈ Sector(c, r, ϕ-2π, ∆ϕ)
    @test b ≈ Sector(c, r, ϕ+∆ϕ, -∆ϕ)
    @test b ≈ Sector(c, r, ϕ+∆ϕ-2π, -∆ϕ)
    @test hash(b) == hash(deepcopy(b))
    @test [0,0] ∈ b  # origin
    @test c ∈ b  # center
    @test (s = one⁻; s.*(∆vₛ+c) ∈ b && s.*(∆vₑ+c) ∈ b)  # slightly outside of end points of arc
    @test (s = one⁺; s.*(∆vₑ+c) ∉ b && s.*(∆vₑ+c) ∉ b)  # slightly inside of end points of arc

    @test surfpt_nearby(c,b) ≈ (c, -[cos(ϕ₀), sin(ϕ₀)])  # handle point at center properly
    @test all([(p = c+0.1r.*[cos(θ),sin(θ)]; surfpt_nearby(p,b) ≈ (c, normalize(p-c))) for θ = (ϕ+π, ϕ+∆ϕ+π, ϕ₀+π, ϕ+π+0.1∆ϕ, ϕ+∆ϕ+π-0.1∆ϕ)])  # outside center

    @test all([(p₀ = ∆vₛ/2 + c; p = s.*p₀; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1]) ≈ ϕ && n ≈ [cos(θ₁), sin(θ₁)]) for s = (one⁻⁻,1,one⁺⁺)])  # around start side
    @test all([(p₀ = ∆vₑ/2 + c; p = s.*p₀; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1]) ≈ ϕ+∆ϕ && n ≈ [cos(θ₂), sin(θ₂)]) for s = (one⁻⁻,1,one⁺⁺)])  # around end side
    @test all([(α = ϕ+0.7∆ϕ; d̂ = [cos(α),sin(α)]; p = s.*r.*d̂ + c; surfpt_nearby(p,b) ≈ (r.*d̂.+c, d̂)) for s = (one⁻⁻,1,one⁺⁺)])  # around arc

    @test all([(p₀ = ∆v + c; p = one⁺⁺.*p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p-p₀))) for ∆v = ∆vs])  # outside end points of arc
    @test (p = ∆vₛ + c; surfpt_nearby(p,b) ≈ (∆vₛ+c, normalize([cos(ϕ),sin(ϕ)]+[cos(θ₁),sin(θ₁)])))   # at start point of arc
    @test (p = ∆vₑ + c; surfpt_nearby(p,b) ≈ (∆vₑ+c, normalize([cos(ϕ+∆ϕ),sin(ϕ+∆ϕ)]+[cos(θ₂),sin(θ₂)])))   # at end point of arc


    @test bounds(b) ≈ (c, [c[1]+r*cos(ϕ), c[2]+r*sin(ϕ+∆ϕ)])
    @test checkbounds(b)

    @test (∆ = rand(2); translate(b,∆) ≈ Sector(c+∆, r, ϕ, ∆ϕ))
end  # @testset "Sector, acute angle"

@testset "Half pipe" begin
    ϕ = π/6  # start angle
    θ = ϕ-π/2  # outward normal to diameter
    c = 0.01 .* [cos(θ), sin(θ)]  # center of disk
    r = 2  # radius
    ∆ϕ = π  # half disk
    h2 = 1.1

    b = SectoralPrism([c; 0], r, ϕ, ∆ϕ, 2h2)

    ∆vₛ = r * [cos(ϕ), sin(ϕ)]  # start of arc
    ∆vₑ = r * [cos(ϕ+∆ϕ), sin(ϕ+∆ϕ)]  # end of arc
    ∆vs = (∆vₛ, ∆vₑ)

    @test b == deepcopy(b)
    @test b ≈ SectoralPrism([c; 0], r, ϕ+2π, ∆ϕ, 2h2)
    @test b ≈ SectoralPrism([c; 0], r, ϕ-2π, ∆ϕ, 2h2)
    @test b ≈ SectoralPrism([c; 0], r, ϕ+∆ϕ, -∆ϕ, 2h2)
    @test b ≈ SectoralPrism([c; 0], r, ϕ+∆ϕ-2π, -∆ϕ, 2h2)
    @test hash(b) == hash(deepcopy(b))
    @test [0,0,0] ∈ b  # center
    @test [c; 0] ∈ b  # center
    @test all([(s = one⁻; [s.*(∆v+c);t*h2] ∈ b) for t = (-one⁻,0,one⁻), ∆v = ∆vs])  # slightly inside of end points of arc
    @test all([(s = one⁺; [s.*(∆v+c);t*h2] ∉ b) for t = (-one⁻,0,one⁻), ∆v = ∆vs])  # slightly outside of end points of arc

    @test all([(xn = surfpt_nearby([s.*c;t*h2/2], b); xn ≈ ([c;t*h2/2], [cos(θ),sin(θ),0])) for s = (one⁻⁻, 1, one⁺⁺), t = (-1,0,1)])  # handle point at center properly

    @test all([(p₀ = ∆vₛ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ && x[3]==t*h2 && n ≈ [cos(θ),sin(θ),0]) for s = (one⁻⁻,1,one⁺⁺), t = (-0.5,0,0.5)])  # around start side, within height
    @test all([(s = 1; p₀ = ∆vₛ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ && x[3]==t*h2 && n ≈ normalize([cos(θ),sin(θ),t])) for t = (-1,1)])  # at start side, at base
    @test all([(s = one⁻⁻; p₀ = ∆vₛ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); x==p && n==[0,0,t]) for t = (-1,1)])  # inside start side, at base
    @test all([(s = one⁺⁺; p₀ = ∆vₛ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ && x[3]≈t*h2 && n ≈ normalize([cos(θ),sin(θ),0])) for t = (-1,1)])  # outside start side, at base
    @test all([(s = one⁺⁺; p₀ = ∆vₛ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ && x[3]≈sign(t)*h2 && n ≈ normalize(p-x)) for t = (-one⁺⁺,one⁺⁺)])  # outside start side, outside base

    @test all([(p₀ = ∆vₑ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ+∆ϕ-2π && x[3]==t*h2 && n ≈ [cos(θ),sin(θ),0]) for s = (one⁻⁻,1,one⁺⁺), t = (-0.5,0,0.5)])  # around end side, within height
    @test all([(s = 1; p₀ = ∆vₑ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ+∆ϕ-2π && x[3]==t*h2 && n ≈ normalize([cos(θ),sin(θ),t])) for t = (-1,1)])  # at end side, at base
    @test all([(s = one⁻⁻; p₀ = ∆vₑ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); x==p && n==[0,0,t]) for t = (-1,1)])  # inside end side, at base
    @test all([(s = one⁺⁺; p₀ = ∆vₑ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ+∆ϕ-2π && x[3]≈t*h2 && n ≈ normalize([cos(θ),sin(θ),0])) for t = (-1,1)])  # outside end side, at base
    @test all([(s = one⁺⁺; p₀ = ∆vₑ/2 + c; p = [s.*p₀;t*h2]; (x,n) = surfpt_nearby(p,b); atan(x[2]-c[2],x[1]-c[1])≈ϕ+∆ϕ-2π && x[3]≈sign(t)*h2 && n ≈ normalize(p-x)) for t = (-one⁺⁺,one⁺⁺)])  # outside end side, outside base

    @test all([(α = ϕ+0.7∆ϕ; d̂ = [cos(α),sin(α)]; p = [s.*r.*d̂ + c; t*h2]; surfpt_nearby(p,b) ≈ ([r.*d̂.+c; t*h2], [d̂;0])) for s = (one⁻⁻,1,one⁺⁺), t = (-0.5,0,0.5)])  # around arc, within height
    @test all([(s = 1; α = ϕ+0.7∆ϕ; d̂ = [cos(α),sin(α)]; p = [s.*r.*d̂ + c; t*h2]; surfpt_nearby(p,b) ≈ ([r.*d̂.+c; t*h2], normalize([d̂;t]))) for t = (-1,1)])  # at arc, at base
    @test all([(s = one⁻⁻; α = ϕ+0.7∆ϕ; d̂ = [cos(α),sin(α)]; p = [s.*r.*d̂ + c; t*h2]; surfpt_nearby(p,b) ≈ (p, [0,0,t])) for t = (-1,1)])  # inside arc, at base
    @test all([(s = one⁺⁺; α = ϕ+0.7∆ϕ; d̂ = [cos(α),sin(α)]; p = [s.*r.*d̂ + c; t*h2]; surfpt_nearby(p,b) ≈ ([r.*d̂.+c; t*h2], normalize([d̂;0]))) for t = (-1,1)])  # outside arc, at base
    @test all([(s = one⁺⁺; α = ϕ+0.7∆ϕ; d̂ = [cos(α),sin(α)]; p = [s.*r.*d̂ + c; t*h2]; (x,n) = surfpt_nearby(p,b); (x,n) ≈ ([r.*d̂.+c; sign(t)*h2], normalize(p-x))) for t = (-one⁺⁺,one⁺⁺)])  # outside arc, outside base

    @test all([(p = [∆vₛ + c; t*h2]; surfpt_nearby(p,b) ≈ (p, normalize(normalize([cos(ϕ),sin(ϕ),0]+[cos(θ),sin(θ),0])+[0,0,t]))) for t = (-1,1)])  # at start point of arc
    @test all([(p = [∆vₑ + c; t*h2]; surfpt_nearby(p,b) ≈ (p, normalize(normalize([cos(ϕ+∆ϕ),sin(ϕ+∆ϕ),0]+[cos(θ),sin(θ),0])+[0,0,t]))) for t = (-1,1)])  # at start point of arc
    @test all([(p₀ = [∆v + c; t*h2]; p = one⁺⁺.*p₀; surfpt_nearby(p,b) ≈ (p₀, normalize(p-p₀))) for ∆v = ∆vs, t = (-1,1)])  # outside end points of arc

    @test bounds(b) ≈ ([c[1]-r, c[2]+r*sin(ϕ+∆ϕ), -h2], [c[1]+r*cos(ϕ), c[2]+r, h2])
    @test checkbounds(b)

    @test (∆ = rand(3); translate(b,∆) ≈ SectoralPrism([c;0]+∆, r, ϕ, ∆ϕ, 2h2))
end  # @testset "Half pipe"
