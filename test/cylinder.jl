@testset "Cylinder" begin
    c = Cylinder([0,0,0], 0.3, 2.2, [0,0,1])
    @test c == deepcopy(c)
    @test hash(c) == hash(deepcopy(c))
    @test [0.2,0.2,1] ∈ c
    @test SVector(0.2,0.2,1.2) ∉ c
    @test [0.2,0.25,1] ∉ c

    @test ((x,nout) = surfpt_nearby([0,0,0],c); (x≈[0,0,1.1] && nout≈[0,0,1]) || (x≈[0.3,0,0] && nout≈[1,0,0]) || (x≈[0,0.3,0] && nout≈[0,1,0]))  # handle point at center properly
    @test all([(p = [0,0,1.1sz]; surfpt_nearby(p,c) ≈ (p,[0,0,sz])) for sz = (-1,1)])  # handle point on axis properly
    @test all([(p = [0.3sx/√2,0.3sy/√2,1.1sz]; surfpt_nearby(1.1p,c) ≈ (p, normalize(1.1p-p))) for sx = (-1,1), sy = (-1,1), sz = (-1,1)])  # outside corners
    @test all([(p = [0.3sx,0,1.1sz]; surfpt_nearby(1.1p,c) ≈ (p, normalize(1.1p-p))) for sx = (-1,1), sz = (-1,1)])  # outside corners
    @test all([(p = [0,0.3sy,1.1sz]; surfpt_nearby(1.1p,c) ≈ (p, normalize(1.1p-p))) for sy = (-1,1), sz = (-1,1)])  # outside corners
    @test all([(p = [0.3sx/√2,0.3sy/√2,1.1sz]; (x,nout) = surfpt_nearby(p,c); (x≈p && all([sx 0 0; 0 sy 0; 0 0 sz]*nout.≥-10eps()))) for sx = (-1,1), sy = (-1,1), sz = (-1,1)])  # on rims
    @test all([(p = [0.3sx,0,1.1sz]; (x,nout) = surfpt_nearby(p,c); (x≈p && all([sx 0 0; 0 0 sz]*nout.≥-10eps()))) for sx = (-1,1), sz = (-1,1)])  # on rims
    @test all([(p = [0,0.3sy,1.1sz]; (x,nout) = surfpt_nearby(p,c); (x≈p && all([0 sy 0; 0 0 sz]*nout.≥-10eps()))) for sy = (-1,1), sz = (-1,1)])  # on rims
    @test all([(p = [0.3sx/√2/2,0.3sy/√2/2,1.1sz]; surfpt_nearby([p[1],p[2],ρ*p[3]],c) ≈ (p,[0,0,sz])) for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,0,1), sy = (-1,0,1), sz = (-1,1)])  # around bases
    @test all([(p = [0.3sx,0,1.1sz/2]; surfpt_nearby([ρ*p[1],p[2],p[3]],c) ≈ (p,[sx,0,0])) for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,1), sz = (-1,0,1)])  # around side
    @test all([(p = [0,0.3sy,1.1sz/2]; surfpt_nearby([p[1],ρ*p[2],p[3]],c) ≈ (p,[0,sy,0])) for ρ = (one⁻⁻,1,one⁺⁺), sy = (-1,1), sz = (-1,0,1)])  # around side

    @test @inferred(normal([0.1,0.2,0.3], c)) ≈ [0.1,0.2,0] / hypot(0.1,0.2)
    @test normal([0.1,0.2,-1.11], c) == [0,0,-1]
    @test normal([0.31, 0, 0.3], c) == [1,0,0]
    @test bounds(c) ≈ ([-0.3,-0.3,-1.1],[0.3,0.3,1.1])
    @test checkbounds(c)
    @test checkbounds(Cylinder([1,17,44], 0.3, 1.1, [1,-2,3]))

    @test (∆ = rand(3); translate(c,∆) ≈ Cylinder([0,0,0]+∆, 0.3, 2.2, [0,0,1]))
end  # @testset "Cylinder"


@testset "Cylinder, rotated" begin
    ax1 = normalize([1,0,-1])
    ax2 = [0,1,0]
    ax3 = normalize([1,0,1])  # ax3 = ax1 × ax2
    cr = Cylinder([0,0,0], 0.3, 2.2, ax3)
    @test cr == deepcopy(cr)
    @test hash(cr) == hash(deepcopy(cr))

    @test ((x,nout) = surfpt_nearby([0,0,0],cr); (x≈0.3ax1 && nout≈ax1) || (x≈0.3ax2 && nout≈ax2) || (x≈1.1ax3 && nout≈ax3))  # handle point at center properly
    @test all([(p = 1.1s3*ax3; surfpt_nearby(p,cr) ≈ (p,s3*ax3)) for s3 = (-1,1)])  # handle point on axis properly
    @test all([(p = 0.3(s1*ax1+s2*ax2)/√2+1.1s3*ax3; surfpt_nearby(1.1p,cr) ≈ (p, normalize(1.1p-p))) for s1 = (-1,1), s2 = (-1,1), s3 = (-1,1)])  # outside corners
    @test all([(p = 0.3s1*ax1+1.1s3*ax3; surfpt_nearby(1.1p,cr) ≈ (p, normalize(1.1p-p))) for s1 = (-1,1), s3 = (-1,1)])  # outside corners
    @test all([(p = 0.3s2*ax2+1.1s3*ax3; surfpt_nearby(1.1p,cr) ≈ (p, normalize(1.1p-p))) for s2 = (-1,1), s3 = (-1,1)])  # outside corners
    @test all([(p = 0.3(s1*ax1+s2*ax2)/√2+1.1s3*ax3; (x,nout) = surfpt_nearby(p,cr); (x≈p && all([s1*ax1 s2*ax2 s3*ax3]'*nout.≥-10eps()) && norm(nout)≈1)) for s1 = (-1,1), s2 = (-1,1), s3 = (-1,1)])  # on rims
    @test all([(p = 0.3s1*ax1+1.1s3*ax3; (x,nout) = surfpt_nearby(p,cr); (x≈p && all([s1*ax1 s3*ax3]'*nout.≥-10eps()) && norm(nout)≈1)) for s1 = (-1,1), s3 = (-1,1)])  # on rims
    @test all([(p = 0.3s2*ax2+1.1s3*ax3; (x,nout) = surfpt_nearby(p,cr); (x≈p && all([s2*ax2 s3*ax3]'*nout.≥-10eps()) && norm(nout)≈1)) for s2 = (-1,1), s3 = (-1,1)])  # on rims
    @test all([surfpt_nearby(0.3(s1*ax1+s2*ax2)/√2/2+ρ*1.1s3*ax3,cr) ≈ (0.3(s1*ax1+s2*ax2)/√2/2+1.1s3*ax3,s3*ax3) for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,0,1), s2 = (-1,0,1), s3 = (-1,1)])  # around bases
    @test all([surfpt_nearby(ρ*0.3s1*ax1+1.1s3*ax3/2,cr) ≈ (0.3s1*ax1+1.1s3*ax3/2, s1*ax1) for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,1), s3 = (-1,0,1)])  # around side
    @test all([surfpt_nearby(ρ*0.3s2*ax2+1.1s3*ax3/2,cr) ≈ (0.3s2*ax2+1.1s3*ax3/2, s2*ax2) for ρ = (one⁻⁻,1,one⁺⁺), s2 = (-1,1), s3 = (-1,0,1)])  # around side

    @test all([@inferred(normal(0.3(s1*ax1+s2*ax2)/√2/2+ρ*1.1s3*ax3,cr)) ≈ s3*ax3 for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,0,1), s2 = (-1,0,1), s3 = (-1,1)])  # around bases
    @test all([normal(ρ*0.3s2*ax2+1.1s3/2*ax3,cr) ≈ s2*ax2 for ρ = (one⁻⁻,1,one⁺⁺), s2 = (-1,1), s3 = (-1,0,1)])  # around side
    @test all([normal(ρ*0.3s1*ax1+1.1s3/2*ax3,cr) ≈ s1*ax1 for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,1), s3 = (-1,0,1)])  # around side

    @test bounds(cr) ≈ (-[(1.1+0.3)/√2,0.3,(1.1+0.3)/√2], [(1.1+0.3)/√2,0.3,(1.1+0.3)/√2])
    @test checkbounds(cr)

    @test (∆ = rand(3); translate(cr,∆) ≈ Cylinder([0,0,0]+∆, 0.3, 2.2, ax3))
end  # @testset "Cylinder, rotated"
