@testset "Ellipsoid" begin
    e = Ellipsoid([0,0], [1,2])
    @test e == deepcopy(e)
    @test hash(e) == hash(deepcopy(e))
    @test [0.3,2*sqrt(1 - 0.3^2)-0.01] ∈ e
    @test [0.3,2*sqrt(1 - 0.3^2)+0.01] ∉ e

    @test ((x,nout) = surfpt_nearby([0,0],e); (x≈[1,0] && nout≈[1,0]) || (x≈[0,2] && nout≈[0,1]))  # handle point at center properly
    @test all([surfpt_nearby([ρ*1sx,0],e)[1] ≈ [1sx,0] for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])
    @test all([surfpt_nearby([0,ρ*2sy],e)[1] ≈ [0,2sy] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
    @test all([(p = [1sx,0]; surfpt_nearby(p,e) ≈ (p,[sx,0])) for sx = (-1,1)])  # handle point on boundary properly
    @test all([(p = [0,2sy]; surfpt_nearby(p,e) ≈ (p,[0,sy])) for sy = (-1,1)])  # handle point on boundary properly

    @test normal([1.1,0],e) ≈ [1,0]
    @test normal([-1.1,0],e) ≈ [-1,0]
    @test normal([0,2.01],e) ≈ [0,1]
    @test bounds(e) == ([-1,-2],[1,2])
    @test checkbounds(e)
    @test checkbounds(Ellipsoid([0,0], [1,2], [1 1; 1 -1]))

    @test (∆ = rand(2); translate(e,∆) ≈ Ellipsoid([0,0]+∆, [1,2]))

    b = Box([0,0], [2,4])
    eb = Ellipsoid(b)
    @test e == eb
    @test hash(e) == hash(eb)
end  # @testset "Ellipsoid"


@testset "Ellipsoid, rotated" begin
    θ = π/3
    R = [cos(θ) sin(θ); sin(θ) -cos(θ)]
    er = Ellipsoid([0,0], [2,3], R)
    bp = GeometryPrimitives.boundpts(er)

    bp1, bp2 = bp[:,1], bp[:,2]

    # Test the two bounding points are on the ellipsoid perimeter.
    @test er == deepcopy(er)
    @test hash(er) == hash(deepcopy(er))
    @test (one⁻ * bp1 ∈ er) && (one⁻ * bp2 ∈ er)
    @test (one⁺ * bp1 ∉ er) && (one⁺ * bp2 ∉ er)

    @test ((x,nout) = surfpt_nearby([0,0],er); (x≈2*R[:,1] && nout≈R[:,1]) || (x≈3*R[:,2] && nout≈R[:,2]))  # handle point at center properly
    @test all([surfpt_nearby(R*[ρ*2sx,0],er)[1] ≈ R*[2sx,0] for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])
    @test all([surfpt_nearby(R*[0,ρ*3sy],er)[1] ≈ R*[0,3sy] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
    @test all([(p = R*[2sx,0]; surfpt_nearby(p,er) ≈ (p,R*[sx,0])) for sx = (-1,1)])  # handle point on boundary properly
    @test all([(p = R*[0,3sy]; surfpt_nearby(p,er) ≈ (p,R*[0,sy])) for sy = (-1,1)])  # handle point on boundary properly

    # Test the normal vector at the two bounding points are the x- and y-directions.
    @test normal(bp1, er) ≈ [1,0]
    @test normal(bp2, er) ≈ [0,1]

    xmax, ymax = bp1[1], bp2[2]
    @test bounds(er) == ([-xmax, -ymax], [xmax, ymax])
    @test checkbounds(er)

    @test (∆ = rand(2); translate(er,∆) ≈ Ellipsoid([0,0]+∆, [2,3], R))

    br = Box([0,0], 2*[2,3], R)
    ebr = Ellipsoid(br)
    @test er == ebr
    @test hash(er) == hash(ebr)
end  # @testset "Ellipsoid, rotated"
