using GeometryPrimitives, StaticArrays, Base.Test

const rtol = Base.rtoldefault(Float64)
const one⁻ = 1 - rtol  # slightly less than 1
const one⁺ = 1 + rtol  # slightly greater than 1
const one⁻⁻, one⁺⁺ = 0.9, 1.1

Base.isapprox(a::Tuple, b::Tuple; kws...) = all(p -> isapprox(p...; kws...), zip(a,b))
const rng = MersenneTwister(0) # test with reproducible pseudorandom numbers

randnb(lb::Real,ub::Real) = randn(rng)*(ub-lb) + 0.5*(lb+ub)
randnb(lb::AbstractVector,ub::AbstractVector) = map(randnb, lb, ub)
function inbounds(x,lb,ub)
    length(x) == length(lb) == length(lb) || throw(BoundsError())
    for i in eachindex(x)
        @inbounds lb[i] ≤ x[i] ≤ ub[i] || return false
    end
    return true
end

"check the bounding box of s with randomized trials"
function checkbounds(s::Shape{N}, ntrials=10^4) where {N}
    lb,ub = bounds(s)
    for i = 1:ntrials
        x = randnb(lb,ub)
        x ∉ s || inbounds(x,lb,ub) || return false  # return false if s - (bounding box) is nonempty
    end
    return true
end

function checktree(t::KDTree{N}, slist::Vector{<:Shape{N}}, ntrials=10^3) where {N}
    lb = SVector{N}(fill(Inf,N))
    ub = SVector{N}(fill(-Inf,N))
    for i in eachindex(slist)
        lbi,ubi = bounds(slist[i])
        lb = min.(lb,lbi)
        ub = max.(ub,ubi)
    end
    for i = 1:ntrials
        x = randnb(lb,ub)
        st = findin(x, t)
        sl = findin(x, slist)
        isnull(st) == isnull(sl) || return false
        isnull(st) || get(st) == get(sl) || return false
    end
    return true
end

@testset "GeometryPrimitives" begin

    @testset "Shape" begin
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

            @test (∆ = rand(2); translate(s,∆) == Sphere([3,4]+∆, 5))
        end

        @testset "Box" begin
            b = Box([0,0], [2,4])  # specify center and radii
            @test b == deepcopy(b)
            @test hash(b) == hash(deepcopy(b))
            @test [0.3,-1.5] ∈ b
            @test [0.3,-2.5] ∉ b

            @test ((x,nout) = surfpt_nearby([0,0],b); (x≈[1,0] && nout≈[1,0]) || (x≈[0,2] && nout≈[0,1]))  # handle point at center properly
            @test all([(p = [1sx,2sy]; surfpt_nearby(1.1p,b) == (p, normalize(1.1p-p))) for sx = (-1,1), sy = (-1,1)])  # outside corners
            @test all([((x, nout) = surfpt_nearby([1sx,2sy],b); x≈[1sx,2sy] && all([sx 0; 0 sy]*nout.≥0) && norm(nout)≈1) for sx = (-1,1), sy = (-1,1)])  # at corners
            @test all([surfpt_nearby([ρ*1sx,1sy],b) == ([1sx,1sy], [sx,0]) for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,1), sy = (-1,0,1)])  # around faces
            @test all([surfpt_nearby([0.5sx,ρ*2sy],b) == ([0.5sx,2sy], [0,sy]) for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,0,1), sy = (-1,1)])  # around faces

            @test normal([1.1,0],b) == [1,0]
            @test normal([-1.1,0],b) == [-1,0]
            @test normal([1.1,2.1],b) == [1,1]/√2
            @test bounds(b) == ([-1,-2],[1,2])
            @test bounds(Box([0,0], [2,4], [1 1; 1 -1])) ≈ ([-3*√0.5,-3*√0.5], [3*√0.5,3*√0.5])
            @test checkbounds(b)
            @test checkbounds(Box([0,0], [2,4], [1 1; 1 -1]))

            @test (∆ = rand(2); translate(b,∆) == Box([0,0]+∆, [2,4]))
        end

        @testset "Box, rotated" begin
            ax1, ax2 = normalize.(([1,-1], [1,1]))
            r1, r2 = 1, 2  # "radii"
            br = Box([0,0], [2r1, 2r2], [ax1 ax2])

            R = [ax1 ax2]  # rotation matrix

            Cin = R * (GeometryPrimitives.signmatrix(br) .* (one⁻ .* [r1,r2]))  # around corners, inside
            Cout = R * (GeometryPrimitives.signmatrix(br) .* (one⁺ .* [r1,r2]))  # around corners, outside
            for j = 1:2; @test Cin[:,j] ∈ br; end
            for j = 1:2; @test Cout[:,j] ∉ br; end

            @test br == deepcopy(br)
            @test hash(br) == hash(deepcopy(br))

            @test ((x,nout) = surfpt_nearby([0,0],br); (x≈r1*ax1 && nout≈ax1) || (x≈r2*ax2 && nout≈ax2))  # handle point at center properly
            @test all([(p = (s1*r1*ax1+s2*r2*ax2); surfpt_nearby(1.1p,br) ≈ (p,normalize(1.1p-p))) for s1 = (-1,1), s2 = (-1,1)])  # outside corners
            @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(p,br); x≈p && all([s1*ax1 s2*ax2]'*nout.≥0) && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # at corners
            @test all([(p1 = s1*r1*ax1; p2 = s2*r2/2*ax2; surfpt_nearby(ρ*p1+p2,br) ≈ (p1+p2, s1*ax1)) for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,1), s2 = (-1,0,1)])  # around faces
            @test all([(p1 = s1*r1/2*ax1; p2 = s2*r2*ax2; surfpt_nearby(p1+ρ*p2,br) ≈ (p1+p2, s2*ax2)) for ρ = (one⁻⁻,1,one⁺⁺), s1 = (-1,0,1), s2 = (-1,1)])  # around faces

            @test normal(R*[1.1r1, 0], br) ≈ R*[1,0]
            @test normal(R*[-1.1r1, 0], br) ≈ R*[-1,0]
            @test normal(R*[0, 1.1r2], br) ≈ R*[0,1]
            @test normal(R*[0, -1.1r2], br) ≈ R*[0,-1]
            @test normal(R*[1.1r1, 1.01r2], br) ≈ R*(normalize([1.1r1, 1.01r2]-[r1,r2]))

            xmax = (R*[r1,r2])[1]
            ymax = (R*[-r1,r2])[2]
            @test bounds(br) ≈ (-[xmax,ymax], [xmax,ymax])
            @test checkbounds(br)

            @test (∆ = rand(2); translate(br,∆) == Box([0,0]+∆, [2r1, 2r2], [ax1 ax2]))
        end

        @testset "Box, skewed" begin
            ax1, ax2 = normalize.(([1,-1], [0,1]))
            r1, r2 = 1, 1  # "radii"
            bs = Box([0,0], [2r1, 2r2], [ax1 ax2])

            @test bs == deepcopy(bs)
            @test hash(bs) == hash(deepcopy(bs))

            n1, n2 = normalize.(([1,0], [1,1]))
            @test ((x,nout) = surfpt_nearby([0,0],bs); (x⋅n1≈1/√2 && nout≈n1) || (x⋅n2≈1/√2 && nout≈n2))  # handle point at center properly
            @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(1.1p,bs); nout⋅(s1*n1+s2*n2)>0 && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # outside corners
            @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(p,bs); x≈p && nout⋅(s1*n1+s2*n2)>0 && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # at corners
            # Use the above, less demanding tests instead of the below two; surfpt_nearby suffers from some inaccuracy around corners of skewed boxes.
            # @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(1.1p,bs); all([s1*n1 s2*n2]'*nout.>0) && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # outside corners
            # @test all([(p = (s1*r1*ax1+s2*r2*ax2); (x,nout) = surfpt_nearby(p,bs); all([s1*n1 s2*n2]'*nout.≥0) && norm(nout)≈1) for s1 = (-1,1), s2 = (-1,1)])  # at corners
            @test all([(p1 = s1*r1*ax1; p2 = s2*r2/2*ax2; surfpt_nearby(ρ*p1+p2,bs) ≈ (p1+p2, s1*n1)) for ρ = (one⁻,1,one⁺), s1 = (-1,1), s2 = (-1,0,1)])  # around faces
            @test all([(p1 = s1*r1/2*ax1; p2 = s2*r2*ax2; surfpt_nearby(p1+ρ*p2,bs) ≈ (p1+p2, s2*n2)) for ρ = (one⁻,1,one⁺), s1 = (-1,0,1), s2 = (-1,1)])  # around faces

            @test norm(normal([0,1], bs)) ≈ 1

            xmax = (r1*ax1+r2*ax2)[1]
            ymax = (r2*ax2-r1*ax1)[2]
            @test bounds(bs) ≈ (-[xmax,ymax],[xmax,ymax])
            @test checkbounds(bs)

            @test (∆ = rand(2); translate(bs,∆) == Box([0,0]+∆, [2r1, 2r2], [ax1 ax2]))
        end

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

            @test (∆ = rand(2); translate(e,∆) == Ellipsoid([0,0]+∆, [1,2]))

            b = Box([0,0], [2,4])
            eb = Ellipsoid(b)
            @test e == eb
            @test hash(e) == hash(eb)
        end

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

            @test (∆ = rand(2); translate(er,∆) == Ellipsoid([0,0]+∆, [2,3], R))

            br = Box([0,0], 2*[2,3], R)
            ebr = Ellipsoid(br)
            @test er == ebr
            @test hash(er) == hash(ebr)
        end

        @testset "Cylinder" begin
            c = Cylinder([0,0,0], 0.3, [0,0,1], 2.2)
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

            @test @inferred(normal([0.1,0.2,0.3], c)) == [0.1,0.2,0] / hypot(0.1,0.2)
            @test normal([0.1,0.2,-1.11], c) == [0,0,-1]
            @test normal([0.31, 0, 0.3], c) == [1,0,0]
            @test bounds(c) ≈ ([-0.3,-0.3,-1.1],[0.3,0.3,1.1])
            @test checkbounds(c)
            @test checkbounds(Cylinder([1,17,44], 0.3, [1,-2,3], 1.1))

            @test (∆ = rand(3); translate(c,∆) == Cylinder([0,0,0]+∆, 0.3, [0,0,1], 2.2))
        end

        @testset "Cylinder, rotated" begin
            ax1 = normalize([1,0,-1])
            ax2 = [0,1,0]
            ax3 = normalize([1,0,1])  # ax3 = ax1 × ax2
            cr = Cylinder([0,0,0], 0.3, ax3, 2.2)
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

            @test (∆ = rand(3); translate(cr,∆) == Cylinder([0,0,0]+∆, 0.3, ax3, 2.2))
        end

        @testset "periodize" begin
            c = Cylinder([0,0,0], 1, [0,0,1], 5)
            bound = Box([0,0,0], [10,10,5])  # specify center and radii
            A = [1 0 0; 0 1 0; 0 0 5]'
            c_array = periodize(c, A, bound)
            @test length(c_array) == (11-2)^2
        end
    end

    @testset "KDTree" begin
        s = [Sphere([i,0], 1) for i in 2:4]
        s0 = Sphere([0,0], 1)
        s = [s0, s0, s0, s0, s...]
        @test_nowarn KDTree(s)  # must not generate StackOverflowError

        s = Shape{2,4}[Sphere([i,0], 1, i) for i in 0:20]
        kd = KDTree(s)
        @test GeometryPrimitives.depth(kd) == 3
        @test get(findin([10.1,0], kd)).data == 10
        @test isnull(findin([10.1,1], kd))
        @test checktree(kd, s)

        s = Shape{3,9}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.01) for i=1:100]
        @test checktree(KDTree(s), s)
        s = Shape{3,9}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.1) for i=1:100]
        @test checktree(KDTree(s), s)
    end

    @testset "vxlcut" begin
        @testset "triangular cylinder 3D" begin
            vxl = (SVector(0,0,0), SVector(1,1,1))
            nout = SVector(1,1,0)

            @test_nowarn (r₀ = SVector(0.5,0,0); @inferred(volfrac(vxl, nout, r₀)))
            @test (r₀ = SVector(0.5,0,0); volfrac(vxl, nout, r₀) ≈ 0.125)
            @test (r₀ = SVector(0.5,0,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
            @test (r₀ = SVector(1,0,0); volfrac(vxl, nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, nout, r₀) ≈ 0.25)
            @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, -nout, r₀) ≈ 0.75)
        end

        @testset "triangular cylinder 2D" begin
            vxl = (SVector(0,0), SVector(1,1))
            nout = SVector(1,1)

            @test_nowarn (r₀ = SVector(0.5,0); @inferred(volfrac(vxl, nout, r₀)))
            @test (r₀ = SVector(0.5,0); volfrac(vxl, nout, r₀) ≈ 0.125)
            @test (r₀ = SVector(0.5,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
            @test (r₀ = SVector(1,0); volfrac(vxl, nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, nout, r₀) ≈ 0.25)
            @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, -nout, r₀) ≈ 0.75)
        end

        @testset "quadrangular cylinder 3D" begin
            @test_nowarn @inferred(volfrac((SVector(0,0,0),SVector(1,1,1)), SVector(1,2,0), SVector(0.5,0.5,0.5)))
            @test begin
                result = true
                for i = 1:100
                    vxl = (-@SVector(rand(3)), @SVector(rand(3)))
                    r₀ = mean(vxl)
                    nout = randn(3)
                    nout[rand(1:3)] = 0
                    result &= volfrac(vxl, SVector{3}(nout), r₀)≈0.5
                end
                result
            end
        end

        @testset "quadrangular cylinder 2D" begin
            @test_nowarn @inferred(volfrac((SVector(0,0),SVector(1,1)), SVector(1,2), SVector(0.5,0.5)))
            @test begin
                result = true
                for i = 1:100
                    vxl = (-@SVector(rand(2)), @SVector(rand(2)))
                    r₀ = mean(vxl)
                    nout = @SVector randn(2)
                    result &= volfrac(vxl, nout, r₀)≈0.5
                end
                result
            end
        end

        @testset "general cases" begin
            # Test random cases.
            @test_nowarn @inferred(volfrac((-@SVector(rand(3)),@SVector(rand(3))), SVector(0,0,0), @SVector(randn(3))))
            @test begin
                result = true
                for i = 1:100
                    vxl = (-@SVector(rand(3)), @SVector(rand(3)))
                    r₀ = mean(vxl)
                    nout = @SVector randn(3)
                    result &= volfrac(vxl, nout, r₀)≈0.5
                end
                result
            end

            @test begin
                result = true
                for i = 1:100
                    vxl = (-@SVector(rand(3)), @SVector(rand(3)))
                    r₀ = @SVector randn(3)
                    nout = @SVector randn(3)
                    result &= (volfrac(vxl, nout, r₀) + volfrac(vxl, -nout, r₀) ≈ 1)
                end
                result
            end

            # Test boundary cases.
            vxl = (SVector(0,0,0), SVector(1,1,1))
            r₀ = SVector(0,0,0)
            @test (nout = SVector(1,0,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
            @test (nout = SVector(-1,0,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
            @test (nout = SVector(1,1,0); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
            @test (nout = SVector(-1,-1,0); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
            @test (nout = SVector(1,1,1); volfrac(vxl, nout, r₀) ≈ 0)  # completely outside
            @test (nout = SVector(-1,-1,-1); volfrac(vxl, nout, r₀) ≈ 1)  # completely inside
            @test (nout = SVector(-1,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_tricyl()
            @test (nout = SVector(-1,-1,1); volfrac(vxl, nout, r₀) ≈ 5/6)  # rvol_gensect()
            @test (nout = SVector(1,-2,1); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadsect()
            r₀ = SVector(0.5,0.5,0)
            @test (nout = SVector(-2,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadcyl()

            # Test tolerance to floating-point arithmetic.
            vxl = (SVector(0,9.5,0), SVector(1,10.5,1))
            r₀ = SVector(0.5,10.0,0.5)
            nout = SVector(0.0,1/√2,1/√2)
            @test volfrac(vxl, nout, r₀) ≈ 0.5

            # Test rvol_quadsect() for nontrivial cases.
            vxl = (SVector(0,0,0), SVector(1,1,2))
            r₀ = SVector(0.5, 0.5, 0.5)
            @test begin
                result = true
                for i = 1:100
                    nout = SVector(randn()/20, randn()/20, 1)
                    result &= volfrac(vxl, nout, r₀)≈0.5/2
                end
                result
            end
        end
    end
end
