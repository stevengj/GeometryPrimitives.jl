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

            @test surfpt_nearby([3,4],s) == ([8,4],[1,0])  # handle point at center properly
            @test all([surfpt_nearby([3+ρ*sx*5,4],s)[1] ≈ [3+sx*5,4] for ρ = (one⁻⁻,1,one⁺⁺), sx = (-1,1)])
            @test all([surfpt_nearby([3,4+ρ*sy*5],s)[1] ≈ [3,4+sy*5] for ρ = (one⁻⁻,1,one⁺⁺), sy = (-1,1)])
            @test all([surfpt_nearby([3+ρ*sx*5/√2,4+ρ*sy*5/√2],s)[1] ≈ [3+sx*5/√2,4+sy*5/√2]
                       for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1), sy = (-1,1)])
            @test all([(p = [3+sx*5,4]; surfpt_nearby(p,s) ≈ (p,[sx,0])) for sx = (-1,1)])  # handle point on boundary properly
            @test all([(p = [3,4+sy*5]; surfpt_nearby(p,s) ≈ (p,[0,sy])) for sy = (-1,1)])  # handle point on boundary properly

            @test normal([-1,2],s) == normalize([-1,2] - [3,4])
            @test bounds(s) == ([-2,-1],[8,9])
            @test checkbounds(s)
            @test checkbounds(Sphere([1,2,3], 2))

        end

        @testset "Box" begin
            b = Box([0,0], [2,4])  # specify center and radii
            @test b == deepcopy(b)
            @test hash(b) == hash(deepcopy(b))
            @test [0.3,-1.5] ∈ b
            @test [0.3,-2.5] ∉ b

            @test surfpt_nearby([0,0],b) == ([1,0],[1,0])  # handle point at center properly
            @test all([surfpt_nearby([ρ*sx*1,0],b)[1] ≈ [sx*1,0] for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])
            @test all([surfpt_nearby([0,ρ*sy*2],b)[1] ≈ [0,sy*2] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
            @test all([(p = [sx*1,0]; surfpt_nearby(p,b) ≈ (p,[sx,0])) for sx = (-1,1)])  # handle point on boundary properly
            @test all([(p = [0,sy*2]; surfpt_nearby(p,b) ≈ (p,[0,sy])) for sy = (-1,1)])  # handle point on boundary properly
            @test surfpt_nearby([1.1,2.01], b)[1] ≈ [1.1,2]
            @test surfpt_nearby([1.01,2.1], b)[1] ≈ [1,2.1]

            @test normal([1.1,0],b) == [1,0]
            @test normal([-1.1,0],b) == [-1,0]
            @test normal([1.1,2.01],b) == [0,1]
            @test bounds(b) == ([-1,-2],[1,2])
            @test bounds(Box([0,0], [2,4], [1 1; 1 -1])) ≈ ([-3*√0.5,-3*√0.5], [3*√0.5,3*√0.5])
            @test checkbounds(b)
            @test checkbounds(Box([0,0], [2,4], [1 1; 1 -1]))
        end

        @testset "Box, rotated" begin
            ax1, ax2 = [1,-1], [1,1]
            r1, r2 = 1, 2  # "radii"
            br = Box([0,0], [2r1, 2r2], [ax1 -ax2])  # use -ax2 to check signs of axes don't matter

            R = [normalize(ax1) normalize(ax2)]  # rotation matrix

            Cin = R * (GeometryPrimitives.signmatrix(br) .* (one⁻ .* [r1,r2]))  # around corners, inside
            Cout = R * (GeometryPrimitives.signmatrix(br) .* (one⁺ .* [r1,r2]))  # around corners, outside
            for j = 1:4; @test Cin[:,j] ∈ br; end
            for j = 1:4; @test Cout[:,j] ∉ br; end

            @test br == deepcopy(br)
            @test hash(br) == hash(deepcopy(br))

            n1, n2 = normalize.((ax1, ax2))
            @test surfpt_nearby([0,0],br) ≈ (n1,n1)  # handle point at center properly
            @test all([surfpt_nearby(ρ*s1*r1*n1,br)[1] ≈ s1*r1*n1 for ρ = (one⁻⁻,one⁺⁺), s1 = (-1,1)])
            @test all([surfpt_nearby(ρ*s2*r2*n2,br)[1] ≈ s2*r2*n2 for ρ = (one⁻⁻,one⁺⁺), s2 = (-1,1)])
            @test all([(p = s1*r1*n1; surfpt_nearby(p,br) ≈ (p,s1*n1)) for s1 = (-1,1)])  # handle point on boundary properly
            @test all([(p = s2*r2*n2; surfpt_nearby(p,br) ≈ (p,s2*n2)) for s2 = (-1,1)])  # handle point on boundary properly

            @test normal(R*[1.1r1, 0], br) ≈ R*[1,0]
            @test normal(R*[-1.1r1, 0], br) ≈ R*[-1,0]
            @test normal(R*[0, 1.1r2], br) ≈ R*[0,1]
            @test normal(R*[0, -1.1r2], br) ≈ R*[0,-1]
            @test normal(R*[1.1r1, 1.01r2], br) ≈ R*[0,1]

            xmax = (R*[r1,r2])[1]
            ymax = (R*[-r1,r2])[2]
            @test bounds(br) ≈ (-[xmax,ymax], [xmax,ymax])
            @test checkbounds(br)
        end

        @testset "Box, skewed" begin
            ax1, ax2 = normalize.(([1,-1], [0,1]))
            r1, r2 = 1, 1  # "radii"
            bs = Box([0,0], [2r1, 2r2], [-ax1 ax2])  # use -ax1 to check signs of axes don't matter

            @test bs == deepcopy(bs)
            @test hash(bs) == hash(deepcopy(bs))

            n1, n2 = normalize.(([1,0], [1,1]))
            cosθ = (-ax1) ⋅ (-n1)
            @test surfpt_nearby([0,0],bs) ≈ (-cosθ*n1,-n1)  # handle point at center properly
            @test all([surfpt_nearby(s2*(r2*ax2+∆ρ*n2),bs)[1] ≈ s2*r2*ax2 for ∆ρ = (-0.1,0.1), s2 = (-1,1)])
            @test all([surfpt_nearby(s1*(r1*ax1+∆ρ*n1),bs)[1] ≈ s1*r1*ax1 for ∆ρ = (-0.1,0.1), s1 = (-1,1)])
            @test all([(p = s1*r1*ax1; surfpt_nearby(p,bs) ≈ (p,s1*n1)) for s1 = (-1,1)])  # handle point on boundary properly
            @test all([(p = s2*r2*ax2; surfpt_nearby(p,bs) ≈ (p,s2*n2)) for s2 = (-1,1)])  # handle point on boundary properly

            @test norm(normal([0,1], bs)) ≈ 1

            xmax = (r1*ax1+r2*ax2)[1]
            ymax = (r2*ax2-r1*ax1)[2]
            @test bounds(bs) ≈ (-[xmax,ymax],[xmax,ymax])
            @test checkbounds(bs)
        end

        @testset "Ellipsoid" begin
            e = Ellipsoid([0,0], [1,2])
            @test e == deepcopy(e)
            @test hash(e) == hash(deepcopy(e))
            @test [0.3,2*sqrt(1 - 0.3^2)-0.01] ∈ e
            @test [0.3,2*sqrt(1 - 0.3^2)+0.01] ∉ e

            @test surfpt_nearby([0,0],e) ≈ ([1,0], [1,0])  # handle point at center properly
            @test all([surfpt_nearby([ρ*sx*1,0],e)[1] ≈ [sx*1,0] for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])
            @test all([surfpt_nearby([0,ρ*sy*2],e)[1] ≈ [0,sy*2] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
            @test all([(p = [sx*1,0]; surfpt_nearby(p,e) ≈ (p,[sx,0])) for sx = (-1,1)])  # handle point on boundary properly
            @test all([(p = [0,sy*2]; surfpt_nearby(p,e) ≈ (p,[0,sy])) for sy = (-1,1)])  # handle point on boundary properly

            @test normal([1.1,0],e) ≈ [1,0]
            @test normal([-1.1,0],e) ≈ [-1,0]
            @test normal([0,2.01],e) ≈ [0,1]
            @test bounds(e) == ([-1,-2],[1,2])
            @test checkbounds(e)
            @test checkbounds(Ellipsoid([0,0], [1,2], [1 1; 1 -1]))

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

            @test surfpt_nearby([0,0],er) ≈ (2*R[:,1], R[:,1])  # handle point at center properly
            @test all([surfpt_nearby(R*[ρ*sx*2,0],er)[1] ≈ R*[sx*2,0] for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])
            @test all([surfpt_nearby(R*[0,ρ*sy*3],er)[1] ≈ R*[0,sy*3] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
            @test all([(p = R*[sx*2,0]; surfpt_nearby(p,er) ≈ (p,R*[sx,0])) for sx = (-1,1)])  # handle point on boundary properly
            @test all([(p = R*[0,sy*3]; surfpt_nearby(p,er) ≈ (p,R*[0,sy])) for sy = (-1,1)])  # handle point on boundary properly

            # Test the normal vector at the two bounding points are the x- and y-directions.
            @test normal(bp1, er) ≈ [1,0]
            @test normal(bp2, er) ≈ [0,1]

            xmax, ymax = bp1[1], bp2[2]
            @test bounds(er) == ([-xmax, -ymax], [xmax, ymax])
            @test checkbounds(er)

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

            @test surfpt_nearby([0,0,0],c) == ([0,0,1.1],[0,0,1])  # handle point at center properly
            @test all([surfpt_nearby([ρ*sx*0.3,0,0],c)[1] ≈ [sx*0.3,0,0] for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])
            @test all([surfpt_nearby([0,ρ*sy*0.3,0],c)[1] ≈ [0,sy*0.3,0] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
            @test all([surfpt_nearby([0,0,ρ*sz*1.1],c)[1] ≈ [0,0,sz*1.1] for ρ = (one⁻⁻,one⁺⁺), sz = (-1,1)])
            @test all([(p = [0,sy*0.3,0]; surfpt_nearby(p,c) ≈ (p,[0,sy,0])) for sy = (-1,1)])  # handle point on boundary properly
            @test all([(p = [0,0,sz*1.1]; surfpt_nearby(p,c) ≈ (p,[0,0,sz])) for sz = (-1,1)])  # handle point on boundary properly

            @test normal([0.1,0.2,-1.3], c) == [0.1,0.2,0] / hypot(0.1,0.2)
            @test normal([0.1,0.2,-1.11], c) == [0,0,-1]
            @test normal([0.31, 0, 0.3], c) == [1,0,0]
            @test bounds(c) ≈ ([-0.3,-0.3,-1.1],[0.3,0.3,1.1])
            @test checkbounds(c)
            @test checkbounds(Cylinder([1,17,44], 0.3, [1,-2,3], 1.1))
        end

        @testset "Cylinder, rotated" begin
            ax1 = normalize([1,0,1])
            ax2 = normalize([1,0,-1])
            cr = Cylinder([0,0,0], 0.3, -ax1, 2.2)  # use -ax1 to make sure sign of axis doesn't matter
            @test cr == deepcopy(cr)
            @test hash(cr) == hash(deepcopy(cr))

            @test surfpt_nearby([0,0,0],cr) ≈ (-1.1ax1,-ax1)  # handle point at center properly
            @test all([(p = s1*1.1*ax1; surfpt_nearby(p,cr) ≈ (p,s1*ax1)) for s1 = (-1,1)])  # handle point on axis properly
            @test all([surfpt_nearby(ρ*s3*1.1ax1,cr)[1] ≈ s3*1.1ax1 for ρ = (one⁻⁻,one⁺⁺), s3 = (-1,1)])
            @test all([surfpt_nearby([0,ρ*sy*0.3,0],cr)[1] ≈ [0,sy*0.3,0] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
            @test all([surfpt_nearby(ρ*sx*0.3*ax2,cr)[1] ≈ sx*0.3*ax2 for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])
            @test all([(p = [0,sy*0.3,0]; surfpt_nearby(p,cr) ≈ (p,[0,sy,0])) for sy = (-1,1)])  # handle point on boundary properly
            @test all([(p = sx*0.3*ax2; surfpt_nearby(p,cr) ≈ (p, sx*ax2)) for sx = (-1,1)])  # handle point on boundary properly

            @test all([normal(ρ*s3*1.1ax1,cr) ≈ s3*ax1 for ρ = (one⁻⁻,one⁺⁺), s3 = (-1,1)])
            @test all([normal([0,ρ*sy*0.3,0],cr) ≈ sy*[0,1,0] for ρ = (one⁻⁻,one⁺⁺), sy = (-1,1)])
            @test all([normal(ρ*sx*0.3*ax2,cr) ≈ sx*ax2 for ρ = (one⁻⁻,one⁺⁺), sx = (-1,1)])

            @test bounds(cr) ≈ (-[(1.1+0.3)/√2,0.3,(1.1+0.3)/√2], [(1.1+0.3)/√2,0.3,(1.1+0.3)/√2])
            @test checkbounds(cr)
        end
    end

    @testset "KDTree" begin
        s = Shape{2,4}[Sphere([i,0], 1, i) for i in 0:20]
        kd = KDTree(s)
        @test GeometryPrimitives.depth(kd) == 4
        @test get(findin([10.1,0], kd)).data == 10
        @test isnull(findin([10.1,1], kd))
        @test checktree(kd, s)
        s = Shape{3,9}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.01) for i=1:100]
        @test checktree(KDTree(s), s)
        s = Shape{3,9}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.1) for i=1:100]
        @test checktree(KDTree(s), s)
    end

    @testset "stability" begin
        b = Box([0,0,0],[1,1,1])
        c = Cylinder([0,0,0],1,[0,0,1],1)
        e = Ellipsoid([0,0,0],[1,1,1])
        s = [Sphere([i,0,0], 1) for i in 0:20]
        svec = [b, c, e, s...]  # long enough to invoke findin(::SVector, ::KDTree) with kd.left and kd.right

        kd = KDTree(svec)

        # The purpose of these tests are not
        @test_nowarn @inferred(surfpt_nearby([0,0,0], get(findin([0,0,0], kd))))
        @test_nowarn @inferred(normal([0,0,0], get(findin([0,0,0], kd))))
        @test_nowarn @inferred([0,0,0] ∈ get(findin([0,0,0], kd)))
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
            @test (nout = SVector(-1,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_tricyl()
            @test (nout = SVector(-1,-1,1); volfrac(vxl, nout, r₀) ≈ 5/6)  # rvol_gensect()
            @test (nout = SVector(1,-2,1); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadsect()
            r₀ = SVector(0.5,0.5,0)
            @test (nout = SVector(-2,1,0); volfrac(vxl, nout, r₀) ≈ 0.5)  # rvol_quadcyl()

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
