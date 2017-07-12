using GeometryPrimitives, StaticArrays, Base.Test

const rtol = Base.rtoldefault(Float64)
const one⁻ = 1 - rtol  # slightly less than 1
const one⁺ = 1 + rtol  # slightly greater than 1

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

function checktree(t::KDTree{N}, slist::Vector{Shape{N}}, ntrials=10^3) where {N}
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
            @test @inferred(ndims(s)) == 2
            @test @inferred([3,9] ∈ s)
            @test [3,9.1] ∉ s
            @test @inferred(normal([-1,2],s)) == normalize([-1,2] - [3,4])
            @test @inferred(bounds(s)) == ([-2,-1],[8,9])
            @test checkbounds(s)
            @test checkbounds(Sphere([1,2,3], 2))
        end

        @testset "Box" begin
            b = Box([0,0], [2,4])  # specify center and radii
            @test b == deepcopy(b)
            @test hash(b) == hash(deepcopy(b))
            @test @inferred([0.3,-1.5] ∈ b)
            @test [0.3,-2.5] ∉ b
            @test @inferred(normal([1.1,0],b)) == [1,0]
            @test normal([-1.1,0],b) == [-1,0]
            @test normal([1.1,2.01],b) == [0,1]
            @test @inferred(bounds(b)) == ([-1,-2],[1,2])
            @test bounds(Box([0,0], [2,4], [1 1; 1 -1])) ≈ ([-3*√0.5,-3*√0.5], [3*√0.5,3*√0.5])
            @test checkbounds(b)
            @test checkbounds(Box([0,0], [2,4], [1 1; 1 -1]))
        end

        @testset "Box, rotated" begin
            ax1, ax2 = [1,-1], [1,1]
            r1, r2 = 1, 2  # "radii"
            br = Box([0,0], [2r1, 2r2], [ax1 ax2])

            R = [normalize(ax1) normalize(ax2)]  # rotation matrix

            Cin = R * (GeometryPrimitives.signmatrix(br) .* (one⁻ .* [r1,r2]))  # around corners, inside
            Cout = R * (GeometryPrimitives.signmatrix(br) .* (one⁺ .* [r1,r2]))  # around corners, outside
            for j = 1:4; @test @inferred(Cin[:,j] ∈ br); end
            for j = 1:4; @test Cout[:,j] ∉ br; end

            @test br == deepcopy(br)
            @test hash(br) == hash(deepcopy(br))
            @test @inferred(normal(R*[1.1r1, 0], br)) ≈ R*[1,0]
            @test normal(R*[-1.1r1, 0], br) ≈ R*[-1,0]
            @test normal(R*[0, 1.1r2], br) ≈ R*[0,1]
            @test normal(R*[0, -1.1r2], br) ≈ R*[0,-1]
            @test normal(R*[1.1r1, 1.01r2], br) ≈ R*[0,1]

            xmax = (R*[r1,r2])[1]
            ymax = (R*[-r1,r2])[2]
            @test @inferred(bounds(br)) ≈ (-[xmax,ymax], [xmax,ymax])
            @test checkbounds(br)
        end

        @testset "Ellipsoid" begin
            e = Ellipsoid([0,0], [1,2])
            @test e == deepcopy(e)
            @test hash(e) == hash(deepcopy(e))
            @test @inferred([0.3,2*sqrt(1 - 0.3^2)-0.01] ∈ e)
            @test [0.3,2*sqrt(1 - 0.3^2)+0.01] ∉ e
            @test @inferred(normal([1.1,0],e)) == [1,0]
            @test normal([-1.1,0],e) == [-1,0]
            @test normal([0,2.01],e) == [0,1]
            @test @inferred(bounds(e)) == ([-1,-2],[1,2])
            @test checkbounds(e)
            @test checkbounds(Ellipsoid([0,0], [1,2], [1 1; 1 -1]))

            b = Box([0,0], [2,4])
            eb = Ellipsoid(b)
            @test e == eb
            @test hash(e) == hash(eb)
        end

        @testset "Ellipsoid, rotated" begin
            θ = π/3
            er = Ellipsoid([0,0], [1,2], [cos(θ) sin(θ); sin(θ) -cos(θ)])
            bp = GeometryPrimitives.boundpts(er)

            bp1, bp2 = bp[:,1], bp[:,2]

            # Test the two bounding points are on the ellipsoid perimeter.
            @test er == deepcopy(er)
            @test hash(er) == hash(deepcopy(er))
            @test (@inferred(one⁻ * bp1 ∈ er)) && (one⁻ * bp2 ∈ er)
            @test (one⁺ * bp1 ∉ er) && (one⁺ * bp2 ∉ er)

            # Test the normal vector at the two bounding points are the x- and y-directions.
            @test @inferred(normal(bp1, er)) ≈ [1,0]
            @test normal(bp2, er) ≈ [0,1]

            xmax, ymax = bp1[1], bp2[2]
            @test @inferred(bounds(er)) == ([-xmax, -ymax], [xmax, ymax])
            @test checkbounds(er)

            br = Box([0,0], [2,4], [cos(θ) sin(θ); sin(θ) -cos(θ)])
            ebr = Ellipsoid(br)
            @test er == ebr
            @test hash(er) == hash(ebr)
        end

        @testset "Cylinder" begin
            c = Cylinder([0,0,0], 0.3, [0,0,1], 2.2)
            @test c == deepcopy(c)
            @test hash(c) == hash(deepcopy(c))
            @test @inferred([0.2,0.2,1] ∈ c)
            @test SVector(0.2,0.2,1.2) ∉ c
            @test [0.2,0.25,1] ∉ c
            @test @inferred(normal([0.1,0.2,-1.3], c)) == [0,0,-1]
            @test normal([0.31, 0, 0.3], c) == [1,0,0]
            @test @inferred(bounds(c)) ≈ ([-0.3,-0.3,-1.1],[0.3,0.3,1.1])
            @test checkbounds(c)
            @test checkbounds(Cylinder([1,17,44], 0.3, [1,-2,3], 1.1))
        end
    end

    @testset "KDTree" begin
        s = Shape{2}[Sphere([i,0], 1, i) for i in 0:20]
        kd = KDTree(s)
        @test GeometryPrimitives.depth(kd) == 4
        @test get(findin([10.1,0], kd)).data == 10
        @test isnull(findin([10.1,1], kd))
        @test checktree(kd, s)
        s = Shape{3}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.01) for i=1:100]
        @test checktree(KDTree(s), s)
        s = Shape{3}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.1) for i=1:100]
        @test checktree(KDTree(s), s)
    end

    @testset "vxlcut" begin
        @testset "triangular cylinder 3D" begin
            vxl = (SVector(0,0,0), SVector(1,1,1))
            nout = SVector(1,1,0)

            @test (r₀ = SVector(0.5,0,0); @inferred(volfrac(vxl, nout, r₀)) ≈ 0.125)
            @test (r₀ = SVector(0.5,0,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
            @test (r₀ = SVector(1,0,0); volfrac(vxl, nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, nout, r₀) ≈ 0.25)
            @test (r₀ = SVector(1,0,0); nout = SVector(1,2,0); volfrac(vxl, -nout, r₀) ≈ 0.75)
        end

        @testset "triangular cylinder 2D" begin
            vxl = (SVector(0,0), SVector(1,1))
            nout = SVector(1,1)

            @test (r₀ = SVector(0.5,0); @inferred(volfrac(vxl, nout, r₀)) ≈ 0.125)
            @test (r₀ = SVector(0.5,0); volfrac(vxl, -nout, r₀) ≈ 0.875)
            @test (r₀ = SVector(1,0); volfrac(vxl, nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0); volfrac(vxl, -nout, r₀) ≈ 0.5)
            @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, nout, r₀) ≈ 0.25)
            @test (r₀ = SVector(1,0); nout = SVector(1,2); volfrac(vxl, -nout, r₀) ≈ 0.75)
        end

        @testset "quadrangular cylinder 3D" begin
            @test begin
                result = true
                for i = 1:100
                    vxl = (-@SVector(rand(3)), @SVector(rand(3)))
                    r₀ = mean(vxl)
                    nout = randn(3)
                    nout[rand(1:3)] = 0
                    result &= @inferred(volfrac(vxl, SVector{3}(nout), r₀)) ≈ 0.5
                end
                result
            end
        end

        @testset "quadrangular cylinder 2D" begin
            @test begin
                result = true
                for i = 1:100
                    vxl = (-@SVector(rand(2)), @SVector(rand(2)))
                    r₀ = mean(vxl)
                    nout = @SVector randn(2)
                    result &= @inferred(volfrac(vxl, nout, r₀)) ≈ 0.5
                end
                result
            end
        end

        @testset "general cases" begin
            # Test random cases.
            @test begin
                result = true
                for i = 1:100
                    vxl = (-@SVector(rand(3)), @SVector(rand(3)))
                    r₀ = mean(vxl)
                    nout = @SVector randn(3)
                    result &= @inferred(volfrac(vxl, nout, r₀)) ≈ 0.5
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
                    result &= @inferred(volfrac(vxl, nout, r₀)) ≈ 0.5/2
                end
                result
            end
        end
    end
end
