using GeometryPrimitives, StaticArrays, Base.Test
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

"check the bounding box of o with randomized trials"
function checkbounds{N}(o::Object{N}, ntrials=10^4)
    lb,ub = bounds(o)
    for i = 1:ntrials
        x = randnb(lb,ub)
        x ∉ o || inbounds(x,lb,ub) || return false
    end
    return true
end

function checktree{N}(t::KDTree{N}, olist::Vector{Object{N}}, ntrials=10^3)
    lb = SVector{N}(fill(Inf,N))
    ub = SVector{N}(fill(-Inf,N))
    for i in eachindex(olist)
        lbi,ubi = bounds(olist[i])
        lb = min(lb,lbi)
        ub = max(ub,ubi)
    end
    for i = 1:ntrials
        x = randnb(lb,ub)
        ot = findin(x, t)
        ol = findin(x, olist)
        isnull(ot) == isnull(ol) || return false
        isnull(ot) || get(ot) == get(ol) || return false
    end
    return true
end

@testset "GeometryPrimitives" begin

    @testset "Object" begin
        @testset "Sphere" begin
            s = Sphere([3,4], 5)
            @test ndims(s) == 2
            @test [3,9] ∈ s
            @test [3,9.1] ∉ s
            @test normal([-1,2],s) == normalize([-1,2] - [3,4])
            @test bounds(s) == ([-2,-1],[8,9])
            @test checkbounds(s)
            @test checkbounds(Sphere([1,2,3], 2))
        end

        @testset "Box" begin
            b = Box([0,0], [2,4])
            @test [0.3,-1.5] ∈ b
            @test [0.3,-2.5] ∉ b
            @test normal([1.1,0],b) == [1,0]
            @test normal([-1.1,0],b) == [-1,0]
            @test normal([1.1,2.01],b) == [0,1]
            @test bounds(b) == ([-1,-2],[1,2])
            @test bounds(Box([0,0], [2,4], [1 1; 1 -1])) ≈ ([-3*√0.5,-3*√0.5], [3*√0.5,3*√0.5])
            @test checkbounds(b)
            @test checkbounds(Box([0,0], [2,4], [1 1; 1 -1]))
        end

        @testset "Ellipsoid" begin
            e = Ellipsoid([0,0], [2,4])
            @test [0.3,2*sqrt(1 - 0.3^2)-0.01] ∈ e
            @test [0.3,2*sqrt(1 - 0.3^2)+0.01] ∉ e
            @test normal([1.1,0],e) == [1,0]
            @test normal([-1.1,0],e) == [-1,0]
            @test normal([0,2.01],e) == [0,1]
            @test bounds(e) == ([-1,-2],[1,2])
            @test checkbounds(e)
            @test checkbounds(Ellipsoid([0,0], [2,4], [1 1; 1 -1]))
        end

        @testset "Cylinder" begin
            c = Cylinder([0,0,0], 0.3, [0,0,1], 2.2)
            @test [0.2,0.2,1] ∈ c
            @test SVector(0.2,0.2,1.2) ∉ c
            @test [0.2,0.25,1] ∉ c
            @test normal([0.1,0.2,-1.3], c) == [0,0,-1]
            @test normal([0.31, 0, 0.3], c) == [1,0,0]
            @test bounds(c) == ([-0.3,-0.3,-1.1],[0.3,0.3,1.1])
            @test checkbounds(c)
            @test checkbounds(Cylinder([1,17,44], 0.3, [1,-2,3], 1.1))
        end
    end

    @testset "KDTree" begin
        o = Object{2}[Sphere([i,0], 1, i) for i in 0:20]
        kd = KDTree(o)
        @test GeometryPrimitives.depth(kd) == 4
        @test get(findin([10.1,0], kd)).data == 10
        @test isnull(findin([10.1,1], kd))
        @test checktree(kd, o)
        o = Object{3}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.01) for i=1:100]
        @test checktree(KDTree(o), o)
        o = Object{3}[Sphere(SVector(randn(rng),randn(rng),randn(rng)), 0.1) for i=1:100]
        @test checktree(KDTree(o), o)
    end

end
