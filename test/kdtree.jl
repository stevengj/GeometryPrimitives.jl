@testset "KDTree" begin
    s = [Ball([i,0], 1) for i in 2:4]
    s0 = Ball([0,0], 1)
    s = [s0, s0, s0, s0, s...]
    @test_nowarn KDTree(s)  # must not generate StackOverflowError

    s = Shape{2,4}[Ball([i,0], 1, i) for i in 0:20]
    kd = KDTree(s)
    @test GeometryPrimitives.depth(kd) == 3
    @test findfirst([10.1,0], kd).data == 10
    @test findfirst([10.1,1], kd) == nothing
    @test checktree(kd, s)

    s = Shape{3,9}[Ball(SVector(randn(rng),randn(rng),randn(rng)), 0.01) for i=1:100]
    @test checktree(KDTree(s), s)
    s = Shape{3,9}[Ball(SVector(randn(rng),randn(rng),randn(rng)), 0.1) for i=1:100]
    @test checktree(KDTree(s), s)
end
