@testset "Truncation type construction      " begin
    # ensure the correct truncation type is obtained from constructor
    @test TruncateSVD(0) isa FullSVD
    randint = rand(1:10)
    @test TruncateSVD(randint) isa TruncateSVD{randint}
end

@testset "Truncation performed correctly    " begin
    # contruct SVD of random matrix and perform truncation
    trunc_length = rand(1:4)
    A = rand(rand(5:10), rand(5:10))
    U, S, V = svd(A)
    A_trunc = ResolventAnalysis.truncate_svd!(svd(A), TruncateSVD(trunc_length))

    # check SVDs are the same as their individually truncated versions
    @test A_trunc.U == U[:, 1:trunc_length]
    @test A_trunc.S == S[1:trunc_length]
    @test A_trunc.Vt == V[:, 1:trunc_length]'

    # check when truncation is not positive
    A_notrunc = ResolventAnalysis.truncate_svd!(svd(A), TruncateSVD(rand(-10:0)))
    @test A_notrunc.U == U
    @test A_notrunc.S == S
    @test A_notrunc.V == V
end
