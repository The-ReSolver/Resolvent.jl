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

@testset "Cholesky decomposition of weights " begin
    # generate weights
    N = rand(3:100)
    ws = chebws(N)

    Z = zeros(N, N)
    @test cholesky(ws)*cholesky(ws)' ≈ [Diagonal(ws) Z            Z            Z;
                                        Z            Diagonal(ws) Z            Z;
                                        Z            Z            Diagonal(ws) Z;
                                        Z            Z            Z            Z]
end

@testset "Matrix SVD at specific mode number" begin
    # generate random flow parameters
    # Re = rand()*200
    # nz = rand(1:100)
    # nt = rand(0:100)
    # fund_freq = rand()*10
    # β = rand()*10
    # N = rand(4:200)
    # ws = chebws(N)
    Re = 10.0
    nz = 1
    nt = 1
    fund_freq = 2π
    β = 1.0
    N = 10
    ws = chebws(N)

    # compute Resolvents
    _, U_sean, S_sean, V_sean = o.quick_example(Re, nz, nt, β, fund_freq, N, 3, nout=4)
    H_me = ResolventAnalysis.resolvent_at_k(nz, nt, ones(N), fund_freq, β, Re, 0.0, chebdiff(N), chebddiff(N))
    SVD_me = svd(H_me, cholesky(ws))

    # @test U_sean ≈ SVD_me.U
    # @test S_sean ≈ SVD_me.S
    # @test V_sean ≈ SVD_me.V
end
