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
    A_trunc = ResolventAnalysis.truncate_svd(svd(A), TruncateSVD(trunc_length))

    # check SVDs are the same as their individually truncated versions
    @test A_trunc.U == U[:, 1:trunc_length]
    @test A_trunc.S == S[1:trunc_length]
    @test A_trunc.Vt == V[:, 1:trunc_length]'

    # check when truncation is not positive
    A_notrunc = ResolventAnalysis.truncate_svd(svd(A), TruncateSVD(rand(-10:0)))
    @test A_notrunc.U == U
    @test A_notrunc.S == S
    @test A_notrunc.V == V
end

@testset "Cholesky decomposition of weights " begin
    # generate weights
    N = rand(3:100)
    ws = chebws(N)

    # compute cholesky decomposition matrices
    L, L_inv = cholesky(ws)

    Z = zeros(N, N)
    @test L*L' ≈   [Diagonal(ws) Z            Z           ;
                    Z            Diagonal(ws) Z           ;
                    Z            Z            Diagonal(ws);]
    
    @test L_inv*L ≈    [I Z Z Z;
                        Z I Z Z;
                        Z Z I Z;]
end

@testset "Resolvent at specific mode number " begin
    # generate random flow parameters
    Re = rand()*200
    kz = rand()*100
    kt = rand()*100
    N = rand(4:200)

    # compute Resolvents
    H0 = o.quick_example(Re, kz, kt, N, 3)
    H_me = ResolventAnalysis.resolvent_at_k(kz, kt, ones(N), Re, 0.0, chebdiff(N), chebddiff(N))

    @test H_me ≈ H0[1:4*N, 1:3*N]
end

@testset "Matrix SVD at specific mode number" begin
    # generate random flow parameters
    Re = rand()*200
    kz = rand()*100
    kt = rand()*100
    N = rand(4:200)
    ws = chebws(N)

    # compute Resolvents
    _, U_sean, S_sean, V_sean = o.quick_example(Re, kz, kt, N, nout=4); S_sean = reshape(S_sean, 4*N)
    H_me = ResolventAnalysis.resolvent_at_k(kz, kt, ones(N), Re, 0.0, chebdiff(N), chebddiff(N))
    SVD_me = svd(H_me, ws)

    @test abs.(SVD_me.U) ≈ abs.(U_sean[1:3*N, 1:length(SVD_me.S)])
    @test SVD_me.S ≈ S_sean[1:length(SVD_me.S)]
    @test abs.(SVD_me.V) ≈ abs.(V_sean[1:3*N, 1:length(SVD_me.S)])
end
