using Test
using Random
using LinearAlgebra

using ResolventAnalysis
using ChebUtils

# # set up the PyCall interface to be able to use Sean's Matlab code directly
using PyCall
o = pyimport("oct2py").octave
o.addpath("./primitive_resolvent/")

@testset "Truncation performed correctly    " begin
    # contruct SVD of random matrix and perform truncation
    trunc_length = rand(1:4)
    A = rand(rand(5:10), rand(5:10))
    U, S, V = svd(A)
    svd_trunc = ResolventAnalysis.truncateSVD(svd(A), trunc_length)

    # check SVDs are the same as their individually truncated versions
    @test svd_trunc.U == U[:, 1:trunc_length]
    @test svd_trunc.S == S[1:trunc_length]
    @test svd_trunc.Vt == V[:, 1:trunc_length]'
end

@testset "Cholesky decomposition of weights " begin
    # generate weights
    N = rand(3:100)
    ws = chebws(N)

    # compute cholesky decomposition matrices
    L, L_inv = ResolventAnalysis.cholesky(ws)

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
    H0 = o.quick_example(Re, kz, kt, N, nout=1)
    H_me = Resolvent(N, chebdiff(N), chebddiff(N))

    @test H_me(kz, kt, ones(N), Re, 0.0) ≈ H0[1:4*N, 1:3*N]
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
    H_me = Resolvent(N, chebdiff(N), chebddiff(N))
    SVD_me = svd(H_me(kz, kt, ones(N), Re, 0.0), ws, 5)

    @test abs.(SVD_me.U) ≈ abs.(U_sean[1:3*N, 1:5])
    @test SVD_me.S ≈ S_sean[1:5]
    @test abs.(SVD_me.V) ≈ abs.(V_sean[1:3*N, 1:5])
end
