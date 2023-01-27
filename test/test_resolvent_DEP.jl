@testset "Resolvent type construction       " begin
    # create array to be used as base
    # Ny = rand(20:100)
    # Nz = rand(2:10)
    # Nt = rand(2:10)
    Ny = 10
    Nz = 3
    Nt = 3
    ω = rand()*10
    β = rand()*10
    Re = rand()*200
    Ro = rand()
    U = zeros(Ny, Nz, Nt)
    dūdy = ones(Ny)
    ws = chebws(Ny)
    Dy = chebdiff(Ny)
    Dy2 = chebddiff(Ny)
    trunc = 5

    # construct resolvents at all the mode numbers with and without truncation
    R = Resolvent(U, ws, dūdy, ω, β, Re, Ro, Dy, Dy2)
    # R_trunc = Resolvent(U, ws, dūdy, ω, β, Re, Ro, Dy, Dy2; trunc=trunc)

    # generate resolvent at a specific mode number and compare to the container type
    kz = rand(1:Nz); kt = rand(1:Nt)
    H_test = ResolventAnalysis.resolvent_at_k(kz, kt, dūdy, ω, β, Re, Ro, Dy, Dy2)
    @test R.res_svd[kz, kt] == H_test
    # @test R_trunc.res_svd[kz, kt] == truncate_svd(svd(H_test, cholesky(ws)...))
end

@testset "Resolvent getter method           " begin
    
end

@testset "Resolvent property getter method  " begin

end

@testset "Resolvent type truncation method  " begin

end
