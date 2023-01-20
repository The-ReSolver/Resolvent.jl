@testset "Resolvent at specific mode number " begin
    Re = rand()*200
    nz = rand(1:100)
    nt = rand(0:100)
    fund_freq = rand()*10
    β = rand()*10
    N = rand(4:200)

    # TODO: talk to Sean about the negated time derivative
    H0 = o.quick_example(Re, nz, nt, β, fund_freq, N, 3)
    H_me = ResolventAnalysis.resolvent_at_k(nz, nt, ones(N), fund_freq, β, Re, 0.0, chebdiff(N), chebddiff(N))

    @test H0 ≈ H_me
end

@testset "Resolvent type construction       " begin
    
end

@testset "Resolvent getter method           " begin
    
end

@testset "Resolvent property getter method  " begin

end
