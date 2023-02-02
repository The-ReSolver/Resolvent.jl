@testset "Converting SVD to mode types      " begin
    # convert to vectors of modes
    my_svd = svd(rand(rand(1:10), rand(1:10)))
    U_modes, S, V_modes = ResolventAnalysis.svd2modes(my_svd, ones(Float64, size(my_svd, 1)))

    @test U_modes isa Vector{ResolventAnalysis.ChannelMode{Float64}}
    @test V_modes isa Vector{ResolventAnalysis.ChannelMode{Float64}}

    U_same = true; V_same = true
    for i in axes(S, 1)
        U_modes[i] == my_svd.U[:, i] ? nothing : U_same = false
        V_modes[i] == my_svd.V[:, i] ? nothing : V_same = false
    end
    @test U_same
    @test V_same
    @test S == svdvals(my_svd)
end
