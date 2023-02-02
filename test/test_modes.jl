@testset "Converting SVD to mode types      " begin
    # convert to vectors of modes
    my_svd = svd(rand(rand(1:10), rand(1:10)))
    U_modes, S, V_modes = ResolventAnalysis.svd2modes(my_svd)

    @test U_modes isa Vector{ResolventAnalysis.ChannelMode{Float64}}
    @test V_modes isa Vector{ResolventAnalysis.ChannelMode{Float64}}

    for i in axes(S, 1)
        @test U_modes[i] == my_svd.U[:, i]
        @test V_modes[i] == my_svd.V[:, i]
    end
    @test S == svdvals(my_svd)
end