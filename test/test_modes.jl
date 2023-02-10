@testset "Abstract interface test function  " begin
    # define a bunch of dummy mode types
    struct ModeA{N, T} <: ResolventAnalysis.AbstractMode{N, T}; end
    struct ModeB{N, T} <: ResolventAnalysis.AbstractMode{N, T}
        mode::Array{T, N}
        ws::Array{T, N}
    end
    struct ModeC{N, T} <: ResolventAnalysis.AbstractMode{N, T}
        mode::Array{T, N}
        ws::Array{T, N}
    end
    LinearAlgebra.dot(::ModeC{N}, ::AbstractArray{<:Any, N}) where {N} = "not a number"
    struct ModeD{T} <: ResolventAnalysis.AbstractMode{1, T}
        mode::Vector{T}
        ws::Vector{T}
    end
    LinearAlgebra.dot(::ModeD, ::AbstractVector) = 1.0

    @test_throws Exception test_modeinterface(ModeA())
    @test_throws MethodError test_modeinterface(ModeB(rand(2, 2), rand(2, 2)))
    @test_throws AssertionError test_modeinterface(ModeC(rand(3, 4, 2), rand(3, 4, 2)))
    @test test_modeinterface(ModeD(rand(5), rand(5))) === nothing
end

@testset "Converting SVD to mode types      " begin
    # convert to vectors of modes
    my_svd = svd(rand(rand(1:10), rand(1:10)))
    U_modes, S, V_modes = ResolventAnalysis.svd2channelmodes(my_svd, ones(Float64, size(my_svd, 1)))

    @test U_modes isa Vector{ChannelMode{Float64}}
    @test V_modes isa Vector{ChannelMode{Float64}}

    U_same = true; V_same = true
    for i in axes(S, 1)
        U_modes[i] == my_svd.U[:, i] ? nothing : U_same = false
        V_modes[i] == my_svd.V[:, i] ? nothing : V_same = false
    end
    @test U_same
    @test V_same
    @test S == svdvals(my_svd)
end

@testset "Channel mode inner-product        " begin

end

@testset "Channel mode projections          " begin

end
