# This file contains the type definition for the Resolvent operator, stored in
# terms of its SVD decomposition.

# TODO: find known test case.
# TODO: test utility functions for the truncation.
# TODO: implement method to reconstruct Resolvent from the SVD object.
# TODO: implement test for Resolvent operator.
# TODO: implement test for it's SVD decomposition.

struct Resolvent{SIZE<:NTuple{3, Int}, TRUNC<:Bool, T<:Real, N<:Int}
    res_svd::Array{SVD, N}

    function Resolvent(Hs::AbstractArray{Matrix{Complex{T}}, N}, trunc::TruncateSVD) where {T, N}
        new{(size(H[s1], 1), size(Hs)...), true, T, N}(truncate_svd!(svd(Hs), trunc))
    end

    function Resolvent(Hs::AbstractArray{Matrix{Complex{T}}, N}, trunc::FullSVD) where {T, N}
        new{(size(H[s1], 1), size(Hs)...), false, T, N}(truncate_svd!(svd(Hs), trunc))
    end
end

function Resolvent(U::AbstractArray{T, 3}, dūdy::Vector{S}, ω::S, β::S, Re::S, Ro::S, Dy::Matrix{S}, Dy2::Matrix{S}, trunc::Int=0) where {T, S}
    # initialise array to hold resolvents
    Hs = [Matrix{Complex{T}}(undef, 4*size(U, 1), 3*size(U, 1)) for kt in 1:size(U, 3), kz in 1:size(U, 2)]

    # loop over mode numbers and generate resolvent at each of them
    for kt in 1:size(U, 3), kz in 1:size(U, 2)
        Hs[kz, kt] = resolvent_at_k(kt, kz, dūdy, ω, β, Re, Ro, Dy, Dy2)
    end

    Resolvent(Hs, TruncateSVD(trunc))
end

# TODO: test these indexing methods
Base.getindex(resolvent::Resolvent, i::Int) = resolvent.res_svd[i] # ! is this first method required?
Base.getindex(resolvent::Resolvent, I...) = resolvent.res_svd[I...]

function Base.getproperty(resolvent::Resolvent{SIZE}, sym::Symbol) where {SIZE}
    if sym == :H
        # initialise matrix of Resolvents
        Hs = Matrix{Matrix{Complex{T}}}(undef, SIZE[2:3]...)

        # loop over mode numbers computing the Resolvent from its SVD
        for kt in 1:SIZE[3], kz in 1:SIZE[2]
            Hs[kz, kt] = resolvent[kz, kt].U*Diagonal(resolvent[kz, kt].S)*resolvent[kz, kt].Vt
        end
    else
        return getfield(resolvent, sym)
    end
end

function resolvent_at_k(kt::Int, kz::Int, dūdy::Vector{T}, ω::T, β::T, Re::T, Ro::T, Dy::Matrix{T}, Dy2::Matrix{T}) where {T}
    # compute wall-normal discretisation size
    Ny = length(dūdy)

    # initialise resolvent matrices
    H_inv = zeros(Complex{T}, 4*Ny, 4*Ny)

    # compute laplacian operator
    Δ = Dy2 - I*(kz*β)^2

    # fill resolvent matrix
    H_inv[1:Ny, 1:Ny] = 1im*kt*ω*I - Δ/Re
    H_inv[1:Ny, (Ny + 1):(2*Ny)] = -Diagonal(dūdy) - I*Ro
    H_inv[(Ny + 1):(2*Ny), 1:Ny] = I*Ro
    H_inv[(Ny + 1):(2*Ny), (Ny + 1):(2*Ny)] = 1im*kt*ω*I - Δ/Re
    H_inv[(Ny + 1):(2*Ny), (3*Ny + 1):end] = Dy
    H_inv[(2*Ny + 1):(3*Ny), (2*Ny + 1):(3*Ny)] = 1im*kt*ω*I - Δ/Re
    H_inv[(2*Ny + 1):(3*Ny), (3*Ny + 1):end] = 1im*kz*β*I
    H_inv[(3*Ny + 1):end, (Ny + 1):(2*Ny)] = -Dy
    H_inv[(3*Ny + 1):end, (2*Ny + 1):(3*Ny)] = -1im*kz*β*I

    # initialise mass matrix
    Z = zeros(Ny, Ny)
    # M = [I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z]

    # apply boundary conditions
    H_inv[1, :] .= 0.0; H_inv[1, 1] = 1.0
    H_inv[Ny:(Ny + 1), :] .= 0.0; H_inv[Ny, Ny] = 1.0; H_inv[Ny + 1, Ny + 1] = 1.0
    H_inv[(2*Ny):(2*Ny + 1), :] .= 0.0; H_inv[2*Ny, 2*Ny] = 1.0; H_inv[2*Ny + 1, 2*Ny + 1] = 1.0
    H_inv[3*Ny, :] .= 0.0; H_inv[3*Ny, 3*Ny] = 1.0
    M[1, :] .= 0.0
    M[Ny:(Ny + 1), :] .= 0.0
    M[(2*Ny):(2*Ny + 1), :] .= 0.0
    M[3*Ny, :] .= 0.0

    # invert resolvent and multiply by mass matrix
    # H = inv(H_inv)*M
    H = inv(H_inv)[:, 1:(3*Ny)]

    return H
end
