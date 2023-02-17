# This file contains the utility functions that are useful for the operation
# of the other code in this module.

const DivideAndConquer = LinearAlgebra.DivideAndConquer
const QRIteration = LinearAlgebra.QRIteration
struct Lanczos <: LinearAlgebra.Algorithm end

struct TruncateSVD{TRUNC}; end
struct FullSVD; end

function TruncateSVD(trunc_length::Integer)
    if trunc_length <= 0
        return FullSVD()
    else
        return TruncateSVD{trunc_length}()
    end
end

truncate_svd(svd::SVD, ::FullSVD) = svd
truncate_svd(svds::Matrix{SVD}, ::FullSVD) = svds
function truncate_svd(svd::SVD, ::TruncateSVD{TRUNC}) where {TRUNC}
    U, S, V = svd
    U_trunc = U[:, 1:TRUNC]
    S_trunc = S[1:TRUNC]
    V_trunc = V[:, 1:TRUNC]
    return SVD(U_trunc, S_trunc, V_trunc')
end
function truncate_svd(svds::Matrix{SVD}, ::TruncateSVD{TRUNC}) where {TRUNC}
    svds_truncated = similar(svds)

    for kt in 1:size(svds, 3), kz in 1:size(svds, 2)
        svds_truncated = truncate_svd(svds[kz, kt], TruncateSVD(TRUNC))
    end

    return svds_truncated
end

function LinearAlgebra.cholesky(ws::AbstractVector)
    sqrt_ws = sqrt.(ws)
    chol = Diagonal(sqrt_ws); chol_inv = Diagonal(1 ./ sqrt_ws)
    Z = zeros(length(ws), length(ws))
    L =    [chol Z    Z    Z;
            Z    chol Z    Z;
            Z    Z    chol Z;]
    L_inv =    [chol_inv Z        Z       ;
                Z        chol_inv Z       ;
                Z        Z        chol_inv;]
    return L, L_inv
end

LinearAlgebra.svd(H::Matrix{Complex{T}}, ws::AbstractVector, nvals::Int=1; alg::LinearAlgebra.Algorithm=DivideAndConquer(), debug::Bool=false) where {T<:Real} = _mysvd(H, cholesky(ws)..., nvals, alg, debug)
function _mysvd(H::Matrix{Complex{T}}, L::Matrix{T}, L_inv::Matrix{T}, nvals::Int, alg::LinearAlgebra.Algorithm, debug::Bool) where {T<:Real}
    mySVD = _mysvd(L*H*L_inv, nvals, alg, debug)

    # convert the singular vectors back to the original space
    LinearAlgebra.mul!(mySVD.U, L_inv, copy(mySVD.U))
    LinearAlgebra.mul!(mySVD.Vt, copy(mySVD.Vt), L_inv)

    return mySVD
end
_mysvd(H::Matrix{Complex{T}}, nvals::Int, ::Lanczos, debug::Bool) where {T<:Real} = ((U, S, V) = tsvd(H, nvals; debug=debug); return SVD(U, S, V'))
function _mysvd(H::Matrix{Complex{T}}, nvals::Int, alg::LinearAlgebra.Algorithm, ::Bool) where {T<:Real}
    # compute svd
    SVD = svd(H; alg)

    # loop over singular values until first zero is reached
    truncate_length = length(SVD.S)
    for (i, si) in enumerate(SVD.S)
        abs(si) < 1e-12 && (truncate_length = i - 1; break)
    end
    nvals < truncate_length ? truncate_length = nvals : nothing

    # truncate all the exactly zero singular values and return
    return truncate_svd(SVD, TruncateSVD(truncate_length))
end

function LinearAlgebra.svd(Hs::Matrix{Matrix{Complex{T}}}, L::Matrix{T}, L_inv::Matrix{T}) where {T}
    SVDs = similar(Hs, LinearAlgebra.SVD)

    for kt in 1:size(Hs, 3), kz in 1:size(Hs, 2)
        SVDs[kz, kt] = LinearAlgebra.svd(Hs[kz, kt], L, L_inv)
    end

    return SVDs
end

function resolvent_at_k(kz, kt, dūdy, ω, β, Re, Ro, Dy, Dy2, ::Type{T}=Float64) where {T}
    # compute wall-normal discretisation size
    Ny = length(dūdy)

    # initialise resolvent matrices
    H_inv = zeros(Complex{T}, 4*Ny, 4*Ny)

    # compute laplacian operator
    Δ = Dy2 - I*(kz*β)^2

    # fill resolvent matrix
    # ! the time derivative is negated here for unknown reasons
    H_inv[1:Ny, 1:Ny] = 1im*kt*ω*I - Δ/Re
    H_inv[1:Ny, (Ny + 1):(2*Ny)] = Diagonal(dūdy) - I*Ro
    H_inv[(Ny + 1):(2*Ny), 1:Ny] = I(Ny)*Ro
    H_inv[(Ny + 1):(2*Ny), (Ny + 1):(2*Ny)] = 1im*kt*ω*I - Δ/Re
    H_inv[(Ny + 1):(2*Ny), (3*Ny + 1):end] = Dy
    H_inv[(2*Ny + 1):(3*Ny), (2*Ny + 1):(3*Ny)] = 1im*kt*ω*I - Δ/Re
    H_inv[(2*Ny + 1):(3*Ny), (3*Ny + 1):end] = 1im*kz*β*I(Ny)
    H_inv[(3*Ny + 1):end, (Ny + 1):(2*Ny)] = -Dy
    H_inv[(3*Ny + 1):end, (2*Ny + 1):(3*Ny)] = -1im*kz*β*I(Ny)

    # initialise mass matrix
    Z = zeros(Ny, Ny)
    M =    [I Z Z;
            Z I Z;
            Z Z I;
            Z Z Z]

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
    H = inv(H_inv)*M

    return H
end
