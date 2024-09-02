# Definitions for the custom SVD of the resolvent matrix

const DivideAndConquer = LinearAlgebra.DivideAndConquer
const QRIteration = LinearAlgebra.QRIteration
struct Lanczos <: LinearAlgebra.Algorithm end
struct Arnoldi <: LinearAlgebra.Algorithm end
struct Adaptive <: LinearAlgebra.Algorithm end

# TODO: remove all this complication please!!!
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


nuc_norm(H::Matrix) = tr(sqrt(H'*H))


LinearAlgebra.svd(H::Matrix{Complex{T}}, ws::AbstractVector, nvals::Union{Nothing, Int}=nothing; alg::LinearAlgebra.Algorithm=DivideAndConquer(), debug::Bool=false) where {T<:Real} = _mysvd(H, cholesky(ws)..., nvals, alg, debug)

function _mysvd(H::Matrix{Complex{T}}, L::Matrix{T}, L_inv::Matrix{T}, nvals::Union{Nothing, Int}, ::Adaptive, debug::Bool) where {T<: Real}
    try
        return _mysvd(H, L, L_inv, nvals, DivideAndConquer(), debug)
    catch end

    try
        return _mysvd(H, L, L_inv, nvals, QRIteration(), debug)
    catch end

    try
        return _mysvd(H, L, L_inv, nvals, Lanczos(), debug)
    catch end

    try
        return _mysvd(H, L, L_inv, nvals, Arnoldi(), debug)
    catch end
    
    throw(LinearAlgebra.LAPACKException)
end

function _mysvd(H::Matrix{Complex{T}}, L::Matrix{T}, L_inv::Matrix{T}, nvals::Union{Nothing, Int}, alg::LinearAlgebra.Algorithm, debug::Bool) where {T<:Real}
    # convert nvals to integer if it is a nothing
    nvals === nothing ? nvals = lastindex(H, 2) : nothing

    # compute svd
    mySVD = _mysvd(L*H*L_inv, nvals, alg, debug)

    # convert the singular vectors back to the original space
    LinearAlgebra.mul!(mySVD.U, L_inv, copy(mySVD.U))
    LinearAlgebra.mul!(mySVD.Vt, copy(mySVD.Vt), L_inv)

    return mySVD
end

_mysvd(H::Matrix{Complex{T}}, nvals::Int, ::Lanczos, debug::Bool) where {T<:Real} = ((U, S, V) = tsvd(H, nvals; debug=debug); return SVD(U, S, V'))

_mysvd(H::Matrix{Complex{T}}, nvals::Int, ::Arnoldi, debug::Bool) where {T<:Real} = (svds(H; nsv=nvals, tol=1e-12, maxiter=1000); return SVD(U, S, V'))

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