# Definitions for the custom SVD of the resolvent matrix

const DivideAndConquer = LinearAlgebra.DivideAndConquer
const QRIteration = LinearAlgebra.QRIteration
struct Lanczos <: LinearAlgebra.Algorithm end
struct Arnoldi <: LinearAlgebra.Algorithm end
struct Adaptive <: LinearAlgebra.Algorithm end



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


truncateSVD(svd::SVD, trunc::Int) = SVD(svd.U[:, 1:trunc], svd.S[1:trunc], svd.V[:, 1:trunc]')
