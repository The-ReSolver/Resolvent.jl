# This file contains the utility functions that are useful for the operation
# of the other code in this module.

struct TruncateSVD{TRUNC}; end
struct FullSVD; end

function TruncateSVD(trunc_length::T) where {T<:Integer}
    if trunc_length <= 0
        return FullSVD()
    else
        return TruncateSVD{trunc_length}()
    end
end

truncate_svd!(svd::SVD, ::FullSVD) = svd
# FIXME: this should be done without having to generate a new SVD object
function truncate_svd!(svd::SVD, ::TruncateSVD{TRUNC}) where {TRUNC}
    U, S, V = svd
    U_trunc = U[:, 1:TRUNC]
    S_trunc = S[1:TRUNC]
    V_trunc = V[:, 1:TRUNC]
    return SVD(U_trunc, S_trunc, V_trunc')
end

function LinearAlgebra.cholesky(ws::AbstractVector)
    chol = Diagonal(sqrt.(ws))
    Z = zeros(length(ws), length(ws))
    return [chol Z Z Z; Z chol Z Z; Z Z chol Z; Z Z Z Z]
end

# ! the in-place operations may not play nice
function LinearAlgebra.svd(H::AbstractMatrix{Complex{T}}, L::AbstractMatrix{T}) where {T<:Real}
    L_inv = inv(L)
    SVD = LinearAlgebra.svd(L*H*L_inv)
    LinearAlgebra.mul!(SVD.U, L_inv, SVD.U)
    LinearAlgebra.mul!(SVD.Vt, SVD.Vt, L_inv)
    return SVD
end
