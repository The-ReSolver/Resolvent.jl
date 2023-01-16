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
