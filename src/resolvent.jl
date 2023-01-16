# This file contains the type definition for the Resolvent operator, stored in
# terms of its SVD decomposition.

# TODO: find known test case.
# TODO: test utility functions for the truncation.
# TODO: implement method to reconstruct Resolvent from the SVD object.
# TODO: implement test for Resolvent operator.
# TODO: implement test for it's SVD decomposition.

struct Resolvent{SIZE<:NTuple{3, Int}, TRUNC<:Bool, T<:Real, N<:Int}
    res_svd::Array{SVD, N}

    function Resolvent(Hs::AbstractArray{Matrix{Complex{T}}, N}, trunc::TruncatedSVD) where {T, N}
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

    Resolvent(Hs, TruncatedSVD(trunc))
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
            Hs[kz, kt] = resolvent.res_svd[kz, kt].U*Diagonal(resolvent.res_svd[kz, kt].S)*resolvent.res_svd[kz, kt].Vt
        end
    else
        return getfield(resolvent, sym)
    end
end

# TODO: create proper way to access SVD at particular mode number (getter method or iteration?)
# TODO: redo this function using `getproperty()` method to allow dot syntax to be used (also would that be good for each part of the decomposition)
function get_resolvent(resolvent::Resolvent{SIZE}) where {SIZE}
    # initialise matrix of Resolvents
    Hs = Matrix{Matrix{Complex{T}}}(undef, SIZE[2:3]...)

    # loop over modes numbers computing the Resolvent from its SVD
    for kt in 1:SIZE[3], kz in 1:SIZE[2]
        Hs[kz, kt] = resolvent.res_svd[kz, kt].U*Diagonal(resolvent.res_svd[kz, kt].S)*resolvent.res_svd[kz, kt].Vt
    end

    return Hs
end
