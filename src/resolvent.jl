# This file contains the type definition for the Resolvent operator, stored in
# terms of its SVD decomposition.

struct Resolvent{SIZE<:NTuple{3, Int}, TRUNC<:Bool, T<:Real, N<:Int}
    res_svd::Array{SVD, N}

    function Resolvent(Hs::AbstractArray{Matrix{Complex{T}}, N}, ws::AbstractVector{T}, trunc::Union{TruncateSVD, FullSVD}) where {T, N}
        new{(size(H[s1], 1), size(Hs)...), typeof(trunc), T, N}(truncate_svd!(svd.(Hs, cholesky(ws)), trunc))
    end

    # function Resolvent(Hs::AbstractArray{Matrix{Complex{T}}, N}, ws::AbstractVector{T}, trunc::FullSVD) where {T, N}
    #     new{(size(H[s1], 1), size(Hs)...), false, T, N}(truncate_svd!(svd.(Hs, cholesky(ws)), trunc))
    # end
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
