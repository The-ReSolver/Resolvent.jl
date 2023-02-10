# This file contains the definitions required for a general mode implementation
# (probably to moved to the interface package at a later date) which defines
# the basis to express flows upon.

# -----------------------------------------------------------------------------
# Abstract interface for modes and projection
# -----------------------------------------------------------------------------

# To properly create a concrete implementation of the abstract modes interface
# the followsing objects have to be defined:
#   - Subtype of AbstractMode with fields for the mode data and the quadrature
#       weights;
#   - LinearAlgebra.dot method for projection of your desired field and a single
#       mode;

# The following are optional:
#   - weights method to extract the quadrature weights array. Default behaviour
#       assumes that the weights are stored as a field of the modes, but
#       overloading this method allows for other approaches.
#   - ModeDimension constructor to determine which dimension to modes operate
#       on

abstract type AbstractMode{N, T} <: AbstractArray{T, N} end

# array interface stuff
Base.IndexStyle(::Type{<:AbstractMode}) = IndexLinear()
Base.parent(m::AbstractMode) = m.mode
Base.size(m::AbstractMode) = size(parent(m))
Base.getindex(m::AbstractMode, i::Int) = parent(m)[i]
Base.setindex!(m::AbstractMode, v, i::Int) = (parent(m)[i] = v)

weights(M::AbstractMode) = M.ws

# * projection of single profile onto single mode
LinearAlgebra.dot(M::AbstractMode{N}, U::AbstractArray{<:Any, N}) where {N} = throw(MethodError(LinearAlgebra.dot, (M, U)))

# * projection of single profile onto set of modes (basis functions)
project!(A::AbstractVector, M::AbstractVector{<:AbstractMode{N}}, U::AbstractArray{<:Any, N}) where {N} = broadcast!(mode->dot(mode, U), A, M)

# * projection of a profile for every wavenumber pair onto a set of modes for each wavenumber pair
function project!(A::AbstractArray, modes::AbstractArray{M}, U::AbstractArray) where {N, M<:AbstractVector{<:AbstractMode{N}}}
    @views for I in CartesianIndices(size(A)[2:end])
        project!(A[:, I], modes[I], U[[Colon() for _ in 1:N]..., I])
    end

    return A
end

function test_modeinterface(mode::AbstractMode)
    # test the weights extraction is well defined
    try
        weights(mode)
    catch e
        throw(e)
    end

    # test the dot product produces a well defined output
    @assert dot(mode, similar(parent(mode))) isa Number "Inner-product does not output single number as output!"
end

# * the only way to improve this is to provide some way of knowing which dimensions
# * to slice over without assumption, maybe this could be stored as some sort of ModeDimension method
# struct ModeDimensions{N}; dims::NTuple{N, Int}; end
# function project!(A::AbstractArray, modes::AbstractArray{M}, U::AbstractArray{<:Any, NU}, dims::ModeDimensions{NU}=ModeDimensions(N, NU)) where {NU, N, M<:AbstractVector{<:AbstractMode{N}}} end
# * The point of the above code is to provide a fallback behaviour for the projection
# * that assumed the dimensions over which the mode applies are the first `N` dimensions
# * (so the first `N` dimensions should be sliced and the rest indexed). By overloading
# * the construction of the ModeDimension object for our concrete type we will be able
# * to modify this method so that it can slice over any number and order of dimensions
# * as desired.

# -----------------------------------------------------------------------------
# Concrete implementation for channel modes
# -----------------------------------------------------------------------------

struct ChannelMode{T<:Number} <: AbstractMode{1, T}
    mode::Vector{T}
    ws::Vector{Float64}
end

function svd2channelmodes(svd::SVD{T}, ws::AbstractArray{<:AbstractFloat, N}) where {T, N}
    # initialise vectors to hold the modes of the decomposition
    U_modes = Vector{ChannelMode{T}}(undef, size(svd, 1))
    V_modes = Vector{ChannelMode{T}}(undef, size(svd, 1))

    # loop over the 
    for i in 1:length(svd.S)
        U_modes[i] = ChannelMode{T}(svd.U[:, i], ws)
        V_modes[i] = ChannelMode{T}(svd.Vt[i, :], ws)
    end

    return U_modes, svd.S, V_modes
end

function LinearAlgebra.dot(mode::ChannelMode, U::AbstractVector)
    sum = 0
    for ny in eachindex(mode)
        sum += weights(mode)[ny]*mode[ny]*U[ny]
    end

    return sum
end

# TODO: test the above projection methods for my concrete implementation one-by-one
