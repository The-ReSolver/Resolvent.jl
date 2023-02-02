# This file contains the definitions required for a general mode implementation
# (probably to moved to the interface package at a later date) which defines
# the basis to express flows upon.

#=
* A mode is simply a function onto which a flow field can be projected onto.

* The goal of any mode type is to provide this basic functionality, not only
* for a single mode, but for collection of modes forming a basis (a simple
* vector would do for this for any flow type).

* This can be expressed abstractly, thus allowing a consistent interface to be
* defined for all the modal analysis techniques used.

* It will also be useful to be able to project modes onto each other, mainly
* for the sake of comparison of the modes generated from different techniques.

* Ultimately this abstract part of this implementation will be moved to the
* interface package.

* Therefore, the following needs to be implemented in this file:
*   - an abstract mode type 🟢
*   - projection methods (for both flow fields and other modes) on this abstract type 
*   - a similar projection method for a collection of modes (inherits from the abstract method so shouldn't need any concrete implementations)
*   - a concrete implementation of this for channel flow (restricts us to modes as a wall-normal profile)
*   - a concrete implementation of the projection method
*   - plotting method for the modes
=#

abstract type AbstractMode{T, N} <: AbstractArray{T, N} end

# array interface stuff
Base.IndexStyle(::Type{<:AbstractMode}) = IndexLinear()
Base.parent(m::AbstractMode) = m.mode
Base.size(m::AbstractMode) = size(parent(m))
Base.getindex(m::AbstractMode, i::Int) = parent(m)[i]
Base.setindex!(m::AbstractMode, v, i::Int) = (parent(m)[i] = v)

function LinearAlgebra.dot(::AbstractMode{<:Any, N}, ::AbstractArray{<:Any, N}) where {N} end
function LinearAlgebra.dot(::AbstractMode{<:Any, N}, ::AbstractArray{<:Any, M}) where {N, M} end

function svd2modes(svd::SVD{T}, ::Type{M}) where {T, M<:AbstractMode{T}}
    # initialise vectors to hold the modes of the decomposition
    U_modes = Vector{M}(undef, size(svd, 1))
    V_modes = Vector{M}(undef, size(svd, 1))

    # loop over the 
    for i in 1:length(svd.S)
        U_modes[i] = M(svd.U[:, i])
        V_modes[i] = M(svd.Vt[i, :])
    end

    return U_modes, svd.S, V_modes
end

"""
    project!(   A::AbstractArray{<:Any, N},
                M::Vector{<:AbstractMode{<:Any, N}},
                U::AbstractArray{<:Any, N}) -> A

Compute the projection of a field onto a set of modes.
"""
function project!(A::AbstractArray{<:Any, N}, M::Vector{<:AbstractMode{<:Any, N}}, U::AbstractArray{<:Any, N}) where {N}
    # project the U onto each mode and assign result to element of A
    for (i, mode) in enumerate(M)
        A[i] = LinearAlgebra.dot(mode, U)
    end

    return A
end

struct ChannelMode{T<:Number} <: AbstractMode{T, 1}
    mode::Vector{T}
end

svd2modes(svd::SVD{T}) where {T} = svd2modes(svd, ChannelMode{T})
