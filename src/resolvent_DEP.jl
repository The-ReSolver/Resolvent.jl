# This file contains the type definition for the Resolvent operator, stored in
# terms of its SVD decomposition.

# ! All the stuff in this file is stupid, and a unnecessary complicated layer on
# ! on top of something relatively simple: a finite set of modes with associated
# ! gains.

# ! Deprecate this disaster and just keep it simple (stupid!).



# FIXME: the array-of-arrays approach to this is cursed and should be abandoned!!!
struct Resolvent{SIZE, T, N}
    res_svd::Array{SVD, N}

    function Resolvent(Hs::AbstractArray{Matrix{Complex{T}}, N}, ws, trunc) where {T, N}
        println("2")
        new{(size(Hs[1], 1), size(Hs)...), T, N}(truncate_svd(svd(Hs, cholesky(ws)...), trunc))
    end
end

function Resolvent(U::AbstractArray, ws::AbstractVector, dūdy::Vector, ω::Real, β::Real, Re::Real, Ro::Real, Dy::AbstractMatrix, Dy2::AbstractMatrix, ::Type{T}=Float64; trunc::Int=0) where {T}
    # check the inputs are valid
    (ω >= 0 && β > 0 && Re > 0 && 0 <= Ro <= 1) || throw(ArgumentError("Flow parameters are not valid!"))
    size(Dy) == size(Dy2) == (size(U, 1), size(U, 1)) || throw(ArgumentError("Differentiation matrices are not compatible sizes!"))

    # initialise array to hold resolvents
    # TODO: use similar here instead as it is cleaner and easier to understand what the fuck is going on
    Hs = [Matrix{Complex{T}}(undef, 4*size(U, 1), 3*size(U, 1)) for _ in 1:size(U, 3), _ in 1:size(U, 2)]

    # loop over mode numbers and generate resolvent at each of them
    for kt in 1:size(U, 3), kz in 1:size(U, 2)
        Hs[kz, kt] .= resolvent_at_k(kz, kt, dūdy, ω, β, Re, Ro, Dy, Dy2)
    end

    println("1")
    Resolvent(Hs, ws, TruncateSVD(trunc))
end

# ! probably breaks the `function barrier` approach that makes value types more efficient
function truncate!(H::Resolvent{SIZE}, trunc::Int) where {SIZE}
    for kt in 1:SIZE[3], kz in 1:SIZE[2]
        H[kz, kt] = truncate_svd(H[kz, kt], TruncateSVD(trunc))
    end

    return H
end

# TODO: test these indexing methods
# Base.getindex(resolvent::Resolvent, i::Int) = resolvent.res_svd[i] # ! is this first method required?
Base.getindex(resolvent::Resolvent, I...) = resolvent.res_svd[I...]

function Base.getproperty(resolvent::Resolvent{SIZE}, sym::Symbol) where {SIZE}
    if sym == :H
        # initialise matrix of Resolvents
        # ! it is generally faster to separate out the initialisation of new arrays from the computations
        # ! done on them into distinct functions (allows compiler optimisations), also for `resolvent_at_k`
        Hs = Matrix{Matrix{Complex{T}}}(undef, SIZE[2:3]...)

        # loop over mode numbers computing the Resolvent from its SVD
        for kt in 1:SIZE[3], kz in 1:SIZE[2]
            Hs[kz, kt] = resolvent[kz, kt].U*Diagonal(resolvent[kz, kt].S)*resolvent[kz, kt].Vt
        end
    else
        return getfield(resolvent, sym)
    end
end
