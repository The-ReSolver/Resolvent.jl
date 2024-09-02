# Definition of the function to compute the resolvent matrix at a specific
# frequency.

function resolvent_at_k(kz, kt, dūdy, Re, Ro, Dy, Dy2, ::Type{T}=Float64) where {T}
    # compute wall-normal discretisation size
    Ny = length(dūdy)

    # initialise resolvent matrices
    H_inv = zeros(Complex{T}, 4*Ny, 4*Ny)

    # compute laplacian operator
    Δ = Dy2 - I*kz^2

    # fill resolvent matrix
    H_inv[1:Ny, 1:Ny] = 1im*kt*I - Δ/Re
    H_inv[1:Ny, (Ny + 1):(2*Ny)] = Diagonal(dūdy) - I*Ro
    H_inv[(Ny + 1):(2*Ny), 1:Ny] = I(Ny)*Ro
    H_inv[(Ny + 1):(2*Ny), (Ny + 1):(2*Ny)] = 1im*kt*I - Δ/Re
    H_inv[(Ny + 1):(2*Ny), (3*Ny + 1):end] = Dy
    H_inv[(2*Ny + 1):(3*Ny), (2*Ny + 1):(3*Ny)] = 1im*kt*I - Δ/Re
    H_inv[(2*Ny + 1):(3*Ny), (3*Ny + 1):end] = 1im*kz*I(Ny)
    H_inv[(3*Ny + 1):end, (Ny + 1):(2*Ny)] = -Dy
    H_inv[(3*Ny + 1):end, (2*Ny + 1):(3*Ny)] = -1im*kz*I(Ny)

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

# TODO: use this in the code to see how it works
# TODO: could rework to do a lot of precomputation
# TODO: do a check to see if a special LU factorisation can be performed
struct Resolvent{Ny}
    H::Array{ComplexF64, 4}
    H_inv::Array{ComplexF64, 4}
    M::Matrix{Float64}
    Δ::Matrix{ComplexF64}
    Dy::Matrix{Float64}
    Dy2::Matrix{Float64}

    function Resolvent(Ny::Int, Dy, Dy2)
        Z = zeros(Ny, Ny)
        M =    [I Z Z;
                Z I Z;
                Z Z I;
                Z Z Z]
        M[1, :] .= 0.0
        M[Ny:(Ny + 1), :] .= 0.0
        M[(2*Ny):(2*Ny + 1), :] .= 0.0
        M[3*Ny, :] .= 0.0
        new{Ny}(zeros(ComplexF64, 4*Ny, 3*Ny), zeros(ComplexF64, 4*Ny, 4*Ny), M, zeros(ComplexF64, Ny, Ny), Dy, Dy2)
    end
end

function (f::Resolvent{Ny})(kz, kt, dūdy, Re, Ro) where {Ny}
    # define laplacian operator
    f.Δ .= f.Dy2 .- I*kz^2

    # fill resolvent matrix
    f.H_inv[1:Ny, 1:Ny] .= 1im*kt*I .- Δ./Re
    f.H_inv[1:Ny, (Ny + 1):(2*Ny)] .= Diagonal(dūdy) .- I.*Ro
    f.H_inv[(Ny + 1):(2*Ny), 1:Ny] .= I(Ny).*Ro
    f.H_inv[(Ny + 1):(2*Ny), (Ny + 1):(2*Ny)] .= 1im.*kt.*I .- Δ./Re
    f.H_inv[(Ny + 1):(2*Ny), (3*Ny + 1):end] .= Dy
    f.H_inv[(2*Ny + 1):(3*Ny), (2*Ny + 1):(3*Ny)] .= 1im.*kt.*I .- Δ./Re
    f.H_inv[(2*Ny + 1):(3*Ny), (3*Ny + 1):end] .= 1im.*kz.*I(Ny)
    f.H_inv[(3*Ny + 1):end, (Ny + 1):(2*Ny)] .= -Dy
    f.H_inv[(3*Ny + 1):end, (2*Ny + 1):(3*Ny)] .= -1im.*kz.*I(Ny)

    # apply boundary conditions
    f.H_inv[1, :] .= 0.0
    f.H_inv[1, 1] = 1.0
    f.H_inv[Ny:(Ny + 1), :] .= 0.0
    f.H_inv[Ny, Ny] = 1.0
    f.H_inv[Ny + 1, Ny + 1] = 1.0
    f.H_inv[(2*Ny):(2*Ny + 1), :] .= 0.0
    f.H_inv[2*Ny, 2*Ny] = 1.0
    f.H_inv[2*Ny + 1, 2*Ny + 1] = 1.0
    f.H_inv[3*Ny, :] .= 0.0
    f.H_inv[3*Ny, 3*Ny] = 1.0

    # invert resolvent and multiply by mass matrix
    f.H .= inv(H_inv)*M
end
