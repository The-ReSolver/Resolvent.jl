# Definition of the function to compute the resolvent matrix at a specific
# frequency.

# TODO: try replacing with in-place operations and views wherever possible in this function
# TODO: make an in-place method so that arrays to store the results can just be reused
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
    # TODO: replace this with a simple slice of the inverse?
    H = inv(H_inv)*M

    return H
end