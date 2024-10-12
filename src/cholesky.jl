# small file to hold Cholesky decomposition functionality

function LinearAlgebra.cholesky(ws::AbstractVector)
    sqrt_ws = sqrt.(ws)
    chol = Diagonal(sqrt_ws); chol_inv = Diagonal(1 ./ sqrt_ws)
    Z = zeros(length(ws), length(ws))
    L =    [chol Z    Z    Z;
            Z    chol Z    Z;
            Z    Z    chol Z;]
    L_inv =    [chol_inv Z        Z       ;
                Z        chol_inv Z       ;
                Z        Z        chol_inv;]
    return L, L_inv
end
