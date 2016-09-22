module QuBase

export init,
    full_matrix,
    exp_value,
    sort_eigs,
    red_rho_1spin


# Functions
# """
#     init(S)
# 
# Initializes calculations by computing Sx, Sy and Sz. All calculations
# are done in the z-basis.
# Returns sparse matrices.
# """
# function init(S::Float64)
#     dim = Int(round(2 * S + 1))
#     S_plus = zeros(dim, dim)
#     for i in range(-S, Int(round(2S)))
#         m = -i - 1
#         ind = Int(i + S + 1)
#         S_plus[ind, ind + 1] = √(S * (S + 1) - m * (m + 1))
#     end
#     S_minus = transpose(S_plus)
# 
#     # Create Sx, Sy and Sz
#     Sx = 0.5 * (S_plus + S_minus)
#     Sy = -0.5im * (S_plus - S_minus)
#     Sz = zeros(dim, dim)
#     for i in S:-1:-S
#         ind = dim - Int(round((i + S)))
#         Sz[ind, ind] = i
#     end
#     return Sx, Sy, Sz
# end
# 
function full_matrix_core(S::SparseMatrixCSC, k::Int64, N::Int64)
    D = size(S)[1]  # Dimensions of the operator/state vector.
    if k == 1
        S_full = kron(S, speye(D^(N - 1)))
    elseif k == 2
        S_full = speye(D)
        S_full = kron(S_full, S)
        S_full = kron(S_full, speye(D^(N - 2)))
    else
        S_full = speye(D^(k - 1))
        S_full = kron(S_full, S)
        S_full = kron(S_full, speye(D^(N - k)))
    end
    return S_full
end

"""
    full_matrix(S, k, N)

Builds the S matrices in an N particle system. Assumes periodic boundary
condition.
"S" could be an operator/state we want to work on. If it is a state, it
must be put in a column vector form. "S" must be sparse.
"k" is the location index of the particle in a particle chain. The first
particle has k=0, the second has k=1 and so on.
Returns a sparse matrix.
"""
function full_matrix(S::Array{Float64, 2}, k::Int64, N::Int64)
    S = sparse(S)
    return full_matrix_core(S, k, N)
end

function full_matrix(S::Array{Complex{Float64}, 2}, k::Int64, N::Int64)
    S = sparse(S)
    return full_matrix_core(S, k, N)
end

function full_matrix(S::SparseMatrixCSC{Float64, Int64}, k::Int64, N::Int64)
    return full_matrix_core(S, k, N)
end

function full_matrix(S::SparseMatrixCSC{Complex{Float64}, Int64}, k::Int64, N::Int64)
    return full_matrix_core(S, k, N)
end

"""
    exp_value(S, psi)

Computes the expected value of an operator with a given state "psi."
"S" is an operator/observable and must be a sparse matrix.
"psi" is a column vector and its sparsity is optional.
"""
function exp_value(S, psi)
    return (ctranspose(psi) * S * psi)[1, 1]
end

"""
    sort_eigs(E, V)

Sort the given set of eigenvectors and eigenvectors by the eigenvalues
from the least to the greatest.
"E" is a list of eigenvalues.
"V" is an eigenvector matrix/array or a list of eigenvectors.
"""
temp = []
function sort_eigs(E::Array, V::Array)
    temp = []
    E_sorted = deepcopy(E)
    for j in 1:size(E)[1]
        for i in 1:size(E)[1]
            if E_sorted[i] > E_sorted[i + 1]
                push!(temp, V_sorted[:, i])
                V_sorted[:, i] = V_sorted[:, i + 1]
                V_sorted[:, i + 1] = pop!(temp)
                push!(temp, E_sorted[i])
                E_sorted[i] = E_sorted[i + 1]
                E_sorted[i + 1] = pop!(temp)
            end
        end
    end
    return E_sorted, V_sorted
end

"""
#    red_rho_1spin(psi, spin)

#Forms a reduced ground state density matrix from a given state "psi."
#Every particle aside from the first will be traced out.
#"psi" must be a column vector. Sparsity is optional. "spin" is the spin
#number. "spin" is optional.
#Returns a sparse matrix.
#"""
#function red_rho_1spin(psi::Array; spin::Float64=0.0)
#    dim = Int(round(2 * spin + 1))
#    if issparse(psi):
#        psi = full(psi)

#    psi_reshaped = 
#end

end
