include("../qubase.jl")

using QuBase

export full_hamiltonian_OB,
    full_hamiltonian_PB,
    nearest_chain_hamiltonian_OB,
    nearest_chain_hamiltonian_PB,
    chain_hamiltonian_OB,
    chain_hamiltonian_PB

"""
    full_operators(N)

Creates lists of full operators for the Hamiltonian.
"""
function full_operators(N)
    Sx, Sy, Sz = sigmax(N), sigmay(N), sigmaz(N)
    Sx_ls = []
    Sy_ls = []
    Sz_ls = []
    for l in 1:N
        push!(Sx_ls, full_matrix(Sx, l, N))
        push!(Sy_ls, full_matrix(Sy, l, N))
        push!(Sz_ls, full_matrix(Sz, l, N))
    end
    return Sx_ls, Sy_ls, Sz_ls
end


"""
    chain_hamiltonaian_OB(N, J, delta, n, chain)

Chain Hamiltonian under open boundary conditions.

"n" is the distance between neighbors. For the chain corresponding to
nearest neighbor interaction, "n" is 1 and for the chain corresponding
to second neighbor interaction, "n" is 2.

"chain" is the upper or lower chain of the system. 1 for the upper chain
and 2 for the lower chain.
"""
function chain_hamiltonian_OB(N::Int64, J::Float64, delta::Float64,
                              n::Int64, chain::Int64)
    Sx_ls, Sy_ls, Sz_ls = full_operators(N)
    H = spzeros(2^N, 2^N)
    for l in chain:n:N - n
        Sxx = Sx_ls[l] * Sx_ls[l + n]
        Syy = Sy_ls[l] * Sy_ls[l + n]
        Szz = Sz_ls[l] * Sz_ls[l + n]
        H += J * (Sxx + Syy + delta * Szz)
    end
    return H
end


"""
    chain_hamiltonaian_PB(N, J, delta, n, chain)

Chain Hamiltonian under periodic boundary conditions.

"n" is the distance between neighbors. For the chain corresponding to
nearest neighbor interaction, "n" is 1 and for the chain corresponding
to second neighbor interaction, "n" is 2.

"chain" is the upper or lower chain of the system. 1 for the upper chain
and 2 for the lower chain.
"""
function chain_hamiltonian_PB(N::Int64, J::Float64, delta::Float64,
                              n::Int64, chain::Int64)
    Sx_ls, Sy_ls, Sz_ls = full_operators(N)
    H = spzeros(2^N, 2^N)
    for l in chain:n:N
        if l + n > N
            ln = l + n - N
        else
            ln = l + n
        end
        Sxx = Sx_ls[l] * Sx_ls[ln]
        Syy = Sy_ls[l] * Sy_ls[ln]
        Szz = Sz_ls[l] * Sz_ls[ln]
        H += J * (Sxx + Syy + delta * Szz)
    end
    return H
end


function nearest_chain_hamiltonian_OB(N::Int64, J1::Float64, delta::Float64)
    return chain_hamiltonian_OB(N, J1, delta, 1, 1)
end


function nearest_chain_hamiltonian_PB(N::Int64, J1::Float64, delta::Float64)
    return chain_hamiltonian_CB(N, J1, delta, 1, 1)
end


function second_nearest_neighbor_chain_hamiltonian_OB(N::Int64, J2::Float64,
                                                  delta::Float64, chain::Int64)
    return chain_hamiltonian_OB(N, J2, delta, 2, chain)
end


function second_nearest_neighbor_chain_hamiltonian_PB(N::Int64, J2::Float64,
                                                  delta::Float64, chain::Int64)
    return chain_hamiltonian_PB(N, J2, delta, 2, chain)
end


"""
    full_hamiltonian_OB(N, J1, J2, delta)

H = Σ^2_(n=1) Σ_l J_n (S^x_l S^x_(l+n) + S^y_l S^y_(l+n) + ΔS^z_l S^z_(l+n))
Open boundary conditions.
"""
function full_hamiltonian_OB(N::Int64, J1::Float64, J2::Float64, delta::Float64)
    zigzag_H = nearest_chain_hamiltonian_OB(N, J1, delta)
    upper_leg = second_nearest_neighbor_chain_hamiltonian_OB(N, J2, delta, 2)
    lower_leg = second_nearest_neighbor_chain_hamiltonian_OB(N, J2, delta, 1)
    return zigzag_H + upper_leg + lower_leg
end


"""
    full_hamiltonian_PB(N, J1, J2, delta)

H = Σ^2_(n=1) Σ_l J_n (S^x_l S^x_(l+n) + S^y_l S^y_(l+n) + ΔS^z_l S^z_(l+n))
Periodic boundary conditions.
"""
function full_hamiltonian_PB(N::Int64, J1::Float64, J2::Float64, delta::Float64)
    zigzag_H = nearest_chain_hamiltonian_PB(N, J1, delta)
    upper_leg = second_nearest_neighbor_chain_hamiltonian_PB(N, J2, delta, 2)
    lower_leg = second_nearest_neighbor_chain_hamiltonian_PB(N, J2, delta, 1)
    return zigzag_H + upper_leg + lower_leg
end
