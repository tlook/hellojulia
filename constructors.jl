module constructors

include("./qutypes.jl")
include("./qumath.jl")
importall Qutypes
importall QuMath

export qzeros, qspzeros, qones, qspones, qeye, qspeye,
    splus, sminus, sigmax, sigmay, sigmaz


"""
    qzeros(m, n)

Creates a zero operator.
"""
function qzeros(m::Int64, n::Int64)
    return DenseOp(zeros(m, n))
end


"""
    qspzeros(m, n)

Creates a sparse zero operator.
"""
function qspzeros(m::Int64, n::Int64)
    return SparseOp(spzeros(m, n))
end


"""
    qones(m, n)

Creates a ones operator.
"""
function qones(m::Int64, n::Int64)
    return DenseOp(ones(m, n))
end


"""
    qspones(m, n)

Creates a sparse ones operator.
"""
function qspones(m::Int64, n::Int64)
    return SparseOp(spones(m, n))
end


"""
    qeye(d)

Creates an identity operator.
"""
function qeye(d::Int64)
    return DenseOp(eye(d))
end


"""
    qspeye(d)

Creates a sparse identity operator.
"""
function qspeye(d::Int64)
    return SparseOp(eye(d))
end


"""
    splus(S)

Creates the raising operator (without ħ) in the z basis
"""
function splus(S::Float64)
    dim = Int(round(2 * S + 1))
    Splus = qzeros(dim, dim)
    for i in range(-S, Int(round(2S)))
        m = -i - 1
        ind = Int(i + S + 1)
        Splus[ind, ind + 1] = √(S * (S + 1) - m * (m + 1))
    end
    return Splus
end


"""
    sminus(S)

Creates the lowering operator (without ħ) in the z basis
"""
function sminus(S::Float64)
    return ctranspose(splus(S))
end


"""
    sigmax(S)

Creates the first Pauli spin operator in the z basis
"""
function sigmax(S::Float64)
    Splus = splus(S)
    Sminus = sminus(S)
    return 0.5 * (Splus + Sminus)
end


"""
    sigmay(S)

Creates the second Pauli spin operator in the z basis
"""
function sigmay(S::Float64)
    Splus = splus(S)
    Sminus = sminus(S)
    return -0.5im * (Splus - Sminus)
end


"""
    sigmaz(S)

Creates the third Pauli spin operator in the z basis
"""
function sigmaz(S::Float64)
    Sz = qzeros(dim, dim)
    for i in S:-1:-S
        ind = dim - Int(round((i + S)))
        Sz[ind, ind] = i
    end
    return Sz
end
