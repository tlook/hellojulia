module QuMath

# TODO: Test/fix the eig and eigs functions
# TODO: Implement expm and ^

include("./qutypes.jl")

export ⊗

importall QuTypes
import Base: +, -, *, /

# Define common mathematical operations for our custom types
# Addition
function +{T<:AbstractState}(a::T, b::T)
    if length(a) ≠ length(b)
        error("Dimension mismatch: a has length $(length(a)) while b has length $(length(b))")
    else
        return T(a.dat + b.dat)
    end
end

function +{T<:AbstractOperator}(a::T, b::T)
    if size(a) ≠ size(b)
        error("Dimension mismatch: a has dimensions $(size(a)[1])x$(size(a)[2]) while b has dimensions $(size(b)[1])x$(size(b)[2])")
    else
        return DenseOp(a.dat + b.dat)
    end
end

# Subtraction
function -{T<:AbstractState}(a::T, b::T)
    if length(a) ≠ length(b)
        error("Dimension mismatch: a has length $(length(a)) while b has length $(length(b))")
    else
        return T(a.dat - b.dat)
    end
end

function -{T<:AbstractOperator}(a::T, b::T)
    if size(a) ≠ size(b)
        error("Dimension mismatch: a has dimensions $(size(a)[1])x$(size(a)[2]) while b has dimensions $(size(b)[1])x$(size(b)[2])")
    else
        return DenseOp(a.dat - b.dat)
    end
end

# Multiplication
# scalar * QuantumType
*{T<:QuantumType}(a::Number, b::T) = T(a * b.dat)

# QuantumType * scalar
*{T<:QuantumType}(a::T, b::Number) = T(a.dat * b)

# bra * ket = scalar
function *(a::DenseBra, b::DenseKet)
    if length(a) ≠ length(b)
        error("Dimension mismatch: a has length $(length(a)) while b has length $(length(b))")
    else
        return (a.dat * b.dat)[1]
    end
end

function *(a::SparseBra, b::SparseKet)
    if length(a) ≠ length(b)
        error("Dimension mismatch: a has length $(length(a)) while b has length $(length(b))")
    else
        return (a.dat * b.dat)[1]
    end
end

# ket * bra = operator
*(a::DenseKet, b::DenseBra) = DenseOp(a.dat * b.dat)
*(a::SparseKet, b::SparseBra) = SparseOp(a.dat * b.dat)

# operator * ket = ket
function *(a::DenseOp, b::DenseKet)
    if size(a)[2] ≠ length(b)
        error("Dimension mismatch: a has dimensions $(size(a)[1])x$(size(a)[2]) while b has length $(length(b))")
    else
        return DenseKet(a.dat * b.dat)
    end
end

function *(a::SparseOp, b::SparseKet)
    if size(a)[2] ≠ length(b)
        error("Dimension mismatch: a has dimensions $(size(a)[1])x$(size(a)[2]) while b has length $(length(b))")
    else
        return SparseKet(a.dat * b.dat)
    end
end

# bra * operator = bra
function *(a::DenseBra, b::DenseOp)
    if length(a) ≠ size(b)[1]
        error("Dimension mismatch: a has length $(length(a)) while b has dimensions $(size(b)[1])x$(size(b)[2])")
    else
        return DenseBra(a.dat * b.dat)
    end
end

function *(a::SparseBra, b::SparseOp)
    if length(a) ≠ size(b)[1]
        error("Dimension mismatch: a has length $(length(a)) while b has dimensions $(size(b)[1])x$(size(b)[2])")
    else
        return SparseBra(a.dat * b.dat)
    end
end

# operator * operator
function *{T<:AbstractOperator}(a::T, b::T)
    if size(a) ≠ size(b)
        error("Dimension mismatch: a has dimensions $(size(a)[1])x$(size(a)[2]) while b has dimensions $(size(b)[1])x$(size(b)[2])")
    else
        return T(a.dat * b.dat)
    end
end

# Division
# QuantumType / scalar
/{T<:QuantumType}(a::T, b::Number) = T(a.dat / b)

# Kron
⊗{T<:QuantumType}(a::T, b::T) = T(kron(a.dat, b.dat))

# Other functions
Base.norm{T<:QuantumType}(a::T) = norm(a.dat)

# Math functions such as eig
Base.eig(op::DenseOp) = eig(op.dat)
Base.eigs(op::SparseOp) = eigs(op.dat)

end
