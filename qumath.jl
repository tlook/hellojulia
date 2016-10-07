export ⊗

import Base: +, -, *, /, ^

# TODO: Test/fix the eig and eigs functions
# TODO: Implement expm

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
*{sT<:Number, T<:QuantumType}(a::sT, b::T) = T{sT}(a * b.dat)

# QuantumType * scalar
*{T<:QuantumType}(a::T, b::Number) = T(a.dat * b)

# bra * ket = scalar
function *(a::Bra, b::Ket)
    if length(a) ≠ length(b)
        error("Dimension mismatch: a has length $(length(a)) while b has length $(length(b))")
    else
        return (a.dat * b.dat)[1]
    end
end

# ket * bra = operator
*(a::Ket, b::Bra) = DenseOp(a.dat * b.dat)

# operator * ket = ket
function *(a::AbstractOperator, b::Ket)
    if size(a)[2] ≠ length(b)
        error("Dimension mismatch: a has dimensions $(size(a)[1])x$(size(a)[2]) while b has length $(length(b))")
    else
        return Ket(a.dat * b.dat)
    end
end

# bra * operator = bra
function *(a::Bra, b::AbstractOperator)
    if length(a) ≠ size(b)[1]
        error("Dimension mismatch: a has length $(length(a)) while b has dimensions $(size(b)[1])x$(size(b)[2])")
    else
        return Bra(a.dat * b.dat)
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

# Power
function ^{T<:AbstractOperator}(a::T, n::Number)
    if size(a)[1] ≠ size(a)[2]
        error("Not a square operator")
    else
        return T(a.dat^n)
    end
end

# Kron
⊗{T<:QuantumType}(a::T, b::T) = T(kron(a.dat, b.dat))

# Other functions
Base.norm{T<:QuantumType}(a::T) = norm(a.dat)

# Math functions such as eig
Base.eig(op::DenseOp) = eig(op.dat)
Base.eigs(op::SparseOp) = eigs(op.dat)
