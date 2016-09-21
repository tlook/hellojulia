module QMMath

# TODO: Implement issparse, full, sparse, print, println
# TODO: Ability to convert sparse types to dense types and vice versa
#       without resorting to full and sparse.

include("./quantumtypes.jl")

export ⊗, rawdata

importall QuTypes
import Base: *, +, -

# Define common mathematical operations for our custom types
# Addition
+(a::DenseKet, b::DenseKet) = DenseKet(a.dat + b.dat)
+(a::SparseKet, b::SparseKet) = SparseKet(a.dat + b.dat)
+(a::DenseBra, b::DenseBra) = DenseBra(a.dat + b.dat)

function +(a::DenseOp, b::DenseOp)
    if size(a.dat) ≠ size(b.dat)
        error("Dimension mismatch")
    else
        return DenseOp(a.dat + b.dat)
    end
end

function +(a::SparseOp, b::SparseOp)
    if size(a.dat) ≠ size(b.dat)
        error("Dimension mismatch")
    else
        return SparseOp(a.dat + b.dat)
    end
end

# Subtraction
-(a::DenseKet, b::DenseKet) = DenseKet(a.dat - b.dat)
-(a::SparseKet, b::SparseKet) = SparseKet(a.dat - b.dat)
-(a::DenseBra, b::DenseBra) = DenseBra(a.dat - b.dat)

function -(a::DenseOp, b::DenseOp)
    if size(a.dat) ≠ size(b.dat)
        error("Dimension mismatch")
    else
        return DenseOp(a.dat - b.dat)
    end
end

function -(a::SparseOp, b::SparseOp)
    if size(a.dat) ≠ size(b.dat)
        error("Dimension mismatch")
    else
        return SparseOp(a.dat - b.dat)
    end
end

# Multiplication
# bra * ket = scalar
*(a::DenseBra, b::DenseKet) = a.dat * b.dat
*(a::SparseBra, b::SparseKet) = a.dat * b.dat

# ket * bra = operator
*(a::DenseKet, b::DenseBra) = DenseOp(a.dat * b.dat)
*(a::SparseKet, b::SparseBra) = SparseOp(a.dat * b.dat)

# operator * ket = ket
*(a::DenseOp, b::DenseKet) = DenseKet(a.dat * b.dat)
*(a::SparseOp, b::SparseKet) = SparseKet(a.dat * b.dat)

# bra * operator = bra
*(a::DenseBra, b::DenseOp) = DenseBra(a.dat * b.dat)
*(a::SparseBra, b::SparseOp) = SparseBra(a.dat * b.dat)

# Kron
⊗(a::DenseBra, b::DenseBra) = DenseBra(kron(a.dat, b.dat))
⊗(a::SparseBra, b::SparseBra) = SparseBra(kron(a.dat, b.dat))
⊗(a::DenseKet, b::DenseKet) = DenseKet(kron(a.dat, b.dat))
⊗(a::SparseKet, b::SparseKet) = SparseKet(kron(a.dat, b.dat))
⊗(a::DenseOp, b::DenseOp) = DenseOp(kron(a.dat, b.dat))
⊗(a::SparseOp, b::SparseOp) = SparseOp(kron(a.dat, b.dat))

# Other functions
Base.ctranspose{T<:AbstractOperator}(a::T) = T(ctranspose(a.dat))
Base.ctranspose(a::DenseBra) = DenseKet(ctranspose(a.dat))
Base.ctranspose(a::SparseBra) = SparseKet(ctranspose(a.dat))
Base.ctranspose(a::DenseKet) = DenseBra(ctranspose(a.dat))
Base.ctranspose(a::SparseKet) = SparseBra(ctranspose(a.dat))
Base.length{T<:QuantumType}(a::T) = length(a.dat)
Base.size{T<:QuantumType}(a::T) = size(a.dat)
Base.ndims{T<:QuantumType}(a::T) = ndims(a.dat)
Base.norm{T<:QuantumType}(a::T) = norm(a.dat)
Base.in{T<:QuantumType}(a::T) = in(a.dat)
Base.:(==){T<:QuantumType}(a::T, b::T) = a.dat==b.dat # Add bases comparison
Base.copy{T<:QuantumType}(a::T) = T(copy(a.dat))
rawdata{T<:QuantumType}(a::T) = a.dat

# Math functions such as eig
Base.eig(op::DenseOp) = eig(op.dat)
Base.eigs(op::SparseOp) = eigs(op.dat)

end
