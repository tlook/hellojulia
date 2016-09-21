module QuTypes

export QuantumType, AbstractState, AbstractOperator,
    Bra, Ket,
    DenseOp, DenseBra, DenseKet,
    SparseOp, SparseBra, SparseKet

# Abstract type definitions
abstract QuantumType
abstract AbstractState <: QuantumType
abstract AbstractOperator <: QuantumType

abstract Bra <: AbstractState
abstract Ket <: AbstractState

# Concrete type definitions
type DenseOp{T<:Number} <: AbstractOperator
    dat::Array{T, 2}
    DenseOp(A) = new(A)
end

type SparseOp{T<:Number} <: AbstractOperator
    dat::SparseMatrixCSC{T, Int64}
    SparseOp(A) = new(A)
end

type DenseBra{T<:Number} <: Bra
    dat::Array{T, 2}
    function DenseBra(ψ)
        if length(size(ψ)) ≠ 2 || size(ψ)[1] ≠ 1
            error("A DenseBra must be a row vector.")
        end
        new(ψ)
    end
end

type SparseBra{T<:Number} <: Bra
    dat::SparseMatrixCSC{T, Int64}
    function SparseBra(ψ)
        if length(size(ψ)) ≠ 2 || size(ψ)[1] ≠ 1
            error("A SparseBra must be a row vector.")
        end
        new(ψ)
    end
end

type DenseKet{T<:Number} <: Ket
    dat::Array{T, 1}
    function DenseKet(ψ)
        if length(size(ψ)) ≠ 1
            error("A DenseKet must be a column vector.")
        end
        new(ψ)
    end
end

type SparseKet{T<:Number} <: Ket
    dat::SparseMatrixCSC{T, Int64}
    function SparseKet(ψ)
        if length(size(ψ)) ≠ 1
            error("A SparseKet must be a column vector.")
        end
        new(ψ)
    end
end

# Constructor functions
SparseOp{T<:Number}(mat::SparseMatrixCSC{T, Int64}) = SparseOp{T}(mat)
DenseOp{T<:Number}(mat::Array{T, 2}) = DenseOp{T}(mat)
SparseBra{T<:Number}(vec::SparseMatrixCSC{T, Int64}) = SparseBra{T}(vec)
DenseBra{T<:Number}(vec::Array{T, 2}) = DenseBra{T}(vec)
SparseKet{T<:Number}(vec::SparseMatrixCSC{T, Int64}) = SparseKet{T}(vec)
DenseKet{T<:Number}(vec::Array{T, 1}) = DenseKet{T}(vec)
end
