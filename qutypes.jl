module QuTypes

# TODO: Make returns of slicing conform to my type
# TODO: Create promotion functions/definitions

export QuantumType, AbstractState, AbstractOperator,
    Bra, Ket,
    DenseOp, DenseBra, DenseKet,
    SparseOp, SparseBra, SparseKet,
    rawdata

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
SparseKet{T<:Number}(vec::SparseVector{T, Int64}) = SparseKet{T}(vec)
DenseKet{T<:Number}(vec::Array{T, 1}) = DenseKet{T}(vec)

# Functions regarding sparsity
Base.sparse(a::DenseBra) = SparseBra(sparse(a.dat))
Base.sparse(a::DenseKet) = SparseKet(sparse(a.dat))
Base.sparse(a::DenseOp) = SparseOp(sparse(a.dat))
Base.full(a::SparseBra) = DenseBra(full(a.dat))
Base.full(a::SparseKet) = DenseKet(full(a.dat))
Base.full(a::SparseOp) = DenseOp(full(a.dat))
SparseOp(a::DenseOp) = sparse(a)
SparseBra(a::DenseBra) = sparse(a)
SparseKet(a::DenseKet) = sparse(a)
DenseOp(a::SparseOp) = full(a)
DenseBra(a::SparseBra) = full(a)
DenseKet(a::SparseKet) = full(a)

# Convert Bra to Ket and vice versa
Base.ctranspose{T<:AbstractOperator}(a::T) = T(ctranspose(a.dat))
Base.ctranspose(a::DenseBra) = DenseKet(ctranspose(a.dat))
Base.ctranspose(a::SparseBra) = SparseKet(ctranspose(a.dat))
Base.ctranspose(a::DenseKet) = DenseBra(ctranspose(a.dat))
Base.ctranspose(a::SparseKet) = SparseBra(ctranspose(a.dat))

# Other redefinitions
Base.length{T<:QuantumType}(a::T) = length(a.dat)
Base.size{T<:QuantumType}(a::T) = size(a.dat)
Base.ndims{T<:QuantumType}(a::T) = ndims(a.dat)
Base.issparse{T<:QuantumType}(a::T) = issparse(a.dat)
Base.in{T<:QuantumType}(a::T) = in(a.dat)
Base.:(==){T<:QuantumType}(a::T, b::T) = a.dat==b.dat # Add bases comparison
Base.getindex{T<:QuantumType}(a::T, ind::Int...) = getindex(a.dat, ind...)
Base.getindex{T<:QuantumType}(a::T, ind::Range...) = getindex(a.dat, ind...)
Base.getindex{T<:QuantumType}(a::T, ind::Vector...) = getindex(a.dat, ind...)
Base.setindex!{T<:QuantumType}(a::T, X::Number, ind::Int...) = setindex!(a.dat, X, ind...)
Base.setindex!{T<:QuantumType}(a::T, X::T, ind::Range...) = setindex!(a.dat, X.dat, ind...)
Base.setindex!{T<:QuantumType}(a::T, X::T, ind::Vector...) = setindex!(a.dat, X.dat, ind...)
Base.copy{T<:QuantumType}(a::T) = T(copy(a.dat))
Base.summary{T<:AbstractOperator}(a::T) = "$(size(a.dat)[1])x$(size(a.dat)[2]) $T"
Base.summary{T<:Bra}(a::T) = "Length $(size(a.dat)[2]) $T"
Base.summary{T<:Ket}(a::T) = "Length $(size(a.dat)[1]) $T"
function Base.show(io::IO, X::QuantumType)
    println(io, summary(X))
    print(io, X.dat)
end
rawdata{T<:QuantumType}(a::T) = a.dat

# Promotion rules

end
