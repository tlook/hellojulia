# TODO: Make returns of slicing conform to my type
# TODO: Create promotion functions/definitions

export QuantumType, AbstractState, AbstractOperator,
    Bra, Ket, DenseOp, SparseOp,
    rawdata

# Abstract type definitions
abstract QuantumType
abstract AbstractState <: QuantumType
abstract AbstractOperator <: QuantumType

# Concrete type definitions
type DenseOp{T<:Number} <: AbstractOperator
    dat::Array{T, 2}
    DenseOp(A) = new(A)
end

type SparseOp{T<:Number} <: AbstractOperator
    dat::SparseMatrixCSC{T, Int64}
    SparseOp(A) = new(A)
end

type Bra{T<:Number} <: AbstractState
    dat::Array{T, 2}
    function Bra(ψ)
        if length(size(ψ)) ≠ 2 || size(ψ)[1] ≠ 1
            error("A Bra must be a row vector.")
        end
        new(ψ)
    end
end

type Ket{T<:Number} <: AbstractState
    dat::Array{T, 1}
    function Ket(ψ)
        if length(size(ψ)) ≠ 1
            error("A DenseKet must be a column vector.")
        end
        new(ψ)
    end
end

# Constructor functions
SparseOp{T<:Number}(mat::SparseMatrixCSC{T, Int64}) = SparseOp{T}(mat)
DenseOp{T<:Number}(mat::Array{T, 2}) = DenseOp{T}(mat)
Bra{T<:Number}(vec::Array{T, 2}) = Bra{T}(vec)
Ket{T<:Number}(vec::Array{T, 1}) = Ket{T}(vec)

# Functions regarding sparsity
Base.sparse(a::DenseOp) = SparseOp(sparse(a.dat))
Base.full(a::SparseOp) = DenseOp(full(a.dat))
SparseOp(a::DenseOp) = sparse(a)
DenseOp(a::SparseOp) = full(a)

# Convert Bra to Ket and vice versa
Base.ctranspose{T<:AbstractOperator}(a::T) = T(ctranspose(a.dat))
Base.ctranspose(a::Bra) = Ket(ctranspose(a.dat)[:, 1])
Base.ctranspose(a::Ket) = Bra(ctranspose(a.dat))

# Other redefinitions:bn
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
    print(io, summary(X))
    println(io, ":")
    print(io, X.dat)
end

rawdata{T<:QuantumType}(a::T) = a.dat

# Promotion rules
#convert{T<:Number}(::
