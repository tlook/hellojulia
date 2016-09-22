# TODO

## quantumtypes.jl
- [x] Create basic types
    - [x] DenseOp
    - [x] SparseOp
    - [x] DenseBra
    - [x] SparseBra
    - [x] DenseKet
    - [x] SparseKet
    
    Food for thought: Do bras and kets need to be sparse on any occasion? How would that change affect everything else?
- [ ] Data fields contained in QuTypes
    - [x] Data array
    - [ ] System length
    - [ ] Spin number
    - [ ] Update functions to work with system length and spin number
- [x] Access raw data contained in QuTypes
    - [x] rawdata 
- [x] Ability to convert between sparse and dense vectors/operators
    - [x] Base.sparse
    - [x] Base.full
    - [x] Getting the type constructors to do conversions from one type to another
- [x] Ability to convert a bra to a ket and vice versa
    - [x] ctranspose
- [ ] Get the following functions to work:
    - [x] Base.length
    - [x] Base.size
    - [x] Base.ndims
    - [x] Base.issparse
    - [ ] Base.getindex
    - [x] Base.setindex!
    - [x] Base.copy
    - [x] Base.summary
    - [x] Base.show
- [ ] Create promotion functions/definitions
    - [ ] Scenario 1: For an operation between two objects, one object contains a different number type than the other. Int ⇒ Float ⇒ Complex
    - [ ] Scenario 2: (Maybe) For a basic mathematical operation (+. -. *, /) between a QuType and a built-in array like type, promote the built-in array type to a QuType

Note: the current getindex doesn't fully work as it returns slices in arrays/sparse matrices instead of one of the QuTypes

## qmmath.jl
- [ ] Create basic arithmetic functions
    - [x] +
    - [x] -
    - [x] *
    - [x] /
    - [ ] ^
    - [x] ⊗
    - [x] Correct error messages to be displayed with dimension mismatch
- [ ] Element wise functions for non-QuTypes x QuTypes
    - [ ] .*
    - [ ] ./
- [ ] Create basic functions specific to the work on spin systems
    - [ ] eig
    - [ ] eigs
    - [ ] eigvals
    - [ ] eigvecs
    - [ ] eigmax
    - [ ] eigmin
    - [ ] expm
- [x] Other essential functions
    - [x] norm

## constructors.jl
- [x] Versions of built-in constructors that generate their corresponding QuTypes directly
    - [x] qzeros
    - [x] qspzeros
    - [x] qones
    - [x] qspones
    - [x] qeye
    - [x] qspeye
- [x] Constructors of basic spin system operators
    - [x] splus
    - [x] sminus
    - [x] sigmax
    - [x] sigmay
    - [x] sigmaz

## quantumbase.jl
- [ ] (List of functions to be completed)
    - [ ] fullmatrix
    - [ ] expval
    - [ ] sorteigs
    - [ ] redrho1spin

Note: all of the above (pre-existing) functions need to be retooled to work with QuTypes
