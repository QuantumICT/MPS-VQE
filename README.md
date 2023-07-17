# MPS-VQE

## Code
* QuantumSpins is a MPS based simulator library written in pure Julia, including quantum gates, quantum circuit, circuit evolution, MPS simulator and other interfaces.
* vqe.jl contains VQE calculation based on MPI parallelization.
* vqe_distributed.jl contains VQE calculation based on Distributed.jl.

## Examples
The parameters of H2, F2, Li2, N2 and other molecule systems as VQE input are provided. param_* stands for the input file of vqe.jl, and QASM based quantum circuits are inputs of vqe_distributed.jl.
To run vqe.jl, 
```
mpiexecjl -n 2 julia vqe.jl param_uccsd_sym_h2_sto3g_1.0
```
To run vqe_distributed.jl,
```
julia -p 2 vqe_distributed.jl
```