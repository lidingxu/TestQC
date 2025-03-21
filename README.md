# Quantum correlations

This is a package of several function allows one to determine if a quantum state is seprable (Quantum Separability) or nonsteerable (Quantum Steerability). 

## Download

The repo contains the submodule kvant which contains some important states and methods for standard quantum information theory routines. Clone inclusively the module as follows:

```bash
$ git clone --recurse-submodules -j8 https://gitlab.com/cn611340/quantum-correlations.git
```

If you forget --recurse-submodules, the submodule is not loaded automatically. 
Then run

```bash
$git submodule update --init
$cd lib/kvant 
$git checkout master 
```

Some dependencies which have to be downloaded possibly are the Julia packages LinearAlgebra, Convex, Combinatorics, MathOptInterface. The Julia package Mosek as solver for convex optimisation problems is also required in the default option. Different solvers that are compatible with Convex.jl can also be used as explained below.

# 1 - Quantum Separability

Method for separability certification for low dimensional system with the polytope adaption technique.
Explanation how the algorithm works is given [here.](https://arxiv.org/abs/2210.10054)



## Startup

Open your julia REPL, navigate to the workspace and load the modules. 

```julia
cd("<YOUR_PATH_TO_REPOSITORY>/quantum-correlations/lab/entanglement/")
include("startup.jl")
using MultiStates
using Entanglement
```
## Minimal example

Consider the standard 2 qubit Bell-state: $$ \ket{\Phi^+} = \frac{1}{\sqrt{2}}(\ket{00} + \ket{11})$$

We can implement this state by first providing its coefficient matrix in the computational basis followed by producing an object of type `MultiState(mat::Matrix{Complex{Float64}}, dims::Array{Int,1})` for which the local dimensions are also specified within the attribute `dims`.

```julia
phi_plus = 1/2 * [1.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 1.0]

>4×4 Matrix{Float64}:
 0.5  0.0  0.0  0.5
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.5  0.0  0.0  0.5

phi_state = MultiState(phi_plus, [2, 2])

>MultiState(ComplexF64[0.5 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.5 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.5 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.5 + 0.0im], [2, 2])

phi_state.mat

>4×4 Matrix{ComplexF64}:
 0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.5+0.0im  0.0+0.0im  0.0+0.0im  0.5+0.0im

phi_state.dims

>2-element Vector{Int64}:
 2
 2

```
The function `EntanglementRobustness` takes as input a MultiState $\rho^{AB}$ and returns a lower bound for the visibility $\chi(\rho^{AB})$ defined by

$$ \chi(\rho^{AB}) = \max\{t: t\rho^{AB} + (1-t)\mathbb{1}/d_{A}d_{B} \in \textnormal{SEP}\}$$

where SEP is the set of separable states.


In Julia: 

```julia
EntanglementRobustness(phi_state)

>after 0 iterations: 0.33333341226620505
0.33333341226620505
```

The `MultiState` model also contains the function `RandomMultiState(dims::Array{Int,1})` for the creation of random (Hilbert-Schmidt distributed) states. `RandomMultiState([3, 4])`, for instance,  generates a random density matrix acting on the Hilbert space $\mathbb{C}^{3}\otimes \mathbb{C}^{4}$

In Julia:

```julia
rho = RandomMultiState([3, 4]);
EntanglementRobustness(rho)
>after 0 iterations: 0.5733033557264642
0.5733033557264642
```




## Advanced Usage 

The wrapper function `EntanglementRobustness` allows for the specification of several optional parameters. The most important are
+ `method::String` : can be "BP" or "PPT" (default `"BP"`)
+ `robustness_type::String` : can be "fixed", "absolute" or "general" (default `"fixed"`)
+ `noise_state::MultiState` : state which is used for mixing (default is not `nothing` which corresponds to white noise)
+ `convergence_recognition::Bool` : algorithm is stopped and value is returned when convergence is achieved (default `true`)
+ `convergence_accuracy::Float64` :  difference threshold for the determination of convergence (default `10^-4`)
+ `solver`: optimizer used for semidefinite programming (default `Mosek.Optimizer`)
+ `nvert::Int64` : number of vertices for the polytope approximation, when a random polytope should be used (default `200`)
+ `qpolytope::Vector{Matrix{ComplexF64}}` : for using a special polytope (default `nothing`)
+ `niter::Int64` : maximal number of iterations (default `5`)

Example for usage:

Let us generate a random $6\times6$-dimensional mixed quantum state.

```julia
rho = RandomMultiState([6, 6])

>MultiState(ComplexF64[0.02839554007403675 + 0.0im -0.0006940850554461725 + 0.00025000213362680605im … -0.005562177450860705 + 0.004619487053373366im 0.005163300878255143 + 0.0005893546427879063im; -0.0006940850554461725 - 0.00025000213362680605im 0.03082528883821383 + 0.0im … -0.0023268222283501 + 0.00839683587797709im -0.0039523366601357315 + 0.002421812552710086im; … ; -0.005562177450860705 - 0.004619487053373366im -0.0023268222283501 - 0.00839683587797709im … 0.035982840013786634 + 0.0im 0.004653518174286961 - 0.0026998505565306243im; 0.005163300878255143 - 0.0005893546427879063im -0.0039523366601357315 - 0.002421812552710086im … 0.004653518174286961 + 0.0026998505565306243im 0.027113898632547172 + 0.0im], [6, 6])
```
We want to specify 400 for the the number of vertices and 10 iterations and the algorithm should not stop after convergence.

```julia
EntanglementRobustness(rho, nvert=400, niter=10, convergence_recognition=false)

>after 0 iterations: 0.5466903200621531
after 1 iterations: 0.5614610524636662
after 2 iterations: 0.5614575353941491
after 3 iterations: 0.5614601765733851
after 4 iterations: 0.561455568191139
after 5 iterations: 0.5614604537465308
after 6 iterations: 0.5614586682079931
after 7 iterations: 0.5614617151388344
after 8 iterations: 0.5614521273931294
after 9 iterations: 0.5614539335065714
0.5614539335065714
```

# 2- Quantum Steerability

To be online soon!!!
