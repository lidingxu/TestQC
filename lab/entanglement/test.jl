#This work for all subdiectories of quantum-correlations
cpath=pwd();
cpath=split(cpath,"quantum-correlations")[1];

push!(LOAD_PATH,cpath*"quantum-correlations/lib/kvant/julia/");
push!(LOAD_PATH,cpath*"quantum-correlations/src/entanglement/");
push!(LOAD_PATH,cpath*"quantum-correlations/src/steerability/");
using Ket
using MultiStates
using Entanglement
print(1- EntanglementRobustness( MultiState(Matrix(state_w(Complex{Float64}, 4)),fill(2,4)), method = "BP", silent = false))

