import Pkg
Pkg.activate("/home/matteo/Projects/PhyloBM/")
using PottsEvolver, FastaIO, TreeTools
include("../utils.jl")

#------PARAMETERS-------
n_samples = 13310
n_sweeps = 200
n_burnin = 100
n_sim = 5
output_root = "/home/matteo/Projects/PhyloBM/DBD/notebooks/mcmc_iid_reweight_randinit"
#-----------------------

println("Nmber of threads: ", Threads.nthreads())

J_reweight, h_reweight = load_potts_parameters("/home/matteo/Projects/PhyloBM/DBD/parameters/reweight/params.dat", aa_alphabet)
potts_reweight = PottsGraph(J_reweight, h_reweight, 1.0, aa_alphabet)
q, L = size(h_reweight)

Threads.@threads for i in 1:n_sim
    parameters_iid_mcmc= SamplingParameters(Teq = n_sweeps*L, burnin = n_burnin*L, step_meaning = :proposed)

    results_iid_mcmc_reweight = mcmc_sample(potts_reweight, n_samples, parameters_iid_mcmc, init=:random_aa)
    write_sequences(results_iid_mcmc_reweight.sequences.data, output_root * "_sim$i.fasta", aa_alphabet)
    println("Simulation $i/$n_sim completed.")
end


