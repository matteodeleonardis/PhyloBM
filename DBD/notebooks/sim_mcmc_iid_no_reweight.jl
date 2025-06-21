import Pkg
Pkg.activate("/home/matteo/Projects/PhyloBM/")
using PottsEvolver, FastaIO, TreeTools
include("../utils.jl")

#------PARAMETERS-------
n_samples = 13310
n_sweeps = 200
n_burnin = 100
n_sim = 5
output_root = "/home/matteo/Projects/PhyloBM/DBD/notebooks/mcmc_iid_no_reweight"
#-----------------------

println("Nmber of threads: ", Threads.nthreads())

J_no_reweight, h_no_reweight = load_potts_parameters("/home/matteo/Projects/PhyloBM/DBD/parameters/no_reweight/params_no_reweight.dat", aa_alphabet)
potts_no_reweight = PottsGraph(J_no_reweight, h_no_reweight, 1.0, aa_alphabet)
q, L = size(h_no_reweight)

wt_fasta = readfasta("/home/matteo/Projects/PhyloBM/DBD/notebooks/DBD_WT.fasta")
wt_int = map(x->aa_alphabet.char_to_index[x], collect(wt_fasta[1][2]))

Threads.@threads for i in 1:n_sim
    parameters_iid_mcmc= SamplingParameters(Teq = n_sweeps*L, burnin = n_burnin*L, step_meaning = :proposed)

    results_iid_mcmc_no_reweight = mcmc_sample(potts_no_reweight, n_samples, parameters_iid_mcmc, init=wt_int)
    write_sequences(results_iid_mcmc_no_reweight.sequences.data, output_root * "_sim$i.fasta", aa_alphabet)
    println("Simulation $i/$n_sim completed.")
end


