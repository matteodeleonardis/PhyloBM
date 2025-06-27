import Pkg
Pkg.activate("/home/matteo/Projects/PhyloBM/")
using PottsEvolver, FastaIO, TreeTools
include("../utils.jl")

#------PARAMETERS-------
n_sim = 15
offset_sim = 5
output_root = "/home/matteo/Projects/PhyloBM/DBD/notebooks/mcmc_tree_accepted"
#-----------------------



println("Nmber of threads: ", Threads.nthreads())

J_reweight, h_reweight = load_potts_parameters("/home/matteo/Projects/PhyloBM/DBD/parameters/reweight/params.dat", aa_alphabet)
potts_reweight = PottsGraph(J_reweight, h_reweight, 1.0, aa_alphabet)
q, L = size(h_reweight)

wt_fasta = readfasta("/home/matteo/Projects/PhyloBM/DBD/notebooks/DBD_WT.fasta")
wt_int = map(x->aa_alphabet.char_to_index[x], collect(wt_fasta[1][2]))

tree = read_tree("/home/matteo/Projects/PhyloBM/DBD/FastTree_out/DBD_tree.nwk")
blm = BranchLengthMeaning(:sweep, :round)
Threads.@threads for i in 1:n_sim
    parameters_tree = SamplingParameters(Teq=0, burnin=0, step_meaning=:accepted, branchlength_meaning=blm)
    results_tree = mcmc_sample(potts_reweight, tree, parameters_tree, init=wt_int)
    write_sequences(results_tree.leaf_sequences.data, output_root * "_sim$(offset_sim + i).fasta", aa_alphabet) 
    println("Simulation $i/$n_sim completed.")
end




