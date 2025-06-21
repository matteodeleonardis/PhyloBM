#load Potts parameter from file
function load_potts_parameters(file_name, aa_alphabet)

    lines = readlines(file_name)
    L = parse(Int, split(lines[end])[2]) + 1
    q = length(aa_alphabet.char_to_index)
    J = zeros(q, q, L, L)
    h = zeros(q, L)
    for line in lines
        fields =split(line)
        if (length(fields) == 6) && (fields[1] == "J")
            _, i, j, ai, aj, val = fields
            J[aa_alphabet.char_to_index[ai[1]],  aa_alphabet.char_to_index[aj[1]], parse(Int, i)+1, parse(Int, j)+1] = parse(Float64, val)
            J[aa_alphabet.char_to_index[aj[1]], aa_alphabet.char_to_index[ai[1]], parse(Int, j)+1,  parse(Int, i)+1] = parse(Float64, val)
        elseif (length(fields) == 4) && (fields[1] == "h")
            _, i, ai, val = fields
            h[aa_alphabet.char_to_index[ai[1]], parse(Int, i)+1] = parse(Float64, val)
        else
            print(line * ": Unrecognized format for parameter.")
            break
        end
    end

    return J, h
end

#Sample the intra-pairwise distance of a set of sequences (Matrix)
function pairwise_dist(seqs::Matrix, n_samples)
    
    sample_idx = randperm(size(seqs,2))[1:n_samples]
    dist_pairwise = zeros(Int, n_samples*(n_samples-1)÷2)
    cnt =1
    for i in 1:n_samples
        for j in 1:i-1
            dist_pairwise[cnt] = sum(seqs[:,sample_idx[i]] .!= seqs[:,sample_idx[j]])
            cnt += 1
        end
    end

    return dist_pairwise
end


function freqs_from_matrix(z; q=20)
    q = max(q, maximum(z))
    L, M = size(z)
    fi = zeros(q, L)
    for m in 1:M
        for i in 1:L
            fi[z[i,m], i] += 1.0
        end
    end
    fi ./= M
    return fi
end

function fasta_to_matrix(fasta, aa_alphabet)
    seqs = map(x->x[2], readfasta(fasta))
    L = length(seqs[1])
    M = length(seqs)
    z = zeros(Int, L, M)
    for m in eachindex(seqs)
        @assert length(seqs[m]) == L
        for i in 1:L
            z[i,m] = aa_alphabet.char_to_index[seqs[m][i]]
        end
    end

    return z
end


function freqs_from_fasta(fasta, aa_alphabet; q=20)
    
    z = fasta_to_matrix(fasta, aa_alphabet)
    return freqs_from_matrix(z, q=q)
end


function create_logo_from_fasta(fasta, aa_alphabet; q=20)

    fi = freqs_from_fasta(fasta, aa_alphabet, q=q)
    q = size(fi, 1)
    logo=logo_from_matrix(fi,String(map(x->aa_alphabet.index_to_char[x], collect(1:q))))
    return logo
end

#Computes the correlation matrix of a set of sequences (Matrix)
function compute_corr_matrix(z, weights=[])
    q = max(20, maximum(z))
    L, M = size(z)

    if length(weights) == 0
        weights = ones(M)
    end
    Meff = sum(weights)

    fi = zeros(q, L)
    for m in 1:M
        for i in 1:L
            fi[z[i,m], i] += weights[m]
        end
    end
    fi ./= Meff
    fi = reshape(fi, q*L, 1)

    fij = zeros(q, L, q, L)
    for m in 1:M
        for i in 1:L
            for j in 1:L
                fij[z[i,m], i, z[j,m], j] += weights[m]
            end
        end
    end
    fij ./= Meff
    fij = reshape(fij, q*L, q*L)

    return fij .- fi*fi'
end

#Computes the cumulatice distrubution (cdf) of a vector of values
function compute_cdf(values)
    N = length(values)
    n_map = countmap(values)
    sorted_values = sort(unique(values)) #ascending
    cdf_values = zeros(length(sorted_values))
    n_cum = 0.0
    for i in eachindex(sorted_values)
        n_cum += n_map[sorted_values[i]]
        cdf_values[i] = (N - n_cum) / N
    end

    return sorted_values, cdf_values
end


function int_to_char(Z, aa_alphabet)
    sequences = Vector{String}(undef, size(Z,2))
    seq_char = Vector{Char}(undef, size(Z,1))
    for j in axes(Z, 2)
        for i in axes(Z, 1)
            seq_char[i] = aa_alphabet.index_to_char[Z[i,j]]
        end
        sequences[j] = String(seq_char)
    end
    return sequences
end


function int_to_char(seq_int, aa_alphabet)
    seq_char = Vector{Char}(undef, length(seq_int))
    for i in eachindex(seq_int)
        seq_char[i] = aa_alphabet.index_to_char[seq_int[i]]
    end
        
    return String(seq_char)
end


function simulate_iid_mc(n_samples, n_sweeps, potts_graph, init_seq)
    L = size(potts_graph.h, 2)
    parameters_mcmc = SamplingParameters(Teq=1, burnin=n_sweeps*L-1, step_meaning=:proposed)
    Z_mcmc = zeros(Int, L, n_samples)
    n_chunk_base = n_samples ÷ Threads.nthreads()
    n_add = n_samples % Threads.nthreads()
    n_chunk = [i <= n_add ? n_chunk_base+1 : n_chunk_base for i in 1:Threads.nthreads()]
    z_chunk = [zeros(Int, L, n_chunk[i]) for i in 1:Threads.nthreads()]
    Threads.@threads for i in 1:Threads.nthreads()
        results_mcmc = [mcmc_sample(potts_graph, 1, parameters_mcmc, init=init_seq) for _ in 1:n_chunk[i]]
        z_chunk[i] = hcat(map(x->x.sequences.data, results_mcmc)...)
    end
    return hcat(z_chunk...)
end


function simulate_iid_mc_and_write(n_samples, n_sweep, potts_graph, init_seq, output_file, aa_alphabet)
    z_mcmc = simulate_iid_mc(n_samples, n_sweep, potts_graph, init_seq)

    FastaWriter(output_file) do fw
        for i_seq in axes(z_mcmc, 2)
            write(fw, ">Seq_$i_seq\n")
            write(fw, int_to_char(z_mcmc[:,i_seq], aa_alphabet))
        end
    end
end


function simulate_tree(potts_graph, tree, init_seq)
    blm = PottsEvolver.BranchLengthMeaning(:sweep, :round)
    parameters_tree = SamplingParameters(Teq=0, burnin=0, step_meaning=:proposed, branchlength_meaning=blm)
    results_tree = mcmc_sample(potts_graph, tree, parameters_tree, init=init_seq)
    return results_tree.leaf_sequences.data
end


function simulate_tree_and_write(potts_graph, tree, init_seq, n_trees, output_file, aa_alphabet)
    z_mcmc = hcat([simulate_tree(potts_graph, tree, init_seq) for _ in 1:n_trees]...)

    FastaWriter(output_file) do fw
        for i_seq in axes(z_mcmc, 2)
            write(fw, ">Seq_$i_seq\n")
            write(fw, int_to_char(z_mcmc[:,i_seq], aa_alphabet))
        end
    end
end


function write_sequences(z, output_file, aa_alphabet)

    FastaWriter(output_file) do fw
        for i_seq in axes(z, 2)
            write(fw, ">Seq_$i_seq\n")
            write(fw, int_to_char(z[:,i_seq], aa_alphabet))
        end
    end
end
