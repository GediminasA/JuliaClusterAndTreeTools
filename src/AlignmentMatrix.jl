
mutable struct Alignment
    names::Array{String,1} # sequence names
    M::Array{DNA,2}
    repr_n:: Array{Int32,1} #number of epresentated species/sequiences - used for
    name_sequence_map:: Dict{String,Int32}

    function Alignment(namesini::Array{String,1}, M::Array{DNA,2})
        repr_n = fill(Int32(1), size(M,1))
        names = fill("", size(M,1))
        Threads.@threads for i in eachindex(namesini)
            parts = split(namesini[i], ";")
            names[i] = parts[1]
            if length(parts) > 1
                repr_n[i] = parse(Int32, split(parts[2], "=")[2])
            else
                repr_n[i] = 1
            end
        end
        return Alignment(names, M, repr_n)
    end

    function Alignment(names::Array{String,1}, M::Array{DNA,2}, repr_n::Array{Int32,1})
        name_sequence_map = Dict{String,Int32}()
        for i in eachindex(names)
            name_sequence_map[names[i]] = i
        end
        return new(names, M, repr_n, name_sequence_map)
    end

    #primary=true  when the initial alignment is beeing read in (as a reference)
    function  Alignment(file_name::String, allowed_symbols::Set{Char}=Set{Char}(['A','T','G','C','N','-']))::Alignment
        sequences = Array{String,1}()
        names = Array{String,1}()
        num = 0
        lengths = Array{Int32,1}()
        skiped = 0
        println(stderr, "Reading aligment $file_name\n")
        FastaReader(file_name) do fr
            for (n, seq) in fr
                seq_ar::Vector{Char} = collect(uppercase(seq))
                proper = true
                ct = 0
                for s in seq_ar
                    ct += 1
                    proper = false
                    for sa in allowed_symbols
                        if s == sa
                            proper = true
                            break
                        end
                    end
                    if !proper
                        seq_ar[ct] = 'N'
                    end
                end
                seq_str = join(seq_ar)
                num += 1
                if mod(num, 1000) == 0
                    print(stderr,"Parsed $num sequences \r")
                end
                push!(names, n)
                push!(sequences, seq_str)
                push!(lengths, length(seq_str))
            end
        end
        len = maximum(lengths)
        #create alignment array
        m=fill(DNA_Gap, num, len)
        for i in 1:num
            m[i,1:length(sequences[i])] = collect(BioSequences.LongDNASeq(sequences[i])) 
            sequences[i]=Base.undef_ref_str
        end
        sequences = nothing
        return Alignment(names, m)
    end
end

#sub alignment  - column slice
function sub_alignment(aln::Alignment, interval::UnitRange{Int64})
    names = aln.names
    repr_n = aln.repr_n
    min = aln.M
    mout = Array{Char}(undef, size(min)[1], length(interval))
    for i in 1:size(min)[1]
        mout[i,:] = min[i, interval]
    end
    return(Alignment(names, mout, repr_n))
end

#sub alignment  - column slice
function sub_alignment!(aln::Alignment, interval::UnitRange{Int64})
    aln.M=aln.M[:,interval]
end

#sub alignment  - subset on particular names
function sub_alignment(aln::Alignment, names::Array{String,1})
    new_names = Array{String,1}() # sequence names
    new_M = Array{Char}(undef, length(names), size(aln.M)[2])
    new_repr_n = Array{Int32,1}() #number of epresentated species/sequiences - used for
    ct = 0
    for n in names 
        ct += 1
        ind = aln.name_sequence_map[n]
        push!(new_names, aln.names[ind])
        new_M[ct, :] = aln.M[ind, :]
        push!(new_repr_n, aln.repr_n[ind])
    end 
    return(Alignment(new_names, new_M, new_repr_n))
end

#sub alignment  - subset on particular names
function sub_alignment!(aln::Alignment, names::Array{String,1})
    new_names = Array{String,1}() # sequence names
    new_repr_n = Array{Int32,1}() 
    ct = 0
    for n in names 
        ct += 1
        ind = aln.name_sequence_map[n]
        push!(new_names, aln.names[ind])
        push!(new_repr_n, aln.repr_n[ind])
    end 
    keep_indexes = [ aln.name_sequence_map[n] for n in names]
    aln.M = aln.M[keep_indexes,:]
    aln.names  = new_names
    aln.repr_n = new_repr_n
    aln.name_sequence_map = Dict{String,Int32}()
    for i in eachindex(aln.names)
        aln.name_sequence_map[aln.names[i]] = i
    end
 end

#sub alignment  - leaves on rows based on iondex
function sub_alignment_rowwise!(aln::Alignment, rows::Array{Int64,1})
    aln.names = aln.names[rows]
    aln.repr_n = aln.repr_n[rows]
    aln.M = aln.M[rows,:]
    aln.name_sequence_map = Dict(zip(aln.names,1:length(aln.names)))
    aln.name_sequence_map = Dict{String,Int32}()
    for i in eachindex(aln.names)
        aln.name_sequence_map[aln.names[i]] = i
    end
end

#sub alignment  - column slice
function sub_alignment_columnwise!(aln::Alignment, columns::Array{Int64,1})
    aln.M=aln.M[:,columns]
end


function write_to_fasta(aln::Alignment, out_f_name::String;write_sizes = false)
    f = GZip.open(out_f_name,"w")
    for i in 1:size(aln.M, 1)
        n = aln.names[i]
        rn = aln.repr_n[i]
        if write_sizes
            outn = "$n;size=$rn"
        else
            outn = "$n"
        end
        s = join(aln.M[i,:])
        write(f, ">$outn\n$s\n")
    end
    close(f)
end

function write_to_fasta(aln::Alignment, interval::UnitRange{Int64}, out_f_name::String;remove_gaps=false)
    f = open(out_f_name, "w")
    for i in 1:size(aln.M, 1)
        n = aln.names[i]
        rn = aln.repr_n[i]
        outn = "$n;size=$rn"
        s = String(aln.M[i,interval])
        if remove_gaps
            s=replace(s, "-"=>"N")
        end
        print(f, ">$outn\n$s\n")
    end
    close(f)
end

"""
    Collect data on Gaps and Ns abundance collumnwise
"""
function get_statistics_on_Gaps_Ns_collumnwise(aln::Alignment)
    m = aln.M
    set_dna = [DNA_N,DNA_Gap]
    out = Array{Array{Int64},1}()
    ctparse = 0
    for clm in eachcol(m)
        cm = countmap(clm)
        N = length(clm)
        cmvals =  map( (x) -> x in keys(cm) ? cm[x] : 0, set_dna)
        push!(cmvals,N)
        push!(out,cmvals)
        ctparse += 1
        if mod(ctparse,1000) == 0
            print("Counted $ctparse columns\r")
        end
    end
    println("Doing transposition...")
    outm = transpose(hcat(out...))
    println("Done. Converting to a data frame...")
    outd = DataFrame(outm,[:N,:Gap,:Total])
    outd.Position = 1:length(outm[:,1])
    println("Done")
    return(outd)
end

"""
    Collect data on Gaps and Ns abundance rowwise
"""
function get_statistics_on_Gaps_Ns_rowwise(aln::Alignment)
    m = aln.M
    set_dna = [DNA_N,DNA_Gap]
    out = Array{Array{Int64},1}()
    ctparse = 0
    for clm in eachrow(m)
        cm = countmap(clm)
        N = length(clm)
        cmvals =  map( (x) -> x in keys(cm) ? cm[x] : 0, set_dna)
        push!(cmvals,N)
        push!(out,cmvals)
        ctparse += 1
        if mod(ctparse,1000) == 0
            print("Counted $ctparse rows\r")
        end
    end
    println("Doing transposition...")
    outm = transpose(hcat(out...))
    println("Done. Converting to a data frame...")
    outd = DataFrame(outm,[:N,:Gap,:Total])
    outd.ID = aln.names
    outd.NR = 1:length(aln.names)
    return(outd)
end
"""
 Cluster a reagion in an alignemnt and replace custered fragments with a cluster representative
"""
function vsearh_cluster!(aln::Alignment, interval::UnitRange{Int64}, id::Float64)
    #create temorary file in dev/shm
    bindir = Conda.bin_dir(:vsearch)
    tmpdir = mktempdir("/dev/shm"; prefix="vsearch_", cleanup=false)
    in_fasta = "$tmpdir/input.fasta"
    out_fasta = "$tmpdir/centroids.fasta"
    clusters = "$tmpdir/centroids.uc"
    write_to_fasta(aln, interval, in_fasta, remove_gaps=true)
    cluster = run(`$bindir/vsearch --quiet --sizein --cluster_size $in_fasta --uc $clusters --id $id`)
    clustermap = get_clustering(clusters)
    for i in eachindex(aln.names)
        aln.M[i,interval] = aln.M[aln.name_sequence_map[clustermap[aln.names[i]]], interval]
    end
end

function get_max_cl_length(aln::Alignment)
    #create temorary file in dev/shm
    bindir = Conda.bin_dir(:vsearch)
    tmpdir = mktempdir("/dev/shm"; prefix="vsearch_", cleanup=true)
    lenin = size(aln.M)[2]
    target = size(aln.M)[1]
    nmb = 0
    ln = 0
    for i in 10:lenin
        clusters = "$tmpdir/$i.centroids.uc"
        in_fasta = "$tmpdir/$i.input.fasta"
        write_to_fasta(aln, 1:i, in_fasta, remove_gaps=true)
        cluster = run(`$bindir/vsearch --quiet --sizein --cluster_size $in_fasta --uc $clusters --id 0.96`)
        clustermap = get_clustering(clusters)
        ln = i
        nmb = length(unique(values(clustermap)))
        if nmb == target
            break
        end
    end
    println("$nmb serquences at $ln length")
end

function vsearh_cluster!(aln::Alignment, window::Int64, step::Int64, id::Float64)
    lenaln = size(aln.M)[2]
    final = lenaln - window + 1
    starts = collect(1:step:final)
    nb = length(starts)
    ct = 0
    for startsmooth in starts
        ct += 1
        startl = startsmooth
        endl = startsmooth + window -1
        if ct == nb
            endl = lenaln
        end
        vsearh_cluster!(aln, startl:endl, id)
    end
end

"Patses uclust file and collects mapping"
function get_clustering(ucinfile::String)
    out = Dict{String,String}()
    curr = " "
    for l in readlines(ucinfile)
        if l[1] == 'H'
            parts = split(l)
            mem, rep = split(parts[9],";")[1], split(parts[10],";")[1]
            out[mem]  =rep
        end
        if l[1] == 'C'
            parts = split(l)
            mem = split(parts[9], ";")[1]
            out[mem] = mem
        end
    end
    return(out)
end

function char_array_entropy(s::String;debug=false)
    return char_array_entropy(collect(s), debug=debug)
end

function char_array_entropy(s::Array{Char,1}, c::Array{Int,1};debug=false)
    if length(s) != length(c)
        error("Number of symbols and their counts do not match")
    end
    #not perfect way...just generate an array
    chars = Array{Char,1}()
    for i in 1:length(s)
        part=fill(s[i], c[i])
        for cc in part 
            push!(chars, cc)
        end 
    end
    return(char_array_entropy(chars, debug=debug))
end

function char_array_entropy(s::Array{Char,1};debug=false)
    c = counter(s)
    n = length(s)
    e=0
    consensus = ' '
    for k in keys(c)
    p = c[k]/n
    e += -1*p*log2(p)
    end
    frac_gap = c['-']/n
    if debug
        println((e, " ", frac_gap))
    end
    consensus = sort(collect(c), by=x->-x[2])
    vals = map((x) -> x[2][1], consensus)
    pvals = maximum(vals)/sum(vals)
    return(e, frac_gap, pvals, consensus[1][1])
end

function aln_conservation(aln::Alignment; ignore_amounts = false)
    entropies = Array{Float64,1}()
    conservation = Array{Float64,1}()
    sz = size(aln.M)
    for i in 1:sz[2]
        chars = collect(join(aln.M[:,i]))
        cnts = aln.repr_n
        if ignore_amounts 
            entr = char_array_entropy(chars)
        else 
            entr = char_array_entropy(chars, Array{Int64,1}(cnts))
        end 
        push!(entropies,entr[1])
        push!(conservation,entr[3])
    end
    return(entropies, conservation)
end

"Identfies sequences which are recombinations of others and has no unique snps"
function eval_recomb(aln::Alignment)
    out_dataframe=DataFrame(Read_name = Array{String,1}(),
                Unique_snps = Array{Int64,1}(),
                Maximum_common_fraction = Array{Float64,1}(),
                Maximum_common_fraction_names = Array{String,1}())
                
    aln_c = aln_conservation(aln)
    m_sub_ini = aln.M[:, aln_c .< 2]
    m_sub = copy(m_sub_ini)
    n_seqs = size(aln.M)[1]
    seqs_length = size(aln.M)[2]
    fractions = Dict{String,Array{Float64,1}}()
    fractions_names = Dict{String,Array{String,1}}()
    for i in 1:n_seqs
        fractions[aln.names[i]] = Array{Float64,1}()
        fractions_names[aln.names[i]] = Array{String,1}()
        for j in 1:n_seqs
            if i != j && (aln.repr_n[i] <  aln.repr_n[j]) 
                #collect fraction of overlaping snps independant for each reference
                replaced,frac = remove_matching(m_sub_ini[i,:], m_sub_ini[j,:])
                push!(fractions[aln.names[i]], frac) 
                push!(fractions_names[aln.names[i]], aln.names[j])
                #collect comparisions
                replaced_acc2,frac2 = remove_matching(m_sub[i,:], m_sub[j,:])
                m_sub[i,:] = replaced_acc2  
            end 
        end 
    end
    out_nonunic = Array{Int64,1}()
    out_max_frac = Array{Float64,1}()
    out_max_frac_names = Array{Float64,1}()
    for i in 1:n_seqs 
        nb_non_unic = 0
        if '-' in m_sub[i,:] 
            nb_non_unic = counter(m_sub[i,:])['-']
        end 
        max_frac_common = 0
        max_frac_common_name = ""
        if length(fractions[aln.names[i]]) > 0
            max_frac_common = maximum(fractions[aln.names[i]])
            max_frac_common_name = join( fractions_names[aln.names[i]][ fractions[aln.names[i]] .== max_frac_common ], " ")
        end
        nb_unic = length(m_sub[i,:]) - nb_non_unic 
        push!(out_dataframe, [aln.names[i], nb_unic, max_frac_common, max_frac_common_name])
    end
    return(out_dataframe)
end

"changes a symbol in act array to '' if in the ref array is the same symbol"
function remove_matching(act::Array{Char,1}, ref::Array{Char,1})
    out = copy(act)
    ct = 0
    for i in 1:length(act)
        if act[i] == ref[i]
            out[i] = '-'
            ct += 1
        end 
    end
    
    return(out, ct/length(act))
end 
