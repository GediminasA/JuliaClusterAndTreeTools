using DataStructures
using CSV 


"Parase genus assignment from a dada results file"
function parse_dada_results(fd::String)
    centroid_genus_map = Dict{String,String}()
    ct = 0
    nact = 0
    totalnr = 0
    genuses = Set{String}()
    genuses_counts = Dict{String,Int64}()
    for l in readlines(fd)
        ct += 1
        if ct > 1
            parts = split(replace(l, "\""=>""), "\t")
            rn = split(parts[1], ";")[1]
            size = parse(Int64, split(split(parts[1], ";")[2], "=")[2])
            totalnr += size
            gn = replace(parts[7], " "=>".")
            gn = replace(gn, "/"=>"-")
            if gn == "NA"
                nact += 1
                gn = "$gn"*"_"*"$nact"
            else 
                if gn in genuses 
                    genuses_counts[gn] += 1
                else 
                    push!(genuses,gn)
                    genuses_counts[gn] = 1
                end 
                gnnr = genuses_counts[gn]
                gn = "$gn"*"_"*"$gnnr"
            end
            centroid_genus_map[rn] = gn
        end
    end
    return(centroid_genus_map)
end 

function fasta_to_dict(sfasta::String;remove_abundance=false)::OrderedDict{String,LongSequence{BioSequences.DNAAlphabet{4}}}
    out = OrderedDict{String,LongSequence{BioSequences.DNAAlphabet{4}}}()        
    record = FASTA.Record()
    reader = FASTA.Reader(open(sfasta, "r"))
    while !eof(reader)
        read!(reader, record)
        id = identifier(record)
        if remove_abundance
            id = split(id, ";")[1]
        end 
        seq = sequence(record)
        out[id]=seq
    end
    close(reader)
    return(out)
end 

"Read salmon output"
function read_salmon(inf::String;remove_abundance=true) #by defaulyt remove ;size= fragment
    din = DataFrame(CSV.File(inf))
    if remove_abundance
        din.Name = map(x->split(x, ";")[1], din.Name)
    end 
    out = Dict(zip(din.Name, din.TPM))
    return(out)
end

function dict_to_fasta(inp::OrderedDict{String,BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}, sfasta::String)
    writefasta(sfasta, inp) 
end 

#test aligned part percentage identity and length fraction compared to the s2 sequence
function check_pair(s1::BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}, s2::BioSequences.LongSequence{BioSequences.DNAAlphabet{4}})
    scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1)
    res = pairalign(LocalAlignment(), s1, s2, scoremodel)
    aln = alignment(res)
    pairs = collect(aln)
    l = count_aligned(aln)
    match = count_matches(aln)
    mismatch = BioAlignments.count_mismatches(aln)
    identityf = round(match/l, digits=3)
    coveragef = Base.maximum([(match + mismatch)/length(s1), ((match + mismatch)/length(s2))])
    coveragef = round(coveragef, digits=3)
    return(identityf, coveragef)
end

#read mapping of gjc file - "tracking of clusterings file"
function parse_mapping(mapfile::String)
    pairs = Array{Tuple{String,String},1}()
    r2c = Dict{String,String}() #read to cluster representative
    c2r = Dict{String,Array{String,1}}() #cluster rep to reads map
    for l in readlines(mapfile)
        parts = split(l)
        push!(pairs, (parts[1], parts[2]))
        r2c[String(parts[1])] = String(parts[2])
    end
    sort!(pairs, lt=(x,y)->isless(x[2], y[2]))
    cur_r2 = ""
    len = length(pairs)
    ct = 0
    cl_ct = 0
    r1s = Array{String,1}()
    cluster_4_work = Array{Tuple{Array{String,1},String},1}()
    @showprogress 1 "Parsing mapping..."  for pair in pairs
        ct += 1
        if ct == 1
            cur_r2 = pair[2]
        end
        if pair[2] != cur_r2 && ct !=1
            cl_ct += 1
            push!(cluster_4_work, (r1s, cur_r2))
            r1s = Array{String,1}()
            cur_r2 = pair[2]
        end
        c2r[cur_r2] = r1s
        push!(r1s, pair[1])
    end
    push!(cluster_4_work, (r1s,cur_r2))
    c2r[cur_r2] = r1s
    return(cluster_4_work, r2c, c2r)
end
