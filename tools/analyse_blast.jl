using ArgParse2
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using DataFramesMeta
using Statistics
using st16SseqJuliaTools
using XAM
using BioTools.BLAST
using CSV
using Artifacts 

function main()
    parser = ArgumentParser(prog = "analyse_alnself.jl",
                        description = "Detect chimeras",
                        epilog = "Removing potential chimeras and repair snps")
    add_argument!(parser, "--inxml", "-x", help = "input blast xml ",type = String)
    add_argument!(parser, "--ref", "-r", help = "reference fasta",type = String)
    add_argument!(parser, "--out", "-o", help = "Fasta file to write out cleaned up sequences",default="chimerasfree.faSTA",type = String)
    add_argument!(parser, "--min-seqid", "-s", help = "Minimum sequence identity in match",type = Float64, default=0.95)
    add_argument!(parser, "--min-covid", "-c", help = "Minimum sequence length fraction in  match",type = Float64, default=0.96)
    #add_argument!(parser, "--output", "-o", help = "File with names",type = String)
    
    maxNcount = 5 #max N counts in ref
    max_match_start_pos = 5

    args = parse_args(parser)
    ref = fasta_to_dict(args.ref)
    xmldata = open(args.inxml) do file
        read(file, String)
    end
    records = split_by_query(readblastXML2(xmldata))
    testid = "2_3209;size=1"
        
    
    #Limits for real match ... mayve should be options
    minf_len = args.min_covid
    minf_id = args.min_seqid # Identity of first 250 fragments of contigs
    pairs = Array{Tuple{String,String},1}()
    chosenseqs = OrderedDict{String,BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}()
    ct = 0
    info = DataFrame(contig = Array{String,1}(),
    coverage=Array{Float64,1}(),
    identity=Array{Float64,1}(),
    chimeric=Array{Bool,1}(),
    gaps=Array{Int64,1}(),
    max_indel_length=Array{Int64,1}(),
    good=Array{Bool,1}(),
    refNcounts=Array{Int64,1}(),
    match_from=Array{Int64,1}(),
    match_to=Array{Int64,1}(),
    totall_length=Array{Int64,1}(),
    match_taxid=Array{String,1}())
    for contig_id in keys(records)
        ct+=1
        matching_hits = get_coverage_maximising_hits(contig_id,records)
        r = matching_hits[1]
        contig_seq = ref[contig_id]
        gaps = r.gaps
        aln = r.alignment.aln
        count_amb = count(isambiguous, r.hit)
        if testid == contig_id
            println(r)
            println(r.hit)
            println(count_amb)
            for h in matching_hits
                println(h)
            end 
        end
        s = first(aln.anchors).seqpos
        e = last(aln.anchors).seqpos
        l=abs(e-s)
        #count count_deletions 
        indel_lengths = Array{Int64,1}()
        push!(indel_lengths,0)
        seqpos, refpos
        if length(aln.anchors) > 1
            for i in 2:length(aln.anchors)
                anc = aln.anchors[i]
                anb = aln.anchors[i-1]
                op = anc.op
                if isdeleteop(op) || isinsertop(op)
                    idl = maximum([abs(anc.refpos - anb.refpos),abs(anc.seqpos - anb.seqpos)])
                    push!(indel_lengths,idl)
                end 
            end 
        end
        contig_length = length(contig_seq)
        covf = round(l/contig_length ,digits=2)
        idf = round(r.identity/l, digits = 2)
        isgood = false
        ischimeric = (length(matching_hits) >1)
        push!(info,[contig_id,covf,idf,ischimeric,gaps,maximum(indel_lengths), isgood, count_amb, r.queryfrom, r.queryto,contig_length,r.hittaxid])
    end
    println("Matching genus to taxids...")
    species = taxid_to_species(info.match_taxid,artifact"taxdump")
    println("...done")
    info.genus = map(x -> split(x)[1],species)
    info.species_wog = map(x -> join(split(x),"_"),species)
    info.species = species
    info.cluster = map(x -> parse(Int64,split(x,"_")[1]),info.contig)
    info.matchlen = abs.(info.match_to .- info.match_from)
    sort!(info,[:cluster,:genus,:species])
    #Prefilter based on id and coverage for search of a typical species 
    infoF  = @where(info, :coverage .>= 0.90 , :identity .>= 0.96, :match_from .<=max_match_start_pos ,:refNcounts.<=maxNcount)
    sub = @select(infoF, :cluster, :genus, :contig,:matchlen)
    cl_s = Dict{Int64,String}() #expected genus per cluster
    gd = groupby(sub, :cluster)
    for cl in unique(sub.cluster)
        ginfo=get(gd, (cluster=cl,), nothing)
        sinfog= groupby(ginfo, :genus)
        s_vs_l = @combine(sinfog, Length = median(:matchlen))
        sort!(s_vs_l,[:Length],rev=true)
        s = s_vs_l.genus[1]
        cl_s[cl] = s
    end

    expected_genus = map(x -> cl_s[x],info.cluster)
    info.Typical_cluster_genus = (info.genus .== expected_genus) 
    
    #test = @where(info, map(x -> occursin("8_246",x),:contig))
    #println(test)
    #exit()

    #good  = @where(info, :coverage .> 0.96 , :identity .> 0.98, .!:chimeric).contig
    good  = @where(info, :gaps .<= 5,:max_indel_length .<=2 , :coverage .>= 0.90 , :identity .>= 0.96, :match_from .<=max_match_start_pos ,:refNcounts.<=maxNcount,:Typical_cluster_genus).contig
    info.good = map(x->(x in good),info.contig)
    CSV.write(args.out*"_blastanalysis.csv",info)
    #collect info about end positions fro write out 
    ends4writeout = Dict(zip(info.contig,info.match_to))
    species4writeout = Dict(zip(info.contig,info.species_wog))
    ct_good = 0
    for k in keys(ref)
        if k  in good
            ct_good += 1
            ko = k*"@"*species4writeout[k]
            chosenseqs[ko] = ref[k][1:ends4writeout[k]]
        end 
    end
    ct_all = length(keys(ref))
    
    println("$ct_all in total sequences, $ct_good sequences were  written to ",args.out)
    dict_to_fasta(chosenseqs,args.out)
end

main()
