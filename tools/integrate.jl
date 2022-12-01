using ArgParse2
using Artifacts 
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
using ProgressMeter

function main()
    parser = ArgumentParser(description = "Integrate clustering counts and assignment")
    add_argument!(parser, "--assembly-clusters", "-a", help = "Assembly clustering info",type = String)
    add_argument!(parser, "--deuplicated-clusters", "-d", help = "Assembly clustering info",type = String)
    add_argument!(parser, "--kraken", "-k", help = "Kraken output of pseudocontigs ",type = String)
    add_argument!(parser, "--kraken2", "-f", help = "Kraken output of fused contigs",type = String)
    add_argument!(parser, "--label", "-l", help = "Label for outout",type = String)
    args = parse_args(parser)
    taxdb=artifact"taxdump"
    
    #Parse data of fused clusters assembly
    fused_taxids = Array{String,1}()
    data_on_fused = DataFrame(Cluster_ID = Array{String,1}(),Nuber_of_contigs=Array{Int64,1}())
    for l in readlines(args.kraken2)
        parts = split(l,"\t")
        tid = parts[3]
        tid = split(split(split(tid,"taxid")[2])[1],")")[1]
        cid = split(parts[2],"@")[1]
        ncont = split(split(parts[2],"@")[2],":")[1]
        push!(data_on_fused,[cid,parse(Int64,ncont)])
        push!(fused_taxids,tid)
    end
    data_on_fused.Genus_based_on_fused_contigs = taxid_to_genus(fused_taxids,taxdb)
    data_on_fused.Species_based_on_fused_contigs = taxid_to_species(fused_taxids,taxdb)


    #Integrate deduplicated and clsutered for assembly clusters
    clusters_cf= args.assembly_clusters  
    clusters_df= args.deuplicated_clusters
    krakenout= args.kraken 
    stem=args.label

    cluster_4_workD,r2cD,c2rD = parse_mapping(clusters_df)
    cluster_4_workC,r2cC,c2rC = parse_mapping(clusters_cf)

    ct = 0
    clusterids = collect(keys(c2rC))
    cl_vs_abund = DataFrame(cluster=Array{String,1}(),cnt=Array{Float64,1}(),size=Array{Int64,1}())
    println(stderr, "Assigning counts")
    @showprogress for c in clusterids
        ct +=length(c2rC[c])
        c_reads=c2rC[c]
        fr = 0
        for cd in keys(c2rD)
            d_reads = c2rD[cd]
            totc = length(d_reads)
            overlapc = length(intersect(d_reads,c_reads))
            fr += overlapc/totc
        end
        push!(cl_vs_abund,[c,fr,length(c_reads)])
    end
    taxids=Array{String,1}()
    clid_spec=Dict{String,String}()
    clids=Array{String,1}()
    println(stderr,"Checking species and genus ssignments...")
    for l in readlines(krakenout)
        id = split(split(l)[2],"_")[1]
        push!(clids,id)
        taxid = split(split(split(l,"taxid")[2])[1],")")[1]
        push!(taxids,taxid)
    end
    species=taxid_to_species(taxids,taxdb)
    genus=taxid_to_genus(taxids,taxdb)
    println(stderr,"...done")
    clid_spec=Dict{String,String}(zip(clids,species))
    clid_genus=Dict{String,String}(zip(clids,genus))
    cl_vs_abund = sort(cl_vs_abund,[:size],rev=true)
    spassignments = Array{String,1}()
    geassignments = Array{String,1}()
    assigned = collect(keys(clid_spec))
    Cluster_ID = Array{String,1}()
    @showprogress for i in  1:nrow(cl_vs_abund)
        is = string(i)
        push!(Cluster_ID,is)
        if is in assigned
            push!(spassignments,clid_spec[is])
            push!(geassignments,clid_genus[is])
        else
            push!(spassignments,"NA")
            push!(geassignments,"NA")
        end
    end
    cl_vs_abund.Species_based_on_pseudocontigs = spassignments
    cl_vs_abund.Genus_based_on_pseudocontigs = geassignments
    cl_vs_abund.Fraction = cl_vs_abund.cnt ./ sum(cl_vs_abund.cnt)
    cl_vs_abund.Sample = fill(split(stem,"/")[end],length(spassignments))
    cl_vs_abund.Cluster_ID = Cluster_ID 
    out = join(cl_vs_abund, data_on_fused, on=:Cluster_ID, kind = :outer) 
    CSV.write(args.label*"_info_on_pseudocontig_and_contigs.csv",out)
    end

main()
