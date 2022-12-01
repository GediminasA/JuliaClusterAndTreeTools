using st16SseqJuliaTools
using ArgParse2
using DataFrames 
using Query 
using CSV 
function julia_main()::Cint
    parser = ArgumentParser( description = "Extract clusters of a proper size and without close chimeras")
    #add_argument!(parser, "-i", "--input-fasta", type = String, help = "Input file")
    add_argument!(parser, "-a", "--cluster-lower-level", type = String, help = "Clustering to group reads for assembly")
    add_argument!(parser, "-b", "--cluster-higher-level", type = String, help = "Clustering to group reads to check 4 chimeras")
    add_argument!(parser, "-t", "--tax-assign", type = String, help = "Taxonomy assignment files - down to genus level. Should match the lower level clustering fasta file ")
    add_argument!(parser, "-m", "--min-cluster-size", type = Float64, help = "Minimum cluster size")
    add_argument!(parser, "-f", "--fasta", type = String, help = "Fasta file matching lower lever of clustering ")
    add_argument!(parser, "-d", "--dir", type = String, help = "Directory for output", default = "grouping_info")
    args = parse_args(parser)
    if !isdir(args.dir)
        mkdir(args.dir)
    end 
    aln_all = Alignment(args.fasta) 
    centroid_genus = parse_dada_results(args.tax_assign)
    cl1  = read_jc_clusters(args.cluster_lower_level)
    #check which cluster has proper size
    proper_sizes = mark_proper_size_clusters(cl1,args.min_cluster_size)
    #get higher level cluster - it should group such close species as E. Coli & Salmonela and their chimeric reads
    cl2   = read_jc_clusters(args.cluster_higher_level)
    recomb_data = DataFrame()
    groupct = 0
    for cl in keys(cl2) 
        ids_for_common_analyses = Array{String,1}()
        for rep in cl2[cl]
                if proper_sizes[rep]
                    push!(ids_for_common_analyses,rep)
                end 
        end 
        if length(ids_for_common_analyses) > 2
            #a case where chimeric cluster can be formed
            groupct += 1
            #println("------------------")
            #for id in ids_for_common_analyses 
            #    println(id," ",length(cl1[id])," ",centroid_genus[id]," ", proper_sizes[id])
            #end 
            sub_aln = sub_alignment(aln_all,ids_for_common_analyses)
            recombs = eval_recomb(sub_aln)
            recombs.Group = fill(groupct,nrow(recombs))
            recomb_data = vcat(recomb_data,recombs)
        end 
    end
    chimof = args.dir*   "/chimeric.csv"
    println(stderr,"Chimeric reads detection info: ", chimof)
    CSV.write(chimof, recomb_data) 
    if nrow(recomb_data) > 0
        chimeric_ids = @from i in recomb_data begin
                                @where i.Unique_snps  < 2
                                @select { i.Read_name}
                                @collect DataFrame 
                        end
    else 
        chimeric_ids = DataFrame(Read_name = [])
    end 
    #filter out detection


    
    #aln = sanitise.Alignment(args.input_fasta)
    #println(aln_conservation(aln,ignore_amounts=true))
    #construct main output data frame
    finald = DataFrame()
    read_names = collect(keys(cl1))
    finald.Centroids = read_names 
    finald.Cluster_name = map( x -> centroid_genus[x], read_names) 
    finald.Genus = map( x -> split(centroid_genus[x],"_")[1], read_names) 
    finald.Size = map(x -> length(x),values(cl1))
    finald.Proper_size  = collect(values(proper_sizes)) 
    finald.Chimeric = map(x -> (x in chimeric_ids.Read_name),read_names)
    allof = args.dir* "/all_clusters.csv"
    println(stderr,"Data on all cluster: ", allof)
    CSV.write(allof, finald)
    chosen_clusters = @from i in finald begin
                            @where i.Proper_size  #&& !i.Chimeric        
                            @select { i.Cluster_name, i.Genus,i.Centroids , i.Size, i.Chimeric}
                            @collect DataFrame 
                        end 
    
    insertcols!(chosen_clusters, 1,:ID => collect(1:nrow(chosen_clusters)))
    chosenof = args.dir*   "/chosen_clusters.csv"
    println(stderr,"Data on chosen cluster: ", chosenof)
    CSV.write(chosenof, chosen_clusters)
    #output contigs
    ct = 0
    for r in eachrow(chosen_clusters)
        ct = r.ID  
        f=open(args.dir*"/$ct","w")
        for id in cl1[r.Centroids]
            println(f,id)
        end 
        close(f)
        f=open(args.dir*"/$ct."*r.Cluster_name,"w")
        println(f,length(cl1[r.Centroids]))
        close(f)
    end


    return 0
end

julia_main()
