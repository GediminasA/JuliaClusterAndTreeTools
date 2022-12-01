using DataFrames
using DataFramesMeta
using CSV

blastfile = ARGS[1]
dadatax = ARGS[2]
outstem = ARGS[3]

blast_fields = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen"
dada_fields = "Contig Kingdom Phylum Class Order Family Genus"
df = CSV.File(blastfile,header = Array{String,1}(split(blast_fields)))
df2 = CSV.File(dadatax,header = Array{String,1}(split(dada_fields)), skipto = 2)
df = DataFrame(df)
df2 = DataFrame(df2)
#ensure that contig names are strings
df.qseqid = string.(df.qseqid)
df2.Contig = string.(df2.Contig)
tt = getindex.(split.(df2.Contig," "),1)

df2.Contig = tt 
contig_genus_map = Dict{String,String}()
contigs = collect(df2.Contig)
genus = df2.Genus
for i in 1:length(contigs)
    contig_genus_map[df2.Contig[i]] = df2.Genus[i]
end
sizes = fill(0,size(df)[1])
species = fill("",size(df)[1])
read_name = fill("",size(df)[1])
genus = fill("",size(df)[1])
genus_species = fill("",size(df)[1])
for i in 1:size(df)[1]
    if occursin("=",df.qseqid[i]) && occursin(";",df.qseqid[i]) 
        sizes[i] = parse(Int64,split(split(df.qseqid[i],"=")[2],";")[1]) #swarm attaches ; to the end of name
    else
        sizes[i] = 1
    end 
    read_name[i] = split(df.qseqid[i],";")[1]
    
    species[i] = String(split(df.sseqid[i],":")[1])
    
    genus[i] = contig_genus_map[df.qseqid[i]]
    genus_species[i] = replace(genus[i]*"||"*species[i],"/"=>"-")
end
df.sizes = sizes
df.species = species
df.read_name = read_name
df.genus = genus
df.genus_species = genus_species


dfout = DataFrame(Genus_Species = Array{String,1}(),
    Maximum_length = Array{Int64,1}(),
    Number_of_sequences = Array{Int64,1}(),
    Maximum_abundance = Array{Int64,1}(), #parses size argument in name, not the real number od sequences
    Total_abundance = Array{Int64,1}(), #parses size argument in name, not the real number od sequences
    Samples = Array{String,1}())
for gs in unique(df.genus_species)
    dfs = @linq df |>
               where(:genus_species .== gs) |>
               orderby(-:sizes)
    dfs_ids = @linq dfs |>
        select(:read_name)
    out = outstem*"_ids_of_"*replace(gs,"/"=>"|")
    CSV.write(out,dfs_ids,header=[""])
    maxl = maximum(dfs.qlen)
    maxnb = maximum(dfs.sizes)
    totalnb = sum(dfs.sizes)
    smp = split(outstem,'/')[end]
    push!(dfout,[gs,maxl,size(dfs)[1],maxnb, totalnb, smp])
end
out = outstem*"_blast_summary.tsv"
CSV.write(out,dfout,delim='\t')
#group by genus
dfout = DataFrame(Genus = Array{String,1}(),
    Maximum_length = Array{Int64,1}(),
    Number_of_sequences = Array{Int64,1}(),
    Maximum_abundance = Array{Int64,1}(), #parses size argument in name, not the real number od sequences
    Total_abundance = Array{Int64,1}(), #parses size argument in name, not the real number od sequences
    Samples = Array{String,1}())
for g in unique(df.genus)
    dfs = @linq df |>
               where(:genus .== g) |>
               orderby(-:sizes)
    dfs_ids = @linq dfs |>
        select(:read_name)
    out = outstem*"_GENUS_ids_of_"*replace(g,"/"=>"|")
    CSV.write(out,dfs_ids,header=[""])
    maxl = maximum(dfs.qlen)
    maxnb = maximum(dfs.sizes)
    totalnb = sum(dfs.sizes)
    smp = split(outstem,'/')[end]
    push!(dfout,[g, maxl, size(dfs)[1], maxnb, totalnb, smp])
end
out = outstem*"_blast_summary_genus.tsv"
CSV.write(out,dfout,delim='\t')

