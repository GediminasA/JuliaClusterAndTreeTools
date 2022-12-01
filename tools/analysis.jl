using BioSequences
using FASTX
using Test
using DataFrames
using CSV 
#tes Alt parsing and repr update saving writing to fasta
using st16SseqJuliaTools

primers = "/home/gediminas/working/primers/primers.fasta"
db_file = "/home/gediminas/Desktop/tm/GRCh37_latest_rna.fna"



#collect primer  st16SseqJuliaTools
primer_inf = Array{Tuple{String,String},1}()
reader = FASTA.Reader(open(primers, "r"))
record = FASTA.Record()

ctp = 0
while !eof(reader)
    global ctp += 1
    read!(reader, record)
    id = identifier(record)
    seq = sequence(record)
    push!(primer_inf,(id,seq))
    if ctp == 100
        break 
    end
end
close(reader)

data = DataFrame(
 ID=Array{String,1}(),   
 SEQ=Array{String,1}(),   
 K=Array{Int64,1}(),   
 K2=Array{Int64,1}(),   
 N=Array{Int64,1}(),   
 NMBBIND=Array{Int64,1}(), 
 TAIL=Array{Float64,1}(), 
 E=Array{Float64,1}(), 
 KMER=Array{String,1}(),
 KMER_COUNT=Array{Int64,1}(),
 KMER_MSIMAT=Array{Int64,1}(),
 TM_DNA_DNA = Array{Float64,1}(),
 TM_DNA_RNA_HS = Array{Float64,1}(),
 TM_DNA_RNA_LS = Array{Float64,1}(),
 TM_DNA_DNA_KMER = Array{Float64,1}(),
)

ctn = 0 
for (k,k2,n) in [(15,3,1),(15,2,1),(15,0,1),(15,3,2),(15,2,2),(15,0,2),(15,3,3),(15,2,3),(15,0,3),(6,3,1),(6,2,1),(6,0,1),(6,3,2),(6,2,2),(6,0,2),(6,3,3),(6,2,3),(6,0,3)]
    kmers_db = KmerData(db_file,k)
    sm = get_significance_model(k,k2,n,kmers_db)
    global ctn += 1
    println("Working with this set", (k,k2,n) )
    for (id,seq) in primer_inf
        parts = []
        push!(parts,id)
        push!(parts,seq)
        push!(parts,k)
        push!(parts,k2)
        push!(parts,n)
        rescal = collect_matches_data(seq,k,k2,n,kmers_db,sm) 
        for r in rescal
            push!(parts,r)
        end 
        push!(data,parts)
    end
    CSV.write("REZ_$ctn.csv",data) 
end 
println(data)