
using st16SseqJuliaTools
using ArgParse
using CSV
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using FASTX 

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--16sr", "-r"
            help = "16S expected region in a read, <start>:<end> "
            arg_type = String
            required = true 
        "--16st", "-t"
            help = "16S expected region in a hmm/alignment, <start>:<end> "
            arg_type = String
            required = true 
        "--min", "-m"
            help = "minimum match length "
            arg_type = Int 
            required = true 
        "--input","-i"
            help = "Input file"
            required = true 
        "--log", "-l"
            help = "Log file to write reads on target"
            arg_type = String
            required = true 
        "--names", "-n"
            help = "Write out names of reads "
            arg_type = String
            required = true 
        "--fastqin", "-q"
            help = "FASTQ.GZi input to filter out reads"
            arg_type = String
            required = true 
        "--fastqout", "-o"
            help = "FATQ.GZ output to filter out reads"
            arg_type = String
            required = true 
    end
    return parse_args(s)
end



function write_fasta(log_file::String, fastq_file::String, fastq_fileo::String, chosen_reads::Set{String})
    logf = open(log_file, "w")
    ct_t = length(chosen_reads)
    ct_p = 0
    ct_all = 0
    #read in fastq and filter
    reader = FASTQ.Reader(GzipDecompressorStream(open(fastq_file)))
    writer = FASTQ.Writer(GzipCompressorStream(open(fastq_fileo,"w")))
    println("Choosing reads with proper 16S fragments in file: $fastq_file")
    record = FASTQ.Record()
    
    while !eof(reader)
        read!(reader, record)
        ct_all += 1
        #println(record)
        id = FASTQ.identifier(record)
        
        if id in chosen_reads 
            ct_p += 1
            write(writer, record) 
            delete!(chosen_reads, id)
        end
        ## Do something.
    end
    
    if ct_p != ct_t
        println(stderr, "Read counts are not equal: $ct_t in csv and $ct_p in the fastq filr")
    end
    frac = ct_p/ct_all
    println(stderr, "On 'target' $frac reads")
    println(logf, "$frac")
    close(writer)
    close(logf)
end


if abspath(PROGRAM_FILE) == @__FILE__
    parsed_args = parse_commandline()
    chosen_reads = extract_rRNA_names(parsed_args["names"], parsed_args["16sr"], parsed_args["16st"], parsed_args["min"], parsed_args["input"])
    write_fasta(parsed_args["log"], parsed_args["fastqin"], parsed_args["fastqout"], chosen_reads)
end