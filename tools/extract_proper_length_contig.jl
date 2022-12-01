using st16SseqJuliaTools
using ArgParse2
using Statistics
using DataStructures
using BioSequences 
using BioAlignments 



function julia_main()::Cint
    parser = ArgumentParser(prog = "",
                            description = "Filter out contigs based on length median",
                            epilog = "")

    add_argument!(parser, "-i", "--infasta", type = String, help = "Input fasta")
    add_argument!(parser, "-o", "--outfasta", type = String, help = "Output fasta")
    add_argument!(parser, "-e", "--length-cutoff", type = Float64, help = "Fraction from maximumlength - shorter contigs will be discarded", default=0.6)
    add_argument!(parser, "-d", "--identity-cutoff", type = Float64, help = "Aligned fraction identity cut off for containment", default=0.97)
    add_argument!(parser, "-v", "--coverage-cutoff", type = Float64, help = "Aligned part length fraction in the shorter sequence", default=0.6)
    add_argument!(parser, "-l", "--length-3endcut", type = Int64, help = "How many 3end bases should be removed before merging", default=20)
    add_argument!(parser,
    "--filter-length",
    action = "store_true",
    default = false,
    help = "Turn on length filtering based on maximumlength") 
    add_argument!(parser,
    "--remove-contained",
    action = "store_true",
    default = false,
    help = "Turn on removal of contained sequences") 
    
    args = parse_args(parser)
    seqs = st16SseqJuliaTools.fasta_to_dict(args.infasta)
    lens = map(x -> length(x),collect(values(seqs)))
    llimit = args.length_cutoff*maximum(lens)
    seqs = sort(seqs, by = x -> length(seqs[x]),rev = true)
    #remove some of 3end
    tmp = OrderedDict{String,BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}()
    for k in keys(seqs)
        tmp[k] = seqs[k][1:end-args.length_3endcut]
    end
    seqs = tmp 
    
    if args.filter_length
        inict = length(seqs)
        seqs = filter(x -> length(x[2]) >= llimit, seqs)
        gotct = length(seqs)
        dct = inict - gotct
        println(stderr,"Removed $dct sequences using length cutoff")
    end
    #len_filtered = filter(x -> len_sorted[x] >= llimit, len_sorted)
    #st16SseqJuliaTools.dict_to_fasta(len_sorted,args.outfasta)
    ks = collect(keys(seqs))
    len_ks = length(ks)
    out = OrderedDict{String,BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}()
    out[ks[1]] = seqs[ks[1]] #longest sequence is always kept
    contained = fill(false,len_ks)
    if args.remove_contained
        for i in 1:len_ks
            for j in i:len_ks
                if i != j
                    identityf,coveragef =  check_pair(seqs[ks[i]], seqs[ks[j]])
                    if round(identityf,digits=2) >= args.identity_cutoff && round(coveragef,digits=2) >= args.coverage_cutoff 
                        contained[j] = true
                        println("aaaa")
                    end
                end
            end
        end
        contained_nr = count(contained)
        println(stderr,"Removed $contained_nr due to containements")
        for i in 1:len_ks
            k = ks[i]
            cnted = contained[i]
            if !cnted
                out[k]=seqs[k]
            end 
        end
        seqs = out 
    end 
    st16SseqJuliaTools.dict_to_fasta(seqs,args.outfasta)
   return 0


end

julia_main()