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

function main()
    parser = ArgumentParser(prog = "analyse_alnself.jl",
                        description = "Detect chimeras",
                        epilog = "Removing potential chimeras and repair snps")
    add_argument!(parser, "--inbam", "-b", help = "input bam ",type = String)
    add_argument!(parser, "--ref", "-r", help = "reference fasta",type = String)
    add_argument!(parser, "--ideal", "-i", help = "ideal contigs - for debugging",type = String,default="")
    add_argument!(parser, "--out", "-o", help = "Fasta file to write out cleaned up sequences",default="chimerasfree.faSTA",type = String)
    add_argument!(parser, "--min-seqid", "-s", help = "Minimum sequence identity in match",type = Float64, default=0.95)
    add_argument!(parser, "--min-covid", "-c", help = "Minimum sequence length fraction in  match",type = Float64, default=0.96)
    #add_argument!(parser, "--output", "-o", help = "File with names",type = String)
    args = parse_args(parser)
    ref = fasta_to_dict(args.ref)
    
    #Limits for real match ... mayve should be options
    minl = 300
    minf_len = args.min_covid
    minf_id = args.min_seqid # Identity of first 250 fragments of contigs
    reader = open(BAM.Reader,args.inbam) #,index = args.inbam*".bai")
    args.ideal == "" ? known_good_contigs = [] : known_good_contigs = collect(keys(fasta_to_dict(args.ideal))) 
    pairs = Array{Tuple{String,String},1}()
    chosenseqs = OrderedDict{String,BioSequences.LongSequence{BioSequences.DNAAlphabet{4}}}()
    ct_good = 0
    ct_aligned = 0
    ct_total = 0
    for record in reader
        ct_total += 1

        if BAM.ismapped(record) && BAM.isprimary(record)  #&& BAM.tempname(record)=="6_674;size=1" 
            ct_aligned += 1 
            rn = BAM.refname(record)
	        tempname = BAM.tempname(record)
            pos = BAM.position(record)
            seq = BAM.sequence(record)
            rpos = BAM.rightposition(record)
            aln = BAM.alignment(record)
            l_total = length(seq)
            l_aligned = abs(rpos-pos)
            mutations = call_variants(record)
            mismct = 0
            for m in mutations
                if m.class == 'M'
                    mismct += 1
                    #mismct += maximum([length(m.reference),length(m.variant)])
                end
                if m.class == 'D'
                    l_aligned -= length(m.reference)
                end
            end
            fl = round(l_aligned/l_total,digits=2)
            fid = round(1 - (length(mutations)/l_aligned),digits=2)
            #println(fl," ",fid," ")
            good = (tempname in known_good_contigs)
            if l_aligned >= minl && ( (fid >= minf_id && fl >= minf_len) || (fid >= 0.99 && fl >= 0.95) )            
                ct_good += 1
                chosenseqs[tempname] = ref[tempname]
            end 
        end 
    end 

    close(reader)
    println("$ct_total in total sequences, $ct_aligned  were found in alignments,  $ct_good sequences were  written to ",args.out)
    dict_to_fasta(chosenseqs,args.out)
end

main()
