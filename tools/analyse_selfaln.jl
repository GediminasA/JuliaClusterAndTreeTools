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
    add_argument!(parser, "--salmon-rez", "-s", help = "Salmon rezults",type = String)
    add_argument!(parser, "--out", "-o", help = "Fasta file to write out cleaned up sequences",default="chimerasfree.faSTA",type = String)
    #add_argument!(parser, "--output", "-o", help = "File with names",type = String)
    args = parse_args(parser)
    salmon_data = read_salmon(args.salmon_rez)
    ref = fasta_to_dict(args.ref)
    
    #Limits for real match ... mayve should be options
    minl = 100
    minid = 0.96
    pretrim_length = 250
    minid_for_close = 0.96 # Identity of first 250 fragments of contigs
    max_unaligned_length_to_consider_conservative_start_extention = 5

    
    costmodel = BioAlignments.CostModel(match=0, mismatch=1, insertion=1, deletion=1);
    reader = open(BAM.Reader,args.inbam) #,index = args.inbam*".bai")
    #known_chimeras = ["1_1868;size=92", "1_6093;size=54", "1_2188;size=22","1_1760;size=8"]
    args.ideal == "" ? known_good_contigs = [] : known_good_contigs = collect(keys(fasta_to_dict(args.ideal))) 

    #intr = ["7_4;size=170", "7_169;size=1344"]
    pairs = Array{Tuple{String,String},1}()
    for record in reader
        if BAM.ismapped(record) #&& SAM.isprimary(record)
            rn = BAM.refname(record)
	        tempname = BAM.tempname(record)
            pos = BAM.position(record)
            seq = BAM.sequence(record)
            rpos = BAM.rightposition(record)
            aln = BAM.alignment(record)
            s = first(aln.anchors).refpos
            e = last(aln.anchors).refpos
            l = abs(s-e)
            qual = BAM.mappingquality(record)
            genus_t = split(tempname,"_")[1]
            genus_r = split(rn,"_")[1]
            #if (tempname != rn) && (genus_t != genus_r) && (genus_t == "1") && (tempname == "1_1160;size=1896") # (tempname == "1_89;size=8") #&& (tempname in known_chimeras) #&& (tempname=="1_1868;size=92")
            t_pos=BioAlignments.ref2seq(aln,pos)[1]
            t_rpos=BioAlignments.ref2seq(aln,rpos)[1]
            rtseq = ref[rn][pos:rpos]
            tseq = ref[tempname][pretrim_length+1:end][t_pos[1]:t_rpos[1]]
            
            #check a match of alignment in bam
            pairalnres = BioAlignments.pairalign(BioAlignments.EditDistance(),rtseq,tseq,costmodel)
            pairaln = BioAlignments.alignment(pairalnres)
            matches = BioAlignments.count_matches(pairaln)
            idfrac = matches/l

            #check if  close sequences 
            t_pref = ref[tempname][1:250]
            r_pref = ref[rn][1:250] 
            pairaln_pref = BioAlignments.pairalign(BioAlignments.EditDistance(),t_pref,r_pref,costmodel)
            pairaln_pref = BioAlignments.alignment(pairaln_pref)
            matches_pref = BioAlignments.count_matches(pairaln_pref)
            idfrac_pref = matches_pref/250

            is_cons5region_extention = (idfrac_pref >= minid_for_close && abs(t_pos+pretrim_length-pos)<= max_unaligned_length_to_consider_conservative_start_extention)  
            if idfrac >= minid && 
                l >= minl &&
                !is_cons5region_extention
                if  (tempname in ["6_12;size=1"]) #&& (rn in["1_809;size=1 "])
                    println("UUUUUUU ",is_cons5region_extention)
                    println(idfrac_pref,"   ",idfrac," AAAAAAA")
                    println("********************************************")
                    println("PAIR: \n>r\n$rtseq\n>t\n$tseq\n---------------")
                    println(">tf")
                    println(ref[tempname][pretrim_length+1:end])
                    println(">rf")
                    println(ref[rn])
                    println(">tf_250")
                    println(ref[tempname][1:250])
                    println(">rf_250")
                    println(ref[rn][1:250])
                    println(">tf_full_$tempname")
                    println(ref[tempname])
                    println(">rf_full_$rn")
                    println(ref[rn])
                    println(tempname," ",rn," ",t_pos," ",t_rpos," ",pos," ",rpos)
                    println("AAAAAAAAAAAAAAAAAAAaaa",idfrac)
                    println(pairaln)
                end 
                push!(pairs,(tempname,rn))
            end
        end
    end
    close(reader)
    #testset = unique(testset)
    #group matches based on template
    tn_rn_map = Dict{String,Array{String,1}}()

    for pair in pairs
        tn,rn =  pair
        if !(tn in keys(tn_rn_map))
            tn_rn_map[tn] = Array{String,1}()
        end
        push!(tn_rn_map[tn],rn)
        s_t = salmon_data[split(tn,";")[1]]
        s_r = salmon_data[split(rn,";")[1]]
        tn_good = (tn in known_good_contigs)
        #println(" $tn $rn $s_t $s_r $tn_good  ")
    end
    #detect if chimeric or not
    chimeras = Array{String,1}()
    intr = Array{String,1}()
    #println("--------------------------------------------------------------")
    for tn in keys(tn_rn_map)
        tn_s = salmon_data[split(tn,";")[1]]
        rn_ss = Array{Float64,1}()
        for rn in tn_rn_map[tn]
            rn_s = salmon_data[split(rn,";")[1]]
            push!(rn_ss, rn_s)
        end
        pot_chim = DataFrame(ref_names = tn_rn_map[tn],ref_salmon = rn_ss)
        sort!(pot_chim,[:ref_salmon],rev=true)
        max_r_s = pot_chim[1,:].ref_salmon
        max_r_id = pot_chim[1,:].ref_names
        chimeric = false
        if tn_s < max_r_s
            chimeric = true
        end
        good_t = (tn in known_good_contigs )
        good_r = (max_r_id in known_good_contigs )
        ct = 0
        println("
        Template: 
            ID: $tn 
            TPM: $tn_s
            good: $good_t 
        Most abundant chimeric partner:
            ID: $max_r_id
            TPM: $max_r_s
            good: $good_r
        ")
         if  chimeric 
            ct += 1
            push!(chimeras,tn)
            if chimeric && good_t  
                println("##############MISCLASSIFICATION##############")
                println("
                Template: 
                    ID: $tn 
                    TPM: $tn_s
                    good: $good_t 
                Most abundant chimeric partner:
                    ID: $max_r_id
                    TPM: $max_r_s
                    good: $good_r
                ")
                ts = ref[tn]
                rs = ref[max_r_id]
                pairalnres = BioAlignments.pairalign(BioAlignments.EditDistance(),ts,rs,costmodel)
                pairaln = BioAlignments.alignment(pairalnres)
                println(pairaln)
                println("Exiting...")
                exit()
            end
        end

    end 
    println("##############RECHECK##############3")
    #println(intr)
    chimeras = unique(chimeras)
    filtered_out = OrderedDict{String,LongSequence{BioSequences.DNAAlphabet{4}}}() 
    ct = 0
    ct_good = 0
    for k in keys(ref)
        if ! (k in  chimeras)
            ct_good += 1
            filtered_out[k]= ref[k]
        else
            ct += 1
        end
    end
    dict_to_fasta(filtered_out,args.out)
    println("$ct chimeras were detected, $ct_good sequences were  written to ",args.out)
end

main()
