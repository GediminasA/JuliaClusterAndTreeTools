using ArgParse2
using DataFrames
using DataStructures
using BioSequences
using CodecZlib
using BioAlignments
using DataFramesMeta
using Statistics
using st16SseqJuliaTools
using ProgressMeter
using XAM

function main()
    parser = ArgumentParser(prog = "analyse_alnself.jl",
                        description = "Detect chimeras",
                        epilog = "Removing potential chimeras and repair snps")
    add_argument!(parser, "--ref", "-r", help = "reference fasta",type = String)
    add_argument!(parser, "--ideal", "-i", help = "ideal contigs - for debugging",type = String,default="")
    add_argument!(parser, "--salmon-rez", "-s", help = "Salmon rezults",type = String)
    add_argument!(parser, "--out", "-o", help = "Fasta file to write out cleaned up sequences",default="chimerasfree.fasta",type = String)
    #add_argument!(parser, "--output", "-o", help = "File with names",type = String)
    args = parse_args(parser)
    salmon_data = read_salmon(args.salmon_rez)
    ref = fasta_to_dict(args.ref)
    length_data=Dict{String,Int64}()
    for k in keys(ref)
        length_data[k] = length(ref[k])
    end 
    
    #Limits for real match ... mayve should be options
    minl = 100
    minid = 0.96
    pretrim_length = 250
    minid_for_close = 0.96 # Identity of first 250 fragments of contigs
    max_unaligned_length_to_consider_conservative_start_extention = 20
    min_dif_in_expression = 0.1
    min_dif_in_length = 0.1

    
    costmodel = BioAlignments.CostModel(match=0, mismatch=1, insertion=1, deletion=1);
    #known_chimeras = ["1_1868;size=92", "1_6093;size=54", "1_2188;size=22","1_1760;size=8"]
    args.ideal == "" ? known_good_contigs = [] : known_good_contigs = collect(keys(fasta_to_dict(args.ideal))) 
    pairs_4chim = Array{Tuple{String,String,Int64,Int64,Float64,Bool},1}()
    pairs_4contain = Array{Tuple{String,String,Int64,Int64,Float64,Float64,Int64},1}()
    aln_costs = AffineGapScoreModel(
        match=5,
        mismatch=-4,
        gap_open=-10,
        gap_extend=-5  
    )
    
    ref_keys = collect(keys(ref))
    @showprogress 1 "Pairwise alignments..."  for k1 in ref_keys
        #TO DO - below there are three consequtive alignments with some similar steps - should be a function there     
        pairs_part_4chim = Array{Tuple{String,String,Int64,Int64,Float64,Bool},1}(undef,length(ref_keys))   
        pairs_part_4contain = Array{Tuple{String,String,Int64,Int64,Float64,Float64,Int64},1}(undef,length(ref_keys))   
        Threads.@threads for i in 1:length(ref_keys)
            k2 = ref_keys[i]
            if k1 != k2
                rn = k1
                tn = tempname = k2
                rs = ref[rn]
                ts = ref[tn]
                ts_wopref = ts[pretrim_length+1:end]
                aln = alignment(pairalign(LocalAlignment(), rs, ts_wopref, aln_costs))
                matches = BioAlignments.count_matches(aln)
                mismatches = BioAlignments.count_mismatches(aln)
                ins = count_insertions(aln)
                del = count_deletions(aln)
                s = pos = first(aln.a.aln.anchors).seqpos+1
                e = rpos = last(aln.a.aln.anchors).seqpos
                l = abs(s-e)+1
                idfrac = matches/l 
                t_pos=first(aln.a.aln.anchors).refpos+1
                t_rpos=last(aln.a.aln.anchors).refpos
                tseq = ref[tempname][pretrim_length+1:end][t_pos:t_rpos]
                rtseq = rs[pos:rpos]
                #check if  close sequences and the detected match could be due to extention of conservative region 
                t_pref = ref[tempname][1:250]
                r_pref = ref[rn][1:250] 
                pairaln_pref = BioAlignments.pairalign(BioAlignments.EditDistance(),t_pref,r_pref,costmodel)
                pairaln_pref = BioAlignments.alignment(pairaln_pref)
                matches_pref = BioAlignments.count_matches(pairaln_pref)
                idfrac_pref = matches_pref/250
                is_cons5region_extention = (idfrac_pref >= minid_for_close && abs(t_pos+pretrim_length-pos)<= max_unaligned_length_to_consider_conservative_start_extention)  
                # Calculate overall alignment
                aln_fulllen = alignment(pairalign(LocalAlignment(), rs, ts, aln_costs))
                matches_fulllen = BioAlignments.count_matches(aln_fulllen)
                mismatches_fulllen = BioAlignments.count_mismatches(aln_fulllen)
                ins_fulllen = count_insertions(aln_fulllen)
                del_fulllen = count_deletions(aln_fulllen)
                s_fulllen = pos_fulllen = first(aln_fulllen.a.aln.anchors).seqpos
                e_fulllen = rpos_fulllen = last(aln_fulllen.a.aln.anchors).seqpos
                t_pos_fulllen=first(aln_fulllen.a.aln.anchors).refpos+1
                t_rpos_fulllen=last(aln_fulllen.a.aln.anchors).refpos
                taln_fulllen = abs(t_pos_fulllen-t_rpos_fulllen)+1
                idfrac_fulllen = matches_fulllen/taln_fulllen 
                
                if  (rn in ["5_20;size=1"]) && (tempname in["5_3133;size=1"])
                     println()
                     println(aln)
                     println()
                     println(aln_fulllen)
                     println(pos_fulllen," ",rpos_fulllen," ",length(rs)," ",length(ts))
                     println((tempname,rn,t_pos_fulllen,t_rpos_fulllen,idfrac_fulllen,taln_fulllen/length(ts),length(rs)))
                end
                if taln_fulllen > 0
                    pairs_part_4contain[i] = (tempname,rn,t_pos_fulllen,t_rpos_fulllen,idfrac_fulllen,taln_fulllen/length(ts),length(rs))
                end
                if idfrac >= minid && 
                    l >= minl
                    if  (rn in ["1_2540;size=1"]) && (tempname in["6_421;size=1"])
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
                        println(aln)
                        println("aln of prefix")
                        println(pairaln_pref)
                        println(idfrac_pref)
                        println(minid_for_close)
                        println(t_pos+pretrim_length-pos)
                        println(max_unaligned_length_to_consider_conservative_start_extention)
                    end 
                    pairs_part_4chim[i] = (tempname,rn,t_pos,t_rpos,idfrac,is_cons5region_extention)
                end
            end
        end
        for ii in 1:length(ref_keys)
            if isassigned(pairs_part_4chim,ii)
                push!(pairs_4chim,pairs_part_4chim[ii])
            end 
            if isassigned(pairs_part_4contain,ii)
                push!(pairs_4contain,pairs_part_4contain[ii])
            end 
        end
    end
    #Collect data for  chimera detection
    tn_rn_map = Dict{String,Array{String,1}}()
    for pair in pairs_4chim
        tn,rn,is_cons5region_extention =  pair[1],pair[2],pair[6]
        if !is_cons5region_extention 
            if !(tn in keys(tn_rn_map)) 
                tn_rn_map[tn] = Array{String,1}()
            end 
            push!(tn_rn_map[tn],rn)
        end
    end
    #detect if chimeric or not
    chimeras = Array{String,1}()
    #println("--------------------------------------------------------------")
    #testf = open("test.csv","w")
    #println(testf,"tn max_r_id  good_t good_r good_both chimeric false_chimeric  tn_s max_r_s difs len_t len_r difl  ")
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
        len_t = length_data[tn]
        len_r = length_data[max_r_id]
        difl = (len_t-len_r)/len_t
        
        chimeric = (tn_s < max_r_s) && (len_t <= len_r)
        
        good_t = (tn in known_good_contigs )
        good_r = (max_r_id in known_good_contigs )
        good_both = good_t && good_r
        difs = abs(max_r_s - tn_s )/tn_s
        false_chimeric = chimeric && good_t
        ct = 0

       # println(testf,"$tn $max_r_id  $good_t $good_r $good_both $chimeric $false_chimeric  $tn_s $max_r_s $difs $len_t $len_r $difl  ")
        if  chimeric 
            println("Template: ID: $tn  TPM: $tn_s  good: $good_t  Most abundant chimeric partner:  ID: $max_r_id   TPM: $max_r_s   good: $good_r")
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
                pairalnres = BioAlignments.pairalign(BioAlignments.EditDistance(),rs,ts[pretrim_length+1:end],aln_costs)
                pairaln = BioAlignments.alignment(pairalnres)
                println(pairaln)
                println("Exiting...")
                exit()
            end
        end

    end
    #close(testf) 
    
    chimeras = unique(chimeras)
    
    #Collect data for  contained detection
    
    contained = Array{String,1}()
    tn_rn_map = Dict{String,Array{Tuple{String,String,Float64,Float64,Int64},1}}()
    for pair in pairs_4contain
        tn,rn,idf,covf,rl =  pair[1],pair[2],pair[5],pair[6],pair[7]
        if !(tn in keys(tn_rn_map)) 
            tn_rn_map[tn] = Array{Tuple{String,String,Float64,Float64,Int64},1}()
        end 
        push!(tn_rn_map[tn],(tn,rn,idf,covf,rl))
    end
    
    #detect contained 
    #merge perfect matches per high coverage
    contained_idf_limit_c = 0.96 # Limid for identity for containment
    contained_covf_limit_c = 0.99 # Limid for coverage for containment
    #merge perfect matches per high identity fragment - merge ideal matches allowing sme artifacts at the ends...
    contained_idf_limit_i = 0.99 # Limid for identity for containment
    contained_covf_limit_i = 0.95 # Limid for coverage for containment
    for ts in tn_rn_map["5_3133;size=1"]
        println(ts)
    end 
    for tn in keys(tn_rn_map)
        tl = length(ref[tn])
        cont = false
        for rd in tn_rn_map[tn]
            if tl <= rd[5]
                if tn == "5_3133;size=1i" && rd[2] == "5_20;size=1"
                    println("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
                    println((round(rd[3],digits=2) >= contained_idf_limit_c && round(rd[4],digits=2) >= contained_covf_limit_c))
                    println((round(rd[3],digits=2), contained_idf_limit_i, round(rd[4],digits=2) ,contained_covf_limit_i))
                    exit()
                end
                if (round(rd[3],digits=2) >= contained_idf_limit_c && round(rd[4],digits=2) >= contained_covf_limit_c) ||
                (round(rd[3],digits=2) >= contained_idf_limit_i && round(rd[4],digits=2) >= contained_covf_limit_i) 
                    if tl == rd[5]
                        #if lengths matches add abundance criteria
                        r_s = salmon_data[split(rd[2],";")[1]]
                        t_s = salmon_data[split(tn,";")[1]]
                        if r_s > t_s
                            cont = true
                            break 
                        end 
                    else 
                        cont = true 
                        break
                    end
                end 
            end
        end
        if cont 
            push!(contained,tn)
        end 
    end 
    filtered_out = OrderedDict{String,LongSequence{BioSequences.DNAAlphabet{4}}}() 
    ct_chimeric = 0
    ct_contained = 0
    ct_good = 0
    for k in keys(ref)
        ischim = (k in chimeras)
        iscont = (k in contained)
        if k == "5_3133;size=1"
            println("AAAAAAAAAAAH ",ischim, " ",iscont)
        end 
        if !ischim && !iscont
            ct_good += 1
            filtered_out[k]= ref[k]
        else
            if ischim 
                ct_chimeric += 1
            end 
            if iscont
                ct_contained += 1
            end 
        end
    end
    dict_to_fasta(filtered_out,args.out)
    println("$ct_chimeric chimeras , $ct_contained contained sequences were detected, $ct_good sequences were  written to ",args.out)
end

main()
