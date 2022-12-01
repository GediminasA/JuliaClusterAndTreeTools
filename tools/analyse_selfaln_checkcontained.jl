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

struct Containment
    template:: String #id of the first sequences
    reference:: String #id of other sequence 
    template_from::Int64
    template_to::Int64
    identity:: Float64 #Aligned part identuty
    coverage:: Float64 #Coverage fraction in trms of template length
    tlength:: Int64 #template length
    rlength:: Int64 #reference length
    endfragid_shorter::Float64
    endfragid_longer::Float64
    aln_fulllen::PairwiseAlignment{LongSequence{DNAAlphabet{4}},LongSequence{DNAAlphabet{4}}} 
    aln_3end_shorter::PairwiseAlignment{LongSequence{DNAAlphabet{4}},LongSequence{DNAAlphabet{4}}} 
    aln_3end_longer::PairwiseAlignment{LongSequence{DNAAlphabet{4}},LongSequence{DNAAlphabet{4}}} 
end

function main()
    parser = ArgumentParser(prog = "analyse_alnself.jl",
                        description = "Detect chimeras",
                        epilog = "Removing potential chimeras and repair snps")
    add_argument!(parser, "--ref", "-r", help = "reference fasta",type = String)
    add_argument!(parser, "--salmon-rez", "-s", help = "Salmon rezults",type = String)
    add_argument!(parser, "--out", "-o", help = "Fasta file to write out cleaned up sequences",default="chimerasfree.fasta",type = String)
    add_argument!(parser, "--fused", "-f", help = "Fasta file to write out cleaned up sequences where contigs of each cluster are joined in one sequence",default="chimerasfree.fasta",type = String)

    
    args = parse_args(parser)
    salmon_data = read_salmon(args.salmon_rez)
    ref = fasta_to_dict(args.ref)
    length_data=Dict{String,Int64}()
    for k in keys(ref)
        length_data[k] = length(ref[k])
    end 
    
    #Limits for real match ... maybe should be options
    endlength1 = 15 #15 #length from end to rtest separately
    endlength2 = 50 #50 #length from end to rtest separately
    endfrag_id_limit = 0.8 #if sequences match general cut-off but the end fragment is dissimilar they are not merged
    contained_idf_limit_c = 0.96 # 0.96 Limid for identity for containment 0.97 #z1 - 7_416;size=1", "4_2136;size=1. Salmonela and coli
    contained_covf_limit_c = 0.99 # Limid for coverage for containment

    
    costmodel = BioAlignments.CostModel(match=0, mismatch=1, insertion=1, deletion=1);
    #known_chimeras = ["1_1868;size=92", "1_6093;size=54", "1_2188;size=22","1_1760;size=8"]
    pairs_4contain = Array{Containment,1}()
    aln_costs = AffineGapScoreModel(
        match=5,
        mismatch=-4,
        gap_open=-10,
        gap_extend=-5  
    )
    aln_costs_wogaps = AffineGapScoreModel(
        match=5,
        mismatch=-4,
        gap_open=-10,
        gap_extend=-5  
    )
    
    ref_keys = collect(keys(ref))
    @showprogress 1 "Pairwise alignments..."  for k1 in ref_keys
        #TO DO - below there are three consequtive alignments with some similar steps - should be a function there     
        pairs_part_4contain = Array{Containment,1}(undef,length(ref_keys))   
        s1 = split(k1,"@")[end]
        Threads.@threads for i in 1:length(ref_keys)
            k2 = ref_keys[i]
            s2 = split(k2,"@")[end]
            #Limit containments analysis to the same species 
            if (k1 != k2) && (s1 == s2)
                rn = k1
                tn = tempname = k2
                rs = ref[rn]
                ts = ref[tn]
                tlen = length(ts)
                rlen = length(rs)

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
                #calculate identity at 3' 
                tendfrag1 = ts[t_rpos_fulllen-endlength1+1:t_rpos_fulllen]
                rendfrag1 = rs[e_fulllen-endlength1+1:e_fulllen]
                tendfrag2 = ts[t_rpos_fulllen-endlength2+1:t_rpos_fulllen]
                rendfrag2 = rs[e_fulllen-endlength2+1:e_fulllen]
                aln_3end1 = alignment(pairalign(GlobalAlignment(), rendfrag1, tendfrag1, aln_costs_wogaps))
                aln_3end2 = alignment(pairalign(GlobalAlignment(), rendfrag2, tendfrag2, aln_costs_wogaps))
                matches_3end1 = BioAlignments.count_matches(aln_3end1)
                matches_3end2 = BioAlignments.count_matches(aln_3end2)
                
                t_pos_3end1=first(aln_3end1.a.aln.anchors).refpos+1
                t_rpos_3end1=last(aln_3end1.a.aln.anchors).refpos
                taln_3end1 = abs(t_pos_3end1-t_rpos_3end1)+1
                idfrac_3end1 = matches_3end1/taln_3end1 
                
                t_pos_3end2=first(aln_3end2.a.aln.anchors).refpos+1
                t_rpos_3end2=last(aln_3end2.a.aln.anchors).refpos
                taln_3end2 = abs(t_pos_3end2-t_rpos_3end2)+1
                idfrac_3end2 = matches_3end2/taln_3end2 
                
                
                if taln_fulllen > 0
                    pairs_part_4contain[i] = Containment(tempname,
                    rn,
                    t_pos_fulllen,
                    t_rpos_fulllen,
                    idfrac_fulllen,
                    taln_fulllen/length(ts),
                    length(ts),
                    length(rs),
                    idfrac_3end1,
                    idfrac_3end2,
                    aln_fulllen,
                    aln_3end1,
                    aln_3end2)
                end
            end
        end
        for ii in 1:length(ref_keys)
            if isassigned(pairs_part_4contain,ii)
                push!(pairs_4contain,pairs_part_4contain[ii])
            end 
        end
    end
    
    #Collect data for  contained detection
    
    contained = Array{String,1}()
    #reformat collected data into dictionary
    tn_rn_map = Dict{String,Array{Containment,1}}()
    for pair in pairs_4contain
        #tn,rn,idf,covf,rl =  pair[1],pair[2],pair[5],pair[6],pair[7]
        if !(pair.template in keys(tn_rn_map)) 
            tn_rn_map[pair.template] = Array{Containment,1}()
        end 
        push!(tn_rn_map[pair.template],pair)
    end
    
    #detect contained 
    
    for tn in keys(tn_rn_map)
        tl = length(ref[tn])
        cont = false
        for rd in tn_rn_map[tn]
            t_s = salmon_data[split(tn,";")[1]]
            r_s = salmon_data[split(rd.reference,";")[1]]
            t_clnr = split(tn,"_")[1]
            r_clnr = split(rd.reference,"_")[1]
            tl_shorter_fraction = abs(rd.rlength-tl)/tl
            if tl <= rd.rlength #&& (t_clnr == r_clnr)
                if (round(rd.identity,digits=2) >= contained_idf_limit_c && 
                    round(rd.coverage,digits=2) >= contained_covf_limit_c &&
                    ( ( (round(rd.endfragid_longer,digits=2) >= endfrag_id_limit) && (round(rd.endfragid_shorter,digits=2) >= endfrag_id_limit)  ) || (tl_shorter_fraction > 0.2) ) ) 
                    if tl == rd.rlength
                        #if lengths matches add abundance criteria
                        if r_s > t_s
                            cont = true
                        end
                    else 
                        cont = true
                    end
                end 
            end
        end
        if cont 
            push!(contained,tn)
        end 
    end 
    filtered_out = OrderedDict{String,LongSequence{BioSequences.DNAAlphabet{4}}}() 
    ct_contained = 0
    ct_good = 0
    for k in keys(ref)
        iscont = (k in contained)
        if  !iscont
            ct_good += 1
            filtered_out[k]= ref[k]
        else
            if iscont
                ct_contained += 1
            end 
        end
    end
    dict_to_fasta(filtered_out,args.out)
    #Join from the same cluster 
    fused_seqs = OrderedDict{String,Array{String,1}}() 
    fused_ids = OrderedDict{String,Array{String,1}}() 
    for k in keys(filtered_out)
        cid = split(k,"_")[1]
        if !(cid in keys(fused_seqs))
            fused_seqs[cid] = Array{String,1}()
            fused_ids[cid] = Array{String,1}()
        end
        push!(fused_seqs[cid],String(filtered_out[k]))
        push!(fused_ids[cid],k)
    end 
    fused_seqnames = OrderedDict{String,LongSequence{DNAAlphabet{4}}}()
    spacer = String(fill('N',300))
    for k in keys(fused_ids) 
        seqname=join(fused_ids[k],";")
        seqname = replace(seqname,";size=1@"=>"|") 
        seqname = replace(seqname,k*"_"=>"") 
        prefix = k*"@"*string(length(fused_ids[k]))*":"
        fused_seqnames[prefix*seqname] = LongSequence{DNAAlphabet{4}}(join(fused_seqs[k],spacer))
    end 
    dict_to_fasta(fused_seqnames,args.fused)

    println(stderr,"$ct_contained contained sequences were detected, $ct_good sequences were  written to ",args.out)
    println(stderr,"Fused sequences are written to ",args.fused, ", name format: <clusterid>@<number of contigs>: ....<contig id|species based on blast>...separated by ;")
end

main()
