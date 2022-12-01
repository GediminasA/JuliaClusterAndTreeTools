using CSV
using DataFrames
using DataFramesMeta
using DataStructures
using CodecZlib
using XAM
using FASTX
using BioSequences
using BioAlignments
using ArgParse2
using st16SseqJuliaTools


"Parse genus asignment - specific for fasta header"
function parse_genus_info(inf::String)
    out = Dict{String,String}()
    df = DataFrame(CSV.File(inf))
    for r in eachrow(df)
	out[string(r.ID)] = split(split(r.Genus,"-")[1],"_")[1]
    end
    return(out)
end



function main()
    parser = ArgumentParser(prog = "analyse_alignment_on_reference.jl",
                        description = "Extract good contigs based on expected contigs of a standard")
    add_argument!(parser, "--inbam", "-i", help = "input sam ",type = String)
    add_argument!(parser, "--output", "-o", help = "File with names",type = String)
    add_argument!(parser, "--ref","-r", help = "Reference",type = String)
    add_argument!(parser, "--contig","-c", help = "Contig",type = String)
    add_argument!(parser, "--genus","-g", help = "Genus",type = String)
    add_argument!(parser, "--contig-fraction","-f", help = "Fraction of a contig that should be aligned",type = Float64,default = 0.7)
    add_argument!(parser, "--contig-identity","-d", help = "Minimum identity of an aligned contig",type = Float64,default = 0.7)
    args = parse_args(parser)
    #Not command line parameters
    min_al_idfrac = args.contig_identity # length of contig aligned
    min_al_lenfrac = args.contig_fraction  # identity fraction
    write_ref_seq = false # write matching reference sequence instead of a contig - for testing
    
    #dat of all matches on reference will be collected for filtering
    matches = DataFrame(Contig_ID = Array{String,1}(),
		        Ref_ID = Array{String,1}(),
		        Contig_length = Array{Int64,1}(),
			Match_length = Array{Int64,1}(),
			Identity = Array{Float64,1}(),
			Coverage_contig = Array{Float64,1}(),
			Coverage_ref = Array{Float64,1}(),
			Contig_start = Array{Int64,1}(), 
			Contig_end = Array{Int64,1}(), 
			Ref_start = Array{Int64,1}(), 
			Ref_end = Array{Int64,1}(),
			Ref_genus = Array{String,1}(),
			Expected_genus = Array{Bool,1}())
			
    out_seqs = OrderedDict{String,LongSequence{DNAAlphabet{4}}}()
    genus_assigns = parse_genus_info(args.genus)
    costmodel = BioAlignments.CostModel(match=0, mismatch=1, insertion=1, deletion=1);
    testseqs =  fasta_to_dict(args.contig)
    ref =  fasta_to_dict(args.ref)
    print(args.inbam)	
    reader = open(BAM.Reader,args.inbam)
    ct = 0 
    
    for record in reader
        if BAM.ismapped(record) #&& SAM.isprimary(record)
            ct += 1
	    pos = BAM.position(record)
            seq = BAM.sequence(record)
            rpos = BAM.rightposition(record)
            aln = XAM.BAM.alignment(record)
            s = first(aln.anchors).refpos
            e = last(aln.anchors).refpos
            l = abs(s-e)
            rn = BAM.refname(record)
	    qual = XAM.BAM.mappingquality(record)
	    tempname = BAM.tempname(record)
	    genus = genus_assigns[split(tempname,"_")[1]]
	    real_length = length(testseqs[BAM.tempname(record)])
	    tempname = BAM.tempname(record)
	    t_pos=BioAlignments.ref2seq(aln,pos)
	    t_rpos=BioAlignments.ref2seq(aln,rpos)
	    rtseq = ref[rn][pos:rpos]
	    tseq = testseqs[tempname][t_pos[1]:t_rpos[1]]
	    ff = findfirst(lowercase(genus),lowercase(rn))
	    expected_g = false
	    if ff != nothing && ff.start == 1
		expected_g = true
	    end




	    
	    #println("PAIR: \n>r\n$rtseq\n>t\n$tseq\n---------------")
	    pairalnres = BioAlignments.pairalign(BioAlignments.EditDistance(),rtseq,tseq,costmodel)
	    pairaln = BioAlignments.alignment(pairalnres)
	    matched = BioAlignments.count_matches(pairaln)
	    idfrac = matched/l
	    lfrac_c = l/real_length
	    lfrac_r = l/length(ref[rn])
	    df_row=[tempname,rn,real_length,l,idfrac,lfrac_c,lfrac_r,t_pos[1],t_rpos[1],pos,rpos,genus,expected_g]
	    push!(matches,df_row)
	    
        end
    end
    close(reader)
    matches_expected = @where(matches, :Expected_genus .== 1)
    ref_c_ds = groupby(matches_expected,:Ref_ID)
    lens =Array{Int64,1}()
    chosen_ids = Array{String,1}()
    for ref_c_d in ref_c_ds
	df = DataFrame(ref_c_d)
	df = @where(df,:Coverage_contig .>= min_al_lenfrac, :Identity .>= min_al_idfrac)
	df = sort(df, :Match_length,rev=true)
	s = df[1,:]
	push!(chosen_ids, s.Contig_ID)
	if write_ref_seq
	    out_seqs[s.Contig_ID] = LongSequence{DNAAlphabet{4}}(ref[s.Ref_ID][s.Ref_start:s.Ref_end])
	else
	    out_seqs[s.Contig_ID] = LongSequence{DNAAlphabet{4}}(testseqs[s.Contig_ID])#[s.Contig_start:s.Contig_end]
	end
	push!(lens,length(testseqs[s.Contig_ID]))
    end
    dict_to_fasta(out_seqs,args.output)#,cut_end = minimum(lens))
    Proper = Array{Bool,1}()
    for i in matches.Contig_ID 
	push!(Proper,(i in chosen_ids))
    end 
    matches.Proper = Proper 
    CSV.write(args.output*".info.csv", matches,delim=';')
end

main()
