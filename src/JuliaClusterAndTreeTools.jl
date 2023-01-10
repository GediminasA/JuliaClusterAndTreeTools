__precompile__() 
module JuliaClusterAndTreeTools
    using Statistics
    using StatsBase 
    using FastaIO
    using AutoHashEquals
    using DataStructures
    using DataFrames
    using Conda
    using BioSequences
    using BioSymbols
    using FASTX
    using JLD2
    using Distributions
    using IntervalSets
    using ProgressMeter
    using Artifacts
    using BioAlignments
    using GZip

    #Alignment analysis stuff 
    export Alignment
    export Alignmentprot
    export read_fasta
    export write_to_fasta
    export sub_alignment
    export sub_alignment!
    export sub_alignment_rowwise!
    export sub_alignment_columnwise!
    export vsearh_cluster!
    export get_max_cl_length
    export char_array_entropy 
    export aln_conservation
    export eval_recomb


    #Export TaxonKitWrap
    export taxid_to_species
    export taxid_to_genus 

    #HMM parsing 
    export extract_rRNA_names 

    #export workflow specific staff
    export parse_mapping
    export parse_dada_results
    export check_pair
    export fasta_to_dict 
    export dict_to_fasta
    export read_salmon

    #export clustering output parsing functinos
    export read_jc_clusters
    export uc2jc
    export mergejc
    export mark_proper_size_clusters


    include("utils.jl"	)
    include("TaxonKitWrap.jl")
    include("HMMutils.jl")
    include("UCparsing.jl")
    include("AlignmentMatrix.jl")
end
