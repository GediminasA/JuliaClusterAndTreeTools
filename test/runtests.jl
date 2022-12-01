module Testsanitise
    using Pkg
    Pkg.activate(".")
    using Test
    using st16SseqJuliaTools

    #include("AlignmentMatrix.jl")
    include("UCparsing.jl")
    include("HMMutils.jl")   
    include("BLASTparsing.jl")     
    include("TaxonKitWrap.jl")  
    include("Utils.jl")

end 
