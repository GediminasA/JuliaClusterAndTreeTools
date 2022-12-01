using st16SseqJuliaTools
using Artifacts
using Test
@testset "taxonkit" begin
    db = artifact"taxdump"
    ids = ["535026","595","34","1350","81852","703612"]
    @test ["Bacillus", "Salmonella", "Myxococcus", "Enterococcus", "NA","Bacillus"] == taxid_to_genus(ids,db)
    @test ["Bacillus subtilis", "Salmonella enterica", "Myxococcus xanthus", "NA", "NA", "Bacillus subtilis"] == taxid_to_species(ids,db)
end