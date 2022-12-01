using st16SseqJuliaTools
using Test
using BioTools.BLAST
using BioSequences
using IntervalSets
@testset "Blast results manipulation" begin 
    playxml = "test/playdata/chimrem/blast_results.xml"
    data = split_by_query(readblastXML(read(playxml, String)))
    @test [10..16, 17..19] == add_new_interval([10..12,15..16,17..19],12..15)
    @test [10..12, 13..14, 15..16, 17..19] == add_new_interval([10..12,15..16,17..19],13..14)
    @test covered_length([10..16,17..19]) == 10
    @test 0.0 == unique_length_fraction([10..12, 13..14, 15..16, 17..19],12..13)
    @test 0.5 == unique_length_fraction([10..12, 13..14, 15..16, 17..19],19..20)
    @test [("gi|1384345818|gb|CP029032.1|", 556866, 556195), ("gi|1384345818|gb|CP029032.1|", 1997488, 1997716)] == map(x -> (x.hitid,x.hitfrom,x.hitto),get_coverage_maximising_hits("1_1828;size=1",data))
end   