using Test
using st16SseqJuliaTools
using Suppressor

@testset "extract_properly_matching_rRNA_names" begin
    names_file, names_fo = mktemp()
    close(names_fo)
    @suppress_err begin
        @testset "Case: BOUNDARIES" begin
            hmm_str = "100366               -          NAME -              664     932     276       1     280       1     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
            @testset "Case: IN" begin
                input_file, input_fo = mktemp()
                println(input_fo, "#COMMENT")
                println(input_fo, "# COMMENT2")
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "0:500", "600:950", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:500", "600:950", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:280", "600:950", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:280", "664:950", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:280", "664:932", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:280", "664:932", 280, input_file)
            end

            @testset "Case: OUT" begin
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(String[]) == extract_rRNA_names(names_file, "2:500", "600:950", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(String[]) == extract_rRNA_names(names_file, "1:279", "600:950", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(String[]) == extract_rRNA_names(names_file, "1:280", "665:932", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(String[]) == extract_rRNA_names(names_file, "1:280", "664:931", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(String[]) == extract_rRNA_names(names_file, "1:280", "664:932", 281, input_file)
            end
            
            @testset "Case: MIX" begin
                hmm_str2 = "100365               -          NAME2 -              500     700     276       1     300       150     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100366", "100365"]) == extract_rRNA_names(names_file, "1:300", "500:932", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "1:300", "500:700", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "150:300", "500:1000", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:300", "500:932", 270, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:280", "600:932", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                println(input_fo, hmm_str)
                close(input_fo)
                @test Set(["100366"]) == extract_rRNA_names(names_file, "1:280", "600:932", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str2)
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "150:300", "500:1000", 150, input_file)

                input_file, input_fo = mktemp()
                println(input_fo, hmm_str)
                println(input_fo, hmm_str2)
                println(input_fo, hmm_str)
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100366", "100365"]) == extract_rRNA_names(names_file, "1:300", "500:932", 150, input_file)
            end

            @testset "Case: POS MIX" begin
                hmm_str2 = "100365               -          NAME2 -              700     500     276       1     300       150     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "150:300", "500:700", 150, input_file)

                hmm_str2 = "100365               -          NAME2 -              500     700     276       1     150       300     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "150:300", "500:700", 150, input_file)

                hmm_str2 = "100365               -          NAME2 -              700     500     276       1     150       300     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "150:300", "500:700", 150, input_file)

                hmm_str2 = "100365               -          NAME2 -              700     500     276       1     150       300     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "300:150", "500:700", 150, input_file)

                hmm_str2 = "100365               -          NAME2 -              700     500     276       1     150       300     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "300:150", "700:500", 150, input_file)

                hmm_str2 = "100365               -          NAME2 -              700     500     276       1     150       300     280    -     1.1e-83  279.3   8.9  reference=Bacillus_subtilis_:2949495-2951495(+)/135-2000_1591 amplicon=1..726 position=1..300"
                input_file, input_fo = mktemp()
                println(input_fo, hmm_str2)
                close(input_fo)
                @test Set(["100365"]) == extract_rRNA_names(names_file, "150:300", "700:500", 150, input_file)
            end
        end
    end
end