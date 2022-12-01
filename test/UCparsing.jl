using Test
using st16SseqJuliaTools
using Suppressor
using DataStructures

@testset "uc2jc" begin
    @testset "Case: C" begin
        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n959\t959\n959\t959\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959\t*\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;AAAAA;AAAAA\t*\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;AAAAA;AAAAA\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n959\t959\n959\t959\n"
    end

    @testset "Case: H" begin
        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n975\t838\n975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=6629\t838\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975\t838\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;cacac;acaac\t838;acac;c;a;ca\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;cacac;acaac\t838;acac;c;a;ca\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975\t838\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=6629\t838\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n975\t838\n975\t838\n975\t838\n975\t838\n"
    end

    @testset "Case: MIXED" begin

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "975\t838\n959\t959\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n975\t838\n959\t959\n975\t838\n959\t959\n975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n959\t959\n959\t959\n959\t959\n975\t838\n975\t838\n975\t838\n975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "S\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "S\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "S\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "S\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "S\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n959\t959\n959\t959\n959\t959\n975\t838\n975\t838\n975\t838\n975\t838\n"

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "X\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "A\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "C\t0\t10115\t*\t*\t*\t*\t*\t959;size=10115\t*\n")
        print(tmp_fo, "G\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "S\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "A\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        print(tmp_fo, "A\t6\t240\t*\t*\t*\t*\t*\t952;size=5124\t*\n")
        print(tmp_fo, "H\t2\t240\t99.6\t+\t0\t0\t240M\t975;size=2856\t838;size=6629\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == "959\t959\n959\t959\n959\t959\n959\t959\n975\t838\n975\t838\n975\t838\n"
    end

    @testset "Case: ANOMALY" begin
        tmp_file, tmp_fo = mktemp()
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == ""

        tmp_file, tmp_fo = mktemp()
        print(tmp_fo, "aCACACACACACA\n")
        close(tmp_fo)
        output = @capture_out begin
            uc2jc(tmp_file)
        end
        @test output == ""

        tmp_file = "AAAAA"
        @test_throws SystemError uc2jc(tmp_file)
    end

end

@testset "mergejc" begin
    @testset "Case: MERGE" begin
        output = ""
        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "4956\t14932\n14980\t14932\n13426\t14932\n")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "14932\t14932\n")
        close(tmp_fo2)
        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '0'
        @test output == "4956\t14932\n14980\t14932\n13426\t14932\n"

        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "4956\t14931\n14980\t14932\n13426\t14932\n")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "14931\t14932\n")
        close(tmp_fo2)

        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '2'
        @test output == "4956\t14932\n"

        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "4956\t14931\n14980\t14932\n13426\t14932\n")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "14931\t14932\n14932\t14933\n")
        close(tmp_fo2)

        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '0'
        @test output == "4956\t14932\n14980\t14933\n13426\t14933\n"

        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "4956\t14931\n14980\t14932\n13426\t14932\n")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "14931\t14932\n14932\t14933\n14930\t14938\n")
        close(tmp_fo2)

        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '0'
        @test output == "4956\t14932\n14980\t14933\n13426\t14933\n"

        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "4956\t14931\n14980\t14932\n13426\t14932\n13450\t14955\n")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "14931\t14932\n14932\t14933\n14930\t14938\n")
        close(tmp_fo2)

        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '1'
        @test output == "4956\t14932\n14980\t14933\n13426\t14933\n"
    end

    @testset "Case: ANOMALY" begin
        output = ""
        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "")
        close(tmp_fo2)

        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '0'
        @test output == ""

        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "4956\t14931\n14980\t14932\n13426\t14932\n13450\t14955\n")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "")
        close(tmp_fo2)

        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '4'
        @test output == ""

        tmp_file1, tmp_fo1 = mktemp()
        print(tmp_fo1, "")
        close(tmp_fo1)
        tmp_file2, tmp_fo2 = mktemp()
        print(tmp_fo2, "14931\t14932\n")
        close(tmp_fo2)

        err = @capture_err begin
            output = @capture_out begin
                mergejc(tmp_file1, tmp_file2)
            end
        end
        @test split(err, "\n")[3][1] == '0'
        @test output == ""

        tmp_file1 = "NO_file"
        tmp_file2 = "NO_file2"
        @suppress_err begin
            @test_throws SystemError mergejc(tmp_file1, tmp_file2)
        end
    end
end

@testset "mark_proper_size_clusters" begin
    test_dict = OrderedDict{String,Array{String,1}}()
    test_dict["Clus1"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus2"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus3"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus4"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus5"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    true_dict = OrderedDict{String,Bool}("Clus1" => 0,"Clus2" => 0,"Clus3" => 0,"Clus4" => 0,"Clus5" => 0)
    @test mark_proper_size_clusters(test_dict, 0.5) == true_dict

    test_dict = OrderedDict{String,Array{String,1}}()
    test_dict["Clus1"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus2"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus3"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus4"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus5"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    true_dict = OrderedDict{String,Bool}("Clus1" => 1,"Clus2" => 1,"Clus3" => 1,"Clus4" => 1,"Clus5" => 1)
    @test mark_proper_size_clusters(test_dict, 0.1) == true_dict

    test_dict = OrderedDict{String,Array{String,1}}()
    test_dict["Clus1"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus2"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus3"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus4"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus5"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    true_dict = OrderedDict{String,Bool}("Clus1" => 1,"Clus2" => 1,"Clus3" => 1,"Clus4" => 1,"Clus5" => 1)
    @test mark_proper_size_clusters(test_dict, 0.2) == true_dict

    test_dict = OrderedDict{String,Array{String,1}}()
    test_dict["Clus1"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus2"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus3"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus4"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    test_dict["Clus5"] = ["ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB", "ABCDABCDAB"]
    true_dict = OrderedDict{String,Bool}("Clus1" => 1,"Clus2" => 0,"Clus3" => 1,"Clus4" => 0,"Clus5" => 1)
    @test mark_proper_size_clusters(test_dict, 0.2) == true_dict
end

@testset "read_jc_clusters" begin
    tmp_file1, tmp_fo1 = mktemp()
    print(tmp_fo1, "4956\t14931\n14980\t14932\n13426\t14932\n13450\t14955\n")
    close(tmp_fo1)
    true_dict = OrderedDict("14932" => ["14980", "13426"],"14931" => ["4956"],"14955" => ["13450"])
    @test read_jc_clusters(tmp_file1) == true_dict

    tmp_file1, tmp_fo1 = mktemp()
    print(tmp_fo1, "4956\t14931\n14980\t14932\n13426\t14932\n13450\t14955\n")
    close(tmp_fo1)
    true_dict = OrderedDict("14932" => ["14980", "13426"],"14931" => ["4956"],"14955" => ["13450"])
    @test read_jc_clusters(tmp_file1) == true_dict

    tmp_file1, tmp_fo1 = mktemp()
    print(tmp_fo1, "4956\t14931\n14980\t14932\n13450\t14955\n")
    close(tmp_fo1)
    true_dict = OrderedDict("14932" => ["14980"],"14931" => ["4956"],"14955" => ["13450"])
    @test read_jc_clusters(tmp_file1) == true_dict

    tmp_file1, tmp_fo1 = mktemp()
    print(tmp_fo1, "4956\t14931\n4957\t14931\n4958\t14931\n4959\t14931\n")
    close(tmp_fo1)
    true_dict = OrderedDict("14931" => ["4956", "4957", "4958", "4959"])
    @test read_jc_clusters(tmp_file1) == true_dict

    tmp_file1, tmp_fo1 = mktemp()
    print(tmp_fo1, "4956\t14931\n4957\t14931\n4958\t14931\n4959\t14931\n4956\t14932\n4957\t14932\n4958\t14932\n")
    close(tmp_fo1)
    true_dict = OrderedDict("14931" => ["4956", "4957", "4958", "4959"], "14932" => ["4956", "4957", "4958"])
    @test read_jc_clusters(tmp_file1) == true_dict

    tmp_file1, tmp_fo1 = mktemp()
    print(tmp_fo1, "4956\t14932\n4957\t14932\n4958\t14932\n4956\t14931\n4957\t14931\n4958\t14931\n4959\t14931\n")
    close(tmp_fo1)
    true_dict = OrderedDict("14931" => ["4956", "4957", "4958", "4959"], "14932" => ["4956", "4957", "4958"])
    @test read_jc_clusters(tmp_file1) == true_dict

end