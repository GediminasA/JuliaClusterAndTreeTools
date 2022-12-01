using st16SseqJuliaTools
if abspath(PROGRAM_FILE) == @__FILE__
    map = read_jc_clusters(ARGS[1])
    sizes = Array{Tuple{String,Int64},1}()
    for k in keys(map)
        push!(sizes,(k,length(map[k])))
    end
    sort!(sizes,by=x->x[2],rev=true)
    for p in map[sizes[1][1]]
        println(p)
    end
end
