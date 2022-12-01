"
Reads clustering from jc/gc  file. Jc/gc - a tsv file were the firts column - cluster members, second - cluster representative
jc - denotes a clustering after one step. gc - denotes clustering tracking several clustering stages ex. rmident followed by swarm
"
function read_jc_clusters(inf::String)
    map = OrderedDict{String,Array{String,1}}()
    pairs = Array{Tuple{String,String},1}()
    for l in readlines(inf)
        parts = split(l)
        push!(pairs, (parts[1], parts[2]))
    end
    #sort to prevent constant check if key exists
    sort!(pairs, by = x -> x[2])
    cur = ""
    for parts in pairs
        if parts[2] != cur
            map[parts[2]] = Array{String,1}()
            cur = parts[2]
        end
        push!(map[parts[2]], parts[1])
    end
    #sort by the size of cluster
    sort!(map, by = x -> length(map[x]), rev=true)
    return(map)
end

"Convert UC file to a simple tabular format"
function uc2jc(ucinfile::String)
    for l in readlines(ucinfile)
        if l[1] == 'H'
            parts = split(l)
            mem, rep = split(parts[9], ";")[1], split(parts[10], ";")[1]
            println("$mem\t$rep")
        end
        if l[1] == 'C'
            parts = split(l)
            mem = split(parts[9], ";")[1]
            println("$mem\t$mem")
        end
    end
end

"join clustering files of two consecutive clusterings. In essence - the cluster representatives are updated."
function mergejc(f1::String, f2::String)
    println(stderr, "Reading new cluster mapping from $f2")
    map_new = Dict{String,String}()
    
    for l in readlines(f2)
        parts = split(l, "\t")
        map_new[split(parts[1], ";")[1]] = split(parts[2], ";")[1]
    end
    println(stderr, "Readed in new mapping")
    skipped = 0

    for l in readlines(f1)
        parts = split(l, "\t")
        org_name = split(parts[1], ";")[1]
        old_cluster = split(parts[2], ";")[1]
        new_cluster = get(map_new, old_cluster, nothing)
        
        if new_cluster == nothing
            skipped += 1
        else
            println("$org_name\t$new_cluster")
        end
    end
    println(stderr, "$skipped cluster were discarded as were not found in the new mapping")
end

"Return a bool array which length is equal to input OrderedDict keys and false mark small cluster"
function mark_proper_size_clusters(cls::OrderedDict{String,Array{String,1}}, minfrac::Float64)
    cls_sizes = map(length, values(cls))
    all_ct = sum(cls_sizes)
    fracs = (cls_sizes ./ all_ct) .>= minfrac
    out = OrderedDict(zip(keys(cls), fracs))
    return(out)
end 