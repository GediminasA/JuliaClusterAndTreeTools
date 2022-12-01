
"Filters HMM outout and filters out desired matches"
function extract_rRNA_names(names_file::String, r16s::String, t16s::String, minl::Int, input::String)
    namesf = open(names_file, "w")

    # get read limit
    r_data = split(r16s, ":")
    r_star = parse(Int64, r_data[1])
    r_end = parse(Int64, r_data[2])
    println(stderr, "Filtering out matches in reads between $r_star:$r_end")
    
    # get red limit
    t_data = split(t16s, ":")
    t_star = parse(Int64, t_data[1])
    t_end = parse(Int64, t_data[2])
    t_star, t_end = sort([t_star, t_end])
    r_star, r_end = sort([r_star, r_end])
    println(stderr, "Filtering out matches in template between $t_star:$t_end")
    chosen_reads = Set{String}()
    # read the hmm output 
    open(input) do f
        line = 1
        while !eof(f)
            x = readline(f)
            x = lstrip(x)
            if x[1] != '#'
                parts = split(x)
                rname = parts[1]
                ts = parse(Int64, parts[5])
                te = parse(Int64, parts[6])
                rs = parse(Int64, parts[9])
                re = parse(Int64, parts[10])
                ts, te = sort([ts, te])
                rs, re = sort([rs, re])

                if ((re-rs+1) >= minl) && (ts >= t_star) && (ts <= t_end) && (te >= t_star) && (te <= t_end) && (rs >= r_star) && (rs <= r_end) && (re >= r_star) && (re <= r_end)
                    println(namesf, rname)
                    push!(chosen_reads, rname)
                end
            end
            line += 1
        end
    end
    close(namesf)
    return chosen_reads
end
