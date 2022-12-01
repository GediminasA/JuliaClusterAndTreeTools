
#Splits  output of BLAST parsing by readblastXML from BioTools to groups matching one query  (queryname - key of a output dictionary) 
function split_by_query(blastres::Array{BLASTResult,1})
    out = Dict{String,Array{BLASTResult,1}}()
    for r in blastres
        if !(r.queryname in keys(out))
            out[r.queryname] = Array{BLASTResult,1}()
        end 
        push!(out[r.queryname], r)
    end
    return(out) 
end

#check the length of a new interval that are not covered by the one in the array
function unique_length_fraction(matches::Array{Interval{:closed,:closed,Int64},1}, interval::Interval{:closed,:closed,Int64})
    difs = Array{Interval{:closed,:closed,Int64},1}()
    unique_length = interval.right - interval.left + 1
    length_all = deepcopy(unique_length)
    for i in matches
        common = intersect(interval, i)
        if !isempty(common)
            unique_length -= common.right - common.left + 1 
        end
    end
    return(unique_length/length_all)
end

#add new interval to interval sets
function add_new_interval(matches::Array{Interval{:closed,:closed,Int64},1}, interval::Interval{:closed,:closed,Int64})
    intersecting = Array{Interval{:closed,:closed,Int64},1}([interval])
    has_intersecting = false
    
    for i in matches
        if !isempty(intersect(interval, i))
            push!(intersecting, i)
            has_intersecting = true
        end
    end
    #get union of intersecring intervals
    if has_intersecting
        expanded = deepcopy(interval)
        for i in intersecting
            expanded = union(expanded,i)
        end
    end
    out = Array{Interval{:closed,:closed,Int64},1}()
    #collect output
    if has_intersecting
        for i in matches
            if isempty(intersect(i, expanded))
                push!(out, i)
            end
        end
        push!(out, expanded)
    else
        expanded = deepcopy(matches)
        push!(expanded, interval)
        out = expanded
    end
    return(sort!(out, by = x-> x.left))
end

function covered_length(intervals::Array{Interval{:closed,:closed,Int64},1})
    out = 0
    for i in intervals
        out += i.right - i.left +1
    end
    return(out)
end

#collects set of hits that extends the coverage. Only hits that overlaps witha certain length limit with higher scoring hits are considered
function get_coverage_maximising_hits(queryid::String, data::Dict{String,Array{BLASTResult,1}}, max_overlap_frac = 0.5)
    hits = data[queryid]
    hits = sort(hits, by = x -> (x.queryto - x.queryfrom +1), rev = true )
    coverge_expanding_hits = Array{BLASTResult,1}([hits[1]])
    coverge_expanding_hits_intervals = [Interval(hits[1].queryfrom,hits[1].queryto)]
    current_length = covered_length([Interval(hits[1].queryfrom, hits[1].queryto)])
    for h in hits[2:end]
        query_interval = Interval(h.queryfrom,h.queryto)
        intervals_with = add_new_interval(coverge_expanding_hits_intervals, query_interval)
        uniqfrac = unique_length_fraction(coverge_expanding_hits_intervals, query_interval)
        new_length = covered_length(intervals_with)
        if new_length > current_length && uniqfrac < max_overlap_frac
            current_length = new_length
            push!(coverge_expanding_hits, h)
            coverge_expanding_hits_intervals = intervals_with
        end
    end
    return(coverge_expanding_hits)
end