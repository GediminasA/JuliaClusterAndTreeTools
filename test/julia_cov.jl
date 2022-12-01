#!/usr/bin/env julia
using LightXML

#=
    Code cov parser
=#

function main(args)

    print_lines = Dict()
    util_lines = Dict()
    uncov_lines = Dict()
    pwd = args[3]
    files_list = args[4:end]
    for agr in files_list
        ct = 1
        open(agr) do file
            print_lines[agr] = []
            util_lines[agr] = []
            uncov_lines[agr] = []
            for ln in eachline(file)
                if occursin("println", ln)
                    push!(print_lines[agr], string(ct))
                elseif occursin(r"(^\s*end\s*$|function|else|using|#|^\s*$|^\")", ln)
                    push!(util_lines[agr], string(ct))
                else
                    push!(uncov_lines[agr], string(ct))
                end
                ct = ct + 1
            end
        end
    end

    xdoc = XMLDocument()
    xroot = create_root(xdoc, "coverage")
    set_attributes(xroot, Dict("version"=>"1.0"))
    packs = new_child(xroot, "packages")
    pack = new_child(packs, "package")
    classes = new_child(pack, "classes")
    open(args[1]) do file
        funk_hit = false
        full_file_name = ""
        lines_c = ""
        for ln in eachline(file)
            line_tmp = split(ln, ":")
            if line_tmp[1] in ["end_of_record", "LH", "LF"]
                if funk_hit
                    for line_uncov in uncov_lines[full_file_name]
                        line = new_child(lines_c, "line")
                        set_attributes(line, Dict("hits"=>0,
                                                  "number"=>line_uncov))
                    end
                end
                funk_hit = false
            elseif line_tmp[2] in files_list
                funk_hit = true
                clas = new_child(classes, "class")
                full_file_name = line_tmp[2]
                file_name = replace(full_file_name, pwd * "/" => "")
                file_name_no_src = replace(file_name, "src/" => "")
                set_attributes(clas, Dict("filename"=>file_name,
                                          "name"=>file_name_no_src))
                new_child(clas, "methods")
                lines_c = new_child(clas, "lines")
            elseif funk_hit
                line_hits = split(line_tmp[2], ",")
                if line_hits[1] in print_lines[full_file_name]
                    hit = 1
                else
                    hit = line_hits[2]
                end
                filter!(x->x!=line_hits[1], uncov_lines[full_file_name])
                line = new_child(lines_c, "line")
                set_attributes(line, Dict("hits"=>hit,
                                          "number"=>line_hits[1]))
            end
        end
    end
    save_file(xdoc, args[2])

end

main(ARGS)
