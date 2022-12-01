
"A wrapper around taxonkit to get matching data dir"
function taxid_to_species(taxids::Array{String,1}, db_dir::String)
    db = db_dir 
    ids = taxids
    outdata=Array{String,1}()

    bindir = Conda.bin_dir(:vsearch)

    tmpdir = mktempdir("/dev/shm"; prefix="taxid_", cleanup=true)
    in_taxid_lis = "$tmpdir/taxidlist.txt"
    lineage_out = "$tmpdir/lineage.txt"
    reformat_out = "$tmpdir/reformat.txt"
    open(in_taxid_lis,"w") do idsf 
        for i in ids
            println(idsf, i)
        end 
    end

    lineage = run(`$bindir/taxonkit lineage --data-dir $db -o $lineage_out  $in_taxid_lis`)
    reformat = read(`$bindir/taxonkit reformat --miss-rank-repl "NA" --data-dir $db $lineage_out`, String)
    for l in split(reformat, "\n")
        if l != "" 
            species = split(l, ";")[end]
            push!(outdata, species)
        end 
    end 
    return(outdata)
end

"A wrapper around taxonkit to get matching data dir"
function taxid_to_genus(taxids::Array{String,1}, db_dir::String)
    db = db_dir 
    ids = taxids
    outdata=Array{String,1}()

    bindir = Conda.bin_dir(:vsearch)

    tmpdir = mktempdir("/dev/shm"; prefix="taxid_", cleanup=true)
    in_taxid_lis = "$tmpdir/taxidlist.txt"
    lineage_out = "$tmpdir/lineage.txt"
    reformat_out = "$tmpdir/reformat.txt"
    open(in_taxid_lis,"w") do idsf 
        for i in ids
            println(idsf,i )
        end 
    end

    lineage = run(`$bindir/taxonkit lineage --data-dir $db -o $lineage_out  $in_taxid_lis`)
    reformat = read(`$bindir/taxonkit reformat --miss-rank-repl "NA" --data-dir $db $lineage_out`, String)
    for l in split(reformat, "\n")
        if l != "" 
            genus = split(l, ";")[end-1]
            push!(outdata, genus)
        end 
    end 
    return(outdata)
end