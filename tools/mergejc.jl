using JuliaClusterAndTreeTools

if abspath(PROGRAM_FILE) == @__FILE__
    mergejc(ARGS[1], ARGS[2])
end
