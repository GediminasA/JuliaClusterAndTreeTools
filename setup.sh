julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.precompile()'
julia --project=. -e 'using Pkg; Pkg.build()'
