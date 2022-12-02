1. Install and setup julia 1.8+, using comands like below :
```bash
    wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.3-linux-x86_64.tar.gz
    tar -xvzf julia-1.8.3-linux-x86_64.tar.gz
    sudo cp -r julia-1.8.3 /opt/
    sudo ln -s /opt/julia-1.8.3/bin/julia /usr/local/bin/julia
    rm julia-1.8.3-linux-x86_64.tar.gz
```
2. Setup required julia environment :
```bash
   cd scripts/julia_modules/st16SseqJuliaTools
   ./setup.sh 
```

3. Run the small scale analysis and a set of tests:
The workflow will download small versions of several databases. In the example below 12 threads are allowed to use.
```bash
    cd testing
    ./run_tests.py  --threads 12
```
