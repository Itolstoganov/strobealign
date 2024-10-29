## Running the uniqueness benchmark
1. Install the required software
```
  conda env create --file environment.yml
  conda activate uniqueness_benchmark
```
2. Download the CHM13 v.2.0 reference genome
3. Run the benchmark script
```
  ./uniqueness_benchmarks.py -r <reference genome> -o <output directory>
```
4. The results table can be found in the `output/stats.tex` file
