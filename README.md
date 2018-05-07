This project is a work-in-progress attempt at reproducing the phylogenetic tree analysis described [here][original paper]. The idea is to produce a graph of phylogenetic relationships between different bacteria using an alignment-free analysis of genomic data.

A corresponding visualization is available [here][viz].

## Files in this repository

- `d2.py` contains Python implementations of the [d2 and d2s statistics][d2].
- `analysis.py` is a script for computing d2 distances between all possible pairs of genomic sequences for a given directory of genomic data.
- `ids.txt` contains the NCBI RefSeq identifiers of the 144 bacteria sequences used in [the original paper][original paper] this repository attempts to reproduce.
- `phylum_data.py` contains some metadata concerning these bacteria.

## Instructions

1. Obtain FASTA (`.fna`) files for a collection of genomic sequences. To use the same sequences analyzed in [this paper][original], you can use the `ids.txt` file in the root of this repository along with the [NCBI Batch Entez tool][ncbi batch entrez]. All considered, this downloads ~400MB of FASTA data.

1. Make sure the [pyfasta][pyfasta] Python library is available in your Python path.

1. Specify the directory of your genomic data with the `DATA_DIR` environment variable, and then run the `analysis.py` file. This will print the computed d2 data in JSON format to stdout. To save this file, you might do something like the following:

   ```
   DATA_DIR=../my_data python analysis.py > out.json
   ```

1. The JSON file can be used with the visualization in [this notebook][viz].

[ncbi batch entrez]: https://www.ncbi.nlm.nih.gov/sites/batchentrez
[pyfasta]: https://github.com/brentp/pyfasta/
[original viz]: http://bioinformatics.org.au/tools/AFnetwork/
[original paper]: https://f1000research.com/articles/5-2789/v2
[viz]: https://beta.observablehq.com/@chnn/producing-phylogenetic-graphs-with-the-tex-d_2-statistic
[d2]: https://www.researchgate.net/publication/40679234_Alignment-Free_Sequence_Comparison_I_Statistics_and_Power
