# KISS Search Program

## Introduction
Welcome to the KISS Search Program, a simple and efficient tool for searching reference sequences against data sequences using k-mers. KISS, which stands for K-Mer Indexing and Sequence Search, provides a fast and memory-efficient solution for identifying similarities between sequences.

## Quick Start
To compile the program and generate the executable (`KISS.out`), simply run the following command in your terminal:

```bash
make
```
## Usage
Once you have compiled the program, you can use it by running the following command:
```bash
./KISS.out -reference <reference_file> -data <data_file>
```
example
```bash
./KISS.out -reference examples/hpvcomplete.fasta -data examples/hpvE6.fasta
```
## Modification Example
You can customize the search parameters by providing additional options. Here's an example:
```bash
./KISS.out -kmer_size 32 -min_matches 3 -seed_skip 30 -cut_off 100 -reference examples/P_falciprum.fna -data examples/P_vivax.fna
```

## Options
- `-kmer_size <value>`: Set the k-mer size (default is 31).
- `-min_matches <value>`: Set the minimum number of matches required (default is 2).
- `-seed_skip <value>`: Set the seed skip value (default is 20).
- `-cut_off <value>`: Set the cutoff value (default is 50).
- `-reference <reference_file>`: Specify the reference sequence file.
- `-data <data_file>`: Specify the data sequence file.
- `-match_all <value>`: match_all indicates whether all matches (true|false), (default is true)
- `-get_sw <value>`: indicates whether Smith Waterman alignment score is required (true|false),(default is true)

##Indexing files
NOTE: the indexer is from [a link](https://github.com/andrewtrotman/KISS), this version of the indexer does have a bug which has been fixed locally, but indexing is done seperate and the the .idx files can be loaded into the examples folder of the programme so they can be called at run time (if you wish)
