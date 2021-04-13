# RESCOT: Restriction Enzyme Set and Combination Optimization Tools
Calculate the rNMP coverage of restriction enzyme (RE) and generate optimized RE sets and RE combinations used for [ribose-seq](https://www.nature.com/articles/nmeth.3259) and the other rNMP capture techniques.

### Dependencies

The following packages are required by RESCOT:

- numpy
- tqdm

### Usage

#### Calculate rNMP coverage

The rNMP coverage is the ratio of rNMP captured region in genome. It can be calculated using reference genome and RE set. 

Command:

```bash
resimulation.py [res] [fasta]
```

Positional parameter:

1. __res__: A text file with RE set information. This file should have 3 columns separated by tab and NO header. The first column is RE name. The second column is the RE pattern, degenerated bases are supported. The third column is the cut position of pattern. Each row represent one RE. There is a sample file named __res.txt__ in the examples.
2. __fasta__: Reference genome sequence in __FASTA__ format

Optional parameter:

1. __-o O__: Output file base name, default = FASTA filename
2. __-l l__: Minimum suitable genome sequence length, default = 81
3. __-L L__: Maximum suitable genome sequence length, default = 481 
4. __--circular CIRCULAR [CIRCULAR ...]__: Circular chromosome/plasmids, default = chrM

Output:

rNMP coverage, RE cut sites, and rNMP missing part.

#### Optimize RE set

Generate optimized RE set using a RE candidate pool and reference genome.

Command:

```bash
optimize_RE_set.py [REpool] [fasta]
```

Positional parameter:

1. __REpool__: RE information file for all RE candidates. Should have the same format as  __res.txt__ in the examples.
2. __fasta__: Reference genome sequence in __FASTA__ format

Optional parameter:

1. __-o O__: Output file base name, default = FASTA filename
2. __-l l__: Minimum suitable genome sequence length, default = 81
3. __-L L__: Maximum suitable genome sequence length, default = 481
4. __-i I__: Number of iterations, default = 100
5. __-n N__: Minimum number of RE in a set, default = 1
6. __-N N__: Maximum number of RE in a set, default = 4
7. __--circular CIRCULAR [CIRCULAR ...]__: Circular chromosome/plasmids, default = chrM
8. __--start N,N,N......__: Index of REs in candidate pool to start optimization with, zero-started index should separated by comma, example:(0,1,3,5)

Output:
RE list, rNMP coverage, RE cut sites, and rNMP missing part for optimized RE set. And the optimization history

#### Optimize RE combination

Generate optimized RE combination using a RE candidate pool and reference genome.

Command:

```bash
optimize_RE_combination.py [REpool] [fasta]
```

Positional parameter:

1. __REpool__: RE information file for all RE candidates. Should have the same format as  __res.txt__ in the examples.
2. __fasta__: Reference genome sequence in __FASTA__ format

Optional parameter:

1. __-o O__: Output file base name, default = FASTA filename
2. __-l l__: Minimum suitable genome sequence length, default = 81
3. __-L L__: Maximum suitable genome sequence length, default = 481
5. __-n N__: Minimum number of RE in a set, default = 1
6. __-N N__: Maximum number of RE in a set, default = 4
6. __-m M__: Number of RE sets in the combination, default=2 
7. __-i I__: Number of random RE sets generated, default = 100
8. __-c C__: Number of iterations for combination, default=200000
9. __--circular CIRCULAR [CIRCULAR ...]__: Circular chromosome/plasmids, default = chrM
10. __--start N,N,N......__: Index of REs in candidate pool to start optimization with, zero-started index should separated by comma, example:(0,1,3,5)

Output:
Random RE sets and optimization history. The optimized RE combination can be found in the optimization history

### License

This software is under GNU GPL v3.0 license

### Contact

If you have any question, please contact me at [pxu64@gatech.edu](mailto:pxu64@gatech.edu).