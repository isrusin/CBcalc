[![license](https://img.shields.io/github/license/mashape/apistatus.svg)]()

# CBcalc - Compositional Bias calculation

## Installation (Linux example)

Download the latest release [here](https://github.com/isrusin/CBcalc/releases).

Unpack:
```
tar -xf cbcalc.tar.gz
```

Install cbclib:
```
python setup.py install --install-platlib=.
```

## Usage
CBcalc is a command line util. For convenient usage, cbcalc.py should have execute permission and its parent folder should be in the PATH environment variable. In the case it could be executed as:

```
cbcalc.py [-s] [-o] [-BMPK] FILE [FILE ...]
```
or
```
cbcalc.py [-s] [-o] [-BMPK] -i PATH
```

Otherwise, you have to use full form:

```
python [path_to_cbcalc]/cbcalc.py ...
```

*Note: CBcalc uses python 2.7*
##### Options (list of sites, output file, methods of calculation)
`-s` is for input list of sites, one-per-line. If ommitted STDIN is used instead.

`-o` is for output file name, STDOUT is default.

A combination of `-B`, `-M`, `-P`, and `-K` flags determines which methods of compositional bias calculation to use and the order of the output columns.

The methods are:
B -- the Bernoulli model-based method;
M -- the maximum order Markov chain model (Mmax)-based method;
P -- the extended Mmax-based method suggested by Pevzner and co-authors;
K -- the non-Markovian method described by Karlin and co-authors.

Default combination is `-MPK`.

##### Arguments (input sequences)
CBcalc has two acting modes. Without `-i`/`--id` option, it treat all positional arguments as names (paths) of fasta files (may be gzipped). The file names (without parent folders and .fasta[.gz] extension) will be used as sequence names in the output.

You can provide a list of sequence IDs as a text file by `-i`/`--id` option. In the case, the only positional argument is a batch for fasta files path. Use `{}` placeholder for sequence ID in the batch. The IDs will be used in the output table as sequence names.

##### Examples:
```
cbcalc.py -s sites.list -o output.tsv input1.fasta input2.fasta.gz
```
Calculate compositional biases for all sites from `sites.list` in fasta files `input1.fasta` and `input2.fasta` with the default set of methods (`M`, `P`, `K`) and write resulted table to `output.tsv`.

```
printf "GATC\nGANNTC" | cbcalc.py -M -K fasta_dir/*.fa.gz
```
Calculate compositional biases for sites `GATC` and `GANNTC` in all files from the `fasta_dir` whose name ends with .fa.gz with `M` and `K` methods and write the resulted table to STDOUT.

```
cbcalc.py -MMMo repeated_columns.tsv fasta_dir/{}.fa -i seq_id.list
```
Calculate compositional biases with `M` method for each file in `fasta\_dir` named X.fa where X is a sequence ID from `acs.list` and write the resulted table with triple repeated "Expected (M)" and "Ratio (M)" columns to `repeated_columns.tsv`, list of sites will be obtained from STDIN.

Please, try `cbcalc.py --help` for some additional details on CBcalc usage.

## Output format
Output is a tab-separated table with a single header line:
```
ID <TAB> Site <TAB> Observed <TAB> (Expected (X) <TAB> Ratio (X) <TAB>) Total
```

__ID__ is a sequence name, __Site__ is a target word (or pattern), __Observed__ is a number of word occurrences in the sequence, __Expected (X)__ is the expected number of word occurrences estimated with method X, __Ratio (X)__ is the ratio of the __Observed__ to the __Expected (X)__, __Total__ is a value very likes sequence length but corrected by word length (the sequence is treated as a linear one). The number and the order of __Expected__ and __Ratio__ columns are determined by the number and the order of correcponding options of the command line.

## Web-interface
[Web-interface](http://mouse.belozersky.msu.ru/tools/cbcalc) could be used for single-sequence requests.

## Requirements
* Python 2.7
* python-dev package (linux)
