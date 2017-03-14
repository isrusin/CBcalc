[![license](https://img.shields.io/github/license/mashape/apistatus.svg)]()

# CBcalc - Compositional Bias calculation

## Installation

### Linux
Download the latest release
[here](https://github.com/isrusin/CBcalc/releases).

Unpack:
```
tar -xf CBcalc.tar.gz
```

Enter CBcalc folder:
```
cd CBcalc/
```

Build cbclib:
```
python setup.py build_ext -i
```

### Other platforms
Only the Linux example is currently available. However, CBcalc seems to
contain no platform-specific code. The only problem is to build `cbclib`
which contains a C-extension module depending on `Python.h` and `zlib`.

## Usage
CBcalc is a command line util. For convenient usage, cbcalc.py should
have execute permission and its parent folder should be in the PATH
environment variable. In the case it could be executed as:

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

#### Options (list of sites, output file, methods of calculation):
`-s` is for whitespace-delimited list of sites. If ommitted STDIN is used
instead.

`-o` is for output file name, STDOUT is default.

A combination of `-B`, `-M`, `-P`, and `-K` flags determines which methods
of compositional bias calculation to use and the order of the output
columns.

The methods are:

`B`, the Bernoulli model-based method;  
`M`, the method based on maximum order Markov chain model (Mmax);  
`P`, the extended Mmax-based method suggested by Pevzner and co-authors;  
`K`, the non-Markovian method described by Burge and co-authors.

Default combination is `-MPK`.

#### Arguments (input sequences):
CBcalc can accept Fasta files in two ways. Without `-i`/`--id` option, it
treats all positional arguments as names (paths) of fasta files (may be
gzipped). The file names (without parent folders and `.fasta[.gz]`
extension) will be used as sequence names in the output table.

Another way is to provide a whitespace-delimited list of sequence IDs as a
text file through `-i`/`--id` option. In the case, the only positional
argument is a batch for path to Fasta files. Use `{}` placeholder for
sequence ID in the batch. The IDs will be used in the output table as
sequence names.

#### Examples:
```
cbcalc.py -s sites.list -o output.tsv input1.fasta input2.fasta.gz
```
Calculate compositional biases for all sites from the `sites.list` in the
fasta files `input1.fasta` and `input2.fasta` with the default set of
methods (`M`, `P`, `K`) and write resulted table to the `output.tsv`.

```
printf "GATC GANNTC" | cbcalc.py -M -K fasta_dir/*.fa.gz
```
Calculate compositional biases for sites `GATC` and `GANNTC` in all files
from the `fasta_dir` whose name ends with `.fa.gz` with `M` and `K` methods
and write the resulted table to STDOUT.

```
cbcalc.py fasta_dir/{}.fa -i seq_id.list
```
Calculate compositional biases with the default set of methods for each
file from the `fasta_dir` named `X.fa` where `X` is a sequence ID from
the `acs.list` and write the resulted table to STDOUT, a list of sites
will be obtained from STDIN.

Please, try `cbcalc.py --help` for some additional details on CBcalc usage.

## Output format
The output is a tab-separated table with the following columns:

`ID`, the sequence name;  
`Site`, the target word or pattern;  
`Observed`, the number of word occurrences in the sequence;  
`Expected (X)`, the expected number of words estimated with the method X;  
`Ratio (X)`, the observed/expected ratio (compositional bias);  
`Total`, a value similar to sequence length but corrected by word length.

The number and the order of `Expected` and `Ratio` columns are
determined by the number and the order of the method options.

## Web-interface
[Web-interface](http://mouse.belozersky.msu.ru/tools/cbcalc) could be used
for single-sequence requests.

## Requirements
* Python 2.7
* python-dev package
