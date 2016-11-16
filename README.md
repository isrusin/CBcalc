[![license](https://img.shields.io/github/license/mashape/apistatus.svg)]()

# CBCalc - Compositional Bias Calculation

## Installation
*Note: only Linux version is currently available.*

Download the last release as [tar.gz](https://github.com/isrusin/cbcalc/releases).

Unpack:
```
tar -xf cbcalc.tar.gz
```

Make C library:
```
cd cbcalc/cbclib
make dll
```

## Usage
CBCalc is a command line util. For convenient usage, cbcalc.py should have execute permission and its parent folder should be in the PATH environment variable. In the case it could be executed as:

```
cbcalc.py [-s] [-o] [METHODS] MODE INPUT
```

##### Options (list of sites, output file, methods of calculation)
`-s` is for input list of sites, one-per-line. If ommitted STDIN is used instead.

`-o` is for output file name specification, STDOUT is default.

`METHODS` is any combination of `-M`, `-P`, and `-K` flags which determines methods of compositional bais calculation (Mmax based, Pevzner's, and Karlin's correspondingly) to use and the order of the output columns. Default combination is `-MPK`.

##### Modes and arguments (input sequences)
CBCalc has two acting modes: "file" mode and "path" mode which use different specification of INPUT argument section. In "file" mode one or several fasta file names should be specified. The file names (with .fasta[.gz] extension ommited) are used as sequence names in the output.

In "path" mode the util takes single path with `{}` placeholder for sequence ID. A list on sequence IDs should be specified with `-i/--id` agrument. The IDs will be used as sequence names in the output.

##### Examples:
```
cbcalc.py -s sites.list -o output.tsv file input1.fasta input2.fasta.gz
```
```
printf "GATC\nGANNTC" | cbcalc.py -M -K file fasta_dir/*.fa.gz | head
```
```
cbcalc.py -MMMo repreated_columns.tsv path fasta_dir/{}.fa -i acs.list
```

Please, try `cbcalc.py -h`, `cbcalc.py file -h`, and `cbcalc.py path -h` for some additional details on CBCalc usage.

## Output format
Output is a tab-separated table with single header line:
```
ID <TAB> Site <TAB> No <TAB> (Xe <TAB> Xr <TAB>) Total
```

__ID__ is a sequence name, __Site__ is a target word (or pattern), __No__ is a number of word occurrences in the sequence, __Xe__ is the expected number of word occurrences estimated with method X, __Xr__ is the ratio of No to Xe, __Total__ is very like sequence length but corrected by word length (the sequence is treated as a linear one). Number and order of Xe and Xr columns are determined by the METHODS section of input arguments.

## Web-interface
For single sequence requests you could use [web-interface](http://mouse.belozersky.msu.ru/tools/cbcalc).

## Requirements
* Python 2.7
* make
* gcc
