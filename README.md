# syncmer
Validation of sycnmers compared to minimizers

## Building from source on Linux
There is a Makefile. It will build the binary, but it does not have  header file dependencies so be careful if you edit headers. You need the gcc compiler. There are no dependencies on third-party libraries. 
The Makefile uses ccache, if you don't have ccache installed then edit the Makefile to delete it.

## Building from source on Windows
There is a Microsoft Visual Studio project file for VS 2017 (v15.4). This will probably work if imported into newer VS versions, I didn't test this.

## Calculate submer metrics on a random string

`syncmer -newbench output.tsv [options]`

Metrics are written to `output.tsv` in tab-separarated text format. Options:

`-k nn or -klo nn -khi nn` Value(s) of k to test (k-mer length).

`-s nn or -slo nn -klo nn` Value(s) of s (syncmer substring length).

`-w nn or -wlo nn -whi nn` Values(s) of w (minimizer window length).

`-d nn` Integer down-sampling parameter.

`-pctid nn or -pctidlo nn -pctidhi nn -pctidinc nn` Mutation rate(s).

`-seqlength nn`  Random sequence length. Default 1000000.

`-pcitidout pctid.tsv` Writes pctid metrics to tabbed file.

`spacingout spacing.tsv` Writes spacing frequencies to tabbed file.

If the lo & hi variants of an option are specified, then all values from lo to hi are tested.

### Examples

`syncmer -newbench results.tsv -k 15 -m 10 -pctidlo 75 -pctidhi 95 -pcitidinc 5`   

`syncmer -newbench results.tsv -k 15 -slo 2 -shi 13 -pctidlo 80 -pctidhi 90 -pcitidinc 10`   

## Whole-genome alignment

`sycnmer -searchsix genome1.fasta -db genome2.fasta -sts ST -k K -w W -t T -log results.txt`   

ST is the submer type: Minimizer1 for minimizers or Syncmer5 for syncmers. 
K is the kmer length, W is the minimizer window length, T is the syncmer substring length (same as -s parameter for -newbench command).

