

# fastjac

A simple program for getting k-mer similarity metrics for two fastq/fasta files, written in Rust.

## Description

This command line program takes as inputs two fastq/fasta file and a k-mer size and outputs the following [k-mer](https://en.wikipedia.org/wiki/K-mer) similarity metrics:

- [Jaccard distance](https://en.wikipedia.org/wiki/Jaccard_index) - the ratio of intersection over union for the two k-mer sets. 
<img src="https://render.githubusercontent.com/render/math?math=\Large \frac{ |A \cap B| }{ |A \cup B| }">




- [Sørensen-Dice coefficient](https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient) -  twice the number of elements common to both sets divided by the sum of the number of elements in each set. 
<img src="https://render.githubusercontent.com/render/math?math=\Large \frac{2 |A \cap B| }{ |A| %2b |B| }">

- Containment - fraction of the k-mers in the query which are found in the reference.
<img src="https://render.githubusercontent.com/render/math?math=\Large \frac{|A \cap B| }{ |A| }">

If the cardinality of the k-mers in **A** is smaller than in **B** then this is the same as the [overlap coefficient](https://en.wikipedia.org/wiki/Overlap_coefficient).

## Install

I provide precompiled binaries for linux only [here](https://github.com/angelovangel/fastkmers/releases/download/v0.1.0/fastjac), but it is simple to compile and run:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/angelovangel/fastjac.git

cd fastjac
cargo build --release

```

The executable file `fastjac` is now under `./target/release/`

## Usage

The main options are `-k` k-mer size to use, `-q` path to query sequence file and `-r` path to reference sequence file. Both `fasta` and `fastq` (also compressed) are supported. Try `fastjac -h` for all options.

The output is a tab-separated table with
the lengths of the k-mer sets, intersection, union, Jaccard distance, Sørensen-Dice coefficient, and containment (fraction of k-mers of query found in reference).

