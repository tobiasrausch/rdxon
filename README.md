[![Build Status](https://travis-ci.org/tobiasrausch/rdxon.svg?branch=master)](https://travis-ci.org/tobiasrausch/rdxon/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/rdxon/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/rdxon.svg)](https://github.com/tobiasrausch/rdxon/releases)
[![GitHub Issues](https://img.shields.io/github/issues/tobiasrausch/rdxon.svg)](https://github.com/tobiasrausch/rdxon/issues)

# FASTQ filtering for rare k-mers

Filter FASTQ files against all 1000 Genomes sequencing data using k-mers. Keep only reads with k-mers missing in 1000 Genomes.

# Installation

`git clone --recursive https://github.com/tobiasrausch/rdxon.git`

`cd rdxon/`

`make all`

# 1000 Genomes k-mer map

Download the 1000 Genomes k-mer maps here: [https://gear.embl.de/data/rdxon/](https://gear.embl.de/data/rdxon/)

# Running

To filter an input FASTQ file against the 1000 Genomes sequencing data simply run

`./src/rdxon -x kmer.x.map -y kmer.y.map -o <output.fq.gz> <input.fq.gz>`

You can also dump all rare k-mers which are absent in 1000 Genomes to a file

`./src/rdxon -x kmer.x.map -y kmer.y.map -u <kmer.gz> -o <output.fq.gz> <input.fq.gz>`

# Acknowledgement

The 1000 Genomes high-coverage data were generated at the New York Genome Center with funds provided by NHGRI Grant 3UM1HG008901-03S1. All cell lines were obtained from the Coriell Institute for Medical Research and from the NIGMS Human Genetic Cell Repository at the Coriell Institute for Medical Research. More information regarding the 1000 Genomes high-coverage data and data reuse is available [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/).
