[![Build Status](https://travis-ci.org/tobiasrausch/rdxon.svg?branch=master)](https://travis-ci.org/tobiasrausch/rdxon/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/rdxon/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/rdxon.svg)](https://github.com/tobiasrausch/rdxon/releases)
[![GitHub Issues](https://img.shields.io/github/issues/tobiasrausch/rdxon.svg)](https://github.com/tobiasrausch/rdxon/issues)

# FASTQ filtering for rare k-mers

# Installation

`git clone --recursive https://github.com/tobiasrausch/rdxon.git`

`cd rdxon/`

`make all`

# Running

`./src/rdxon -x kmer.x.map -y kmer.y.map -o <out.fq> <in.fq>`
