You can build a [rdxon](https://github.com/tobiasrausch/rdxon) singularity container (SIF file) using

`sudo singularity build rdxon.sif rdxon.def`

Once you have built the container you can run analysis using

`singularity exec rdxon.sif rdxon -x kmer.x.map -y kmer.y.map -o out.fq.gz in.fq.gz`
