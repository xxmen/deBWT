#!/bin/bash
d_file=$1
declare -i kmer = $2
declare -i kmer2 = $kmer+2
start = `date +%s`

jellyfish count -m $kmer2 -s 5M -t 10 $d_file
jellyfish dump mer_counts.jf >| mer_counts_dumps.fa

end = `date +%s`
usetime = $((end-start))
echo "jellyfish use: $usetime 's"

(make)

start = $end
./deBWT -k $kmer -d $d_file -m mer_counts_dumps.fa

end=`date +%s`
runtime=$((end-start))
echo "build BWT use : $runtime 's"