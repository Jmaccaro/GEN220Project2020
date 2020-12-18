#!/bin/bash

#To make a fungal database of just the Onygenales proteomes. It is from Jason and contains one outgroup (A. fumigatus)
ln -s /bigdata/stajichlab/shared/projects/Onygenales/SFD/Phylogeny/pep
cp pep/allseq .

#To build an hmm for my ortholog of interest (D-arabinitol dehydrogenase) and search that against all the other Onygenales
module load hmmer/3
hmmbuild ArDh.faa > ArDH.hmm
hmmsearch --domtblout ArDh.domtbl ArDH.hmm allseq > ArDh.hmmsearch

#To generate a list of the top hit names and exclude repeats 
cat ArDh.domtbl | head -55 | awk '!seen[$1]++' | grep -v '^#'| awk '{print$1}' > ArdH.top.hits.names.txt

# To reformat my database 
cat allseq | tr -d '\n' > allseq.1line.fa 
cat allseq.1line.fa | sed 's/>/\'$'\n>/g' > allseq.lineb4.faa
cat allseq.lineb4.faa | sed -E 's/>[A-Z]+\|[A-Z]+\_[0-9]+/&\n/g' > allseq.1line.faa

#To pull sequences for my top hit names into one file 
mkdir ArDhseqs
for names in $(cat ArdH.top.hits.names.txt)
do 
grep -A 1 $names allseq.1line.faa > ArDhseqs/$names.fasta 
done
cat ArDhseqs/*.fasta > allseqs.top.ArDh.faa

#Then I manually cleaned up the lines that had annotation data in the sequences
#To remove repeats that arrose for some reason after pulling sequences. Then added my SCO file. 
grep -A 1 -f ArdH.top.hits.names.txt ArDh.corrected.faa | grep -v '^--' > clean.uniq.ArDhseqs.faa
cat ArDH.noAAP.faa clean.uniq.ArDhseqs.faa > ArDh.clean.withAsco.faa

#To align seqs and build a tree
module load hmmer/3
module load db-pfam
hmmalign ArDH.hmm ArDh.clean.withAsco.faa > clean.ArDh.stalkholm
esl-reformat clustal clean.ArDh.stalkholm > clean.ArDh.stalkholm.clustal
iqtree2 -s clean.ArDh.stalkholm.clustal -nt AUTO
