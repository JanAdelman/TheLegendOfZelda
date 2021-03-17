#!/bin/bash

mkdir BiasAwayResults
cd BiasAwayResults
mkdir BiasAwayDinucleotide
mkdir BiasAwayMononucleotide
cd ..

# synthetic k-mer generation  
#for repetition in {1..25}
#do
#biasaway k -f drosophila_genome.fa > BiasAwayResults/BiasAwayDinucleotide/bias_away_${repetition}.fa -k 2
#done

for repetition in {1..25}
do 
biasaway k -f drosophila_genome.fa > BiasAwayResults/BiasAwayMononucleotide/bias_away_${repetition}.fa -k 1
done 

# sliding window k-mer generation 
