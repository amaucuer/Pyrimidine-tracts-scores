# Pyrimidine-tracts-scores
SPY_script_v1.R calculates a global score for the presence and strongness of pyrimidine tracts at the 3’ end of introns (as described in “U2AF65 assemblies drive sequence‐specific splice site recognition”. EMBO Reports, EMBO Press, 2019, 20 (8)). 

First, sequences files are prepared. As an example, such files can be prepared starting with the “test_SPY” file: (chr, first base to extract, last base to extract, exon coordinates, gene name, strand) (exons coordinates and gene names are optional and can be replaced by anything. This test file corresponds to the exons studied in Tari et al, Embo rep. 2019 figure S14.
The script below (bash) will extract the sequences in file “test_SPY_seq” and then create a matrix with nucleotides “test_SPY_tab “
These two files serve as input for R script “SPY_script” that calculates SPY scores.

                                                                         
bedtools getfasta -fi ~/rnastar/GRCh38.primary_assembly.genome.fa  -bed test_SPY -fo test_seq -s -tab ;
paste test_SPY test_seq | cut -f 4,5,6,8 | sed 1i"exon\tgene\tstrand\tintron_sequence"> test_SPY_seq;

awk 'NR>1 {split($4,a,""); printf $1"\t";for (i=1; i<102-length(a);i++) printf "N\t";for(i=1;i<length(a);++i) printf a[i]"\t"; print a[length(a)]}' test_SPY_seq | sed 1i"exon"> test_SPY_tab  


For any questions or concerns please contact alexandre.maucuer_at_inserm.fr
