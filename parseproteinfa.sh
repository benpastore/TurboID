#!/usr/bin/sh

# parse protein fasta into fasta without isoforms 
module load R

FA=$1
NAME=`basename $FA .protein.fa.gz`

# Make fasta with only the CDS in the header
zcat $FA | awk '{split($0,a," "); print a[1] }' > $name.protein.CDS.fa

# Make a table with Isoform / Sequence name / Wormbase name 
zcat $FA |\
    grep ">" |\
    sed -e 's/>//g' |\
    sed -e 's/wormpep=//g' |\
    sed -e 's/gene=//' |\
    awk -v OFS='\t' '{ print $1,$3 }' > $name.protein.ID

# Linearise fasta
cat $name.protein.CDS.fa |\
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |\
    sed -e 's/>//g' |\
    awk -v OFS='\t' '{print $1,$2,$3=length($2) }' > $name.protein.table

echo """
table = read.delim(\"$name.protein.table\", header = F, sep = '\t')
ID = read.delim(\"$name.protein.ID\", header = F, sep = '\t')

names(table) = c(\"transcript\", \"sequence\", \"length\")
names(ID) = c(\"transcript\",\"wormbase\")

master = merge(table, ID, by.x = \"transcript\")
col_order = c('transcript', 'wormbase', 'sequence', 'length')
master = master[ ,col_order]
write.table(master,\"$name.protein.ID.table.merge\", col.names = F, sep = '\t', row.names = F, quote = F)
""" > mergeTables.R
Rscript mergeTables.R

cat $name.protein.ID.table.merge |\
    sort -k2,2 -k4,4nr |\
    sort -k2,2 -u -s |\
    cut -f2,3,4 > $name.protein.ID.table.merge.unique

cat $name.protein.ID.table.merge.unique |\
    cut -f1,3 > $name.protein.unique.lengths

cat $name.protein.ID.table.merge.unique |\
    cut -f1,2 |\
    awk -v OFS='\t' '{ print ">"$1"\n"$2 }' | fold -w 60 > $name.protein.unique.fa

