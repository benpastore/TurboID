#!/usr/bin/sh

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --mem=32gb
#SBATCH --account=PCON0160

export PATH=$PATH:/fs/ess/PCON0160/bin/iupred2a

data=/fs/ess/PCON0160/projects/TurboID/data/iupred
FA=$data/c_elegans.WS230.protein.unique.fa
OUTPUT=$data/iupred_data

if [ ! -d "$output" ]
then
mkdir $output
fi

cd $data
cat $f | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | sed -e 's/>//g' > protein.table

while read gene seq
do 

echo -e ">$gene\n$seq" > tmp

iupred2a.py tmp long > $output/$gene.iupred.txt

done < protein.table