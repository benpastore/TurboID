# TurboID
TurboID Mass Spectrometry Analysis 

1. Run ./parseproteinfa.sh. 

   To parse the protein fasta sequences. The script will take the longest isoform of protein for each gene. It will also output the wormbase gene ID and protein
   sequence.
  
2. Run ./iupred_analysis.sh
   
   This script finds the per base iupred score for each protein in a given fasta and outputs the results file in the format "protein".iupred.txt
   
3. Run ./process_iupred.py.
   
   This will process output from iupred_analysis.sh. Takes total iupred score for protein and divides by protein length to yeild mean disorder of protein

4. Used the turboid-analysis.R script for statistical analysis and plotting
