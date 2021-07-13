#!/users/PAS1473/benpasto1/.conda/envs/rnaseq_basic/bin/python

import sys 
import os
import pandas as pd
from tqdm import tqdm

data='/users/PAS1473/benpasto1/PCON0160/projects/TurboID/data/iupred/iupred_data'

files = os.listdir(DATA)

n_files = len(files)
print(n_files)

disorders = []
proteins = []

os.chdir(data)

for file in tqdm(files) : 
    scores = []
    with open('{}'.format(file), 'r') as f : 
        for line in f : 
            if not line.startswith("#") :
               
                info = line.strip().split() 
                score = float(info[2])
                position = int(info[0])
                scores.append(score)

        total_score = sum(scores)

        disorder = total_score / position
        protein = file.replace(".iupred.txt", "")

        proteins.append(protein)
        disorders.append(disorder)
    f.close() 

res = pd.DataFrame({'protein':proteins,'disorder':disorders})
res.to_csv("/users/PAS1473/benpasto1/PCON0160/projects/TurboID/data/iupred/average_disorder_proteins.txt", sep = '\t', header = False, index = False)