from rdkit import Chem
from rdkit.Chem import AllChem, QED, Draw, MolStandardize
import pandas as pd

from Formatting import *
from Generation import *

def construct_mol(sequence):
    return Chem.MolFromSequence(sequence)

def main():
    
    headers = ["Sequence", "Score", "MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]

    storage = pd.DataFrame(columns=headers)

    peptideScaffold = "AA*DA*WY"
    buildingBlockList = ["ARN","KMF"]
    #TODO import and format both input values

    
    for sequence in generate_library(peptideScaffold=peptideScaffold,buildingBlockList=buildingBlockList):

        macrocycle = construct_mol(sequence)

        highest_score = 0
        highest_properties = None

        # Now, you can proceed with isomer enumeration
        for isomer in AllChem.EnumerateStereoisomers(macrocycle):
            score = QED.default(isomer)
            if score > highest_score:
                highest_properties = QED.properties(isomer)
                highest_score = score
        
        print(sequence, highest_score, highest_properties)

        new_row = {
            'Sequence': sequence, 
            'Score': highest_score,
            'MW' : highest_properties[0],
            'ALOGP' : highest_properties[1],
            'HBA' : highest_properties[2],
            'HBD' : highest_properties[3],
            'PSA' : highest_properties[4],
            'ROTB' : highest_properties[5],
            'AROM' : highest_properties[6],
            'ALERTS' : highest_properties[7],
            }
        
        
        storage = storage._append(new_row, ignore_index=True)

    storage = storage.sort_values(by="Score", ascending=False)
    
    storage.to_excel("output.xlsx", index=False, float_format="%.20f")

if __name__ == "__main__":
    main()