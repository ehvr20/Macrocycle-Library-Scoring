from rdkit import Chem
from rdkit.Chem import AllChem, QED, Draw, MolStandardize
import pandas as pd

from Formatting import *
from Generation import *

def constructMolecule(sequence):
    """
    Constructs a molecule from the given sequence.

    Args:
        sequence (str): A peptide sequence.

    Returns:
        RDKit Mol object: The molecule constructed from the sequence.
    """
    return Chem.MolFromSequence(sequence)

def main():
    """
    Main function for processing peptide sequences and generating a DataFrame.
    """
    # Define the headers for the DataFrame
    headers = ["Sequence", "Score", "MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]

    # Initialize an empty DataFrame with the defined headers
    storage = pd.DataFrame(columns=headers)

    # Define peptide scaffold and building blocks
    peptideScaffold = ["AA*DA*WY", "AA*YW*WY"]
    buildingBlockList = [["ARN", "KMF"], ["WYF", "KAF"]]
    
    # TODO: Import and format both input values

    for scaffold, building_blocks in zip(peptideScaffold, buildingBlockList):
        for sequence in generate_library(peptideScaffold=scaffold, buildingBlockList=building_blocks):
            macrocycle = constructMolecule(sequence)

            highest_score = 0
            highest_properties = None

            # Enumerate stereoisomers and calculate scores
            for isomer in AllChem.EnumerateStereoisomers(macrocycle):
                score = QED.default(isomer)
                if score > highest_score:
                    highest_properties = QED.properties(isomer)
                    highest_score = score

            print(sequence, highest_score, highest_properties)

            # Create a new row for the DataFrame
            new_row = {
                'Sequence': sequence,
                'Score': highest_score,
                'MW': highest_properties[0],
                'ALOGP': highest_properties[1],
                'HBA': highest_properties[2],
                'HBD': highest_properties[3],
                'PSA': highest_properties[4],
                'ROTB': highest_properties[5],
                'AROM': highest_properties[6],
                'ALERTS': highest_properties[7],
            }

            # Check if the new row is not empty before adding it
            if new_row:
                new_index = len(storage)
                storage.loc[new_index] = new_row

    # Sort the DataFrame by the 'Score' column in descending order
    storage = storage.sort_values(by="Score", ascending=False)

    # Save the DataFrame to an Excel file
    storage.to_excel("output.xlsx", index=False, float_format="%.20f")

if __name__ == "__main__":
    main()
