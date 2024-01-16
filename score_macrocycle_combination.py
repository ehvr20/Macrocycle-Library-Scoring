from rdkit import Chem
from rdkit.Chem import QED
from itertools import product
import pandas as pd
import argparse

def generate_library(peptide_scaffold, building_block_list):
    """
    Generate macrocycle sequences based on peptide scaffold and building blocks.

    Args:
        peptide_scaffold (str): Peptide scaffold with '*' as a placeholder.
        building_block_list (list): List of building blocks for each '*' in the scaffold.

    Yields:
        str: A generated macrocycle sequence.
    """
    count = 0
    char_list = []

    for char in peptide_scaffold:
        if char == "*":
            char_list.append(building_block_list[count])
            count += 1
        else:
            char_list.append([char])

    for string in product(*char_list):
        yield "".join(string)

def read_file(input_file):
    """
    Read input file line by line and yield (scaffold, building blocks) tuples.

    Args:
        input_file (str): Path to the input file.

    Yields:
        tuple: (scaffold, building blocks) extracted from each line.
    """
    with open(input_file, 'r') as file:
        for line in file:
            input_values = line.split()
            yield input_values[0], input_values[1:]

def main(input_file, output_path=None):
    """
    Main function to process input file, generate macrocycles, and output QED scores to a csv file.

    Args:
        input_file (str): Path to the input file.
        output_path (str, optional): Path to the output file. Default is "./outs/macrocycle_scores.csv".
    """
    if output_path is None:
        output_path = "./outs/macrocycle_scores.csv"

    headers = ["Score", "Sequence", "MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]
    database = pd.DataFrame(columns=headers)

    for scaffold, building_block in read_file(input_file):
        for sequence in generate_library(scaffold, building_block):
            macrocycle = Chem.MolFromSequence(sequence)

            if macrocycle is not None:
                Chem.SanitizeMol(macrocycle)
            else:
                print(f"{sequence} is invalid")
                continue

            score = QED.default(macrocycle)
            #score = QED.qed(macrocycle, w=(0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95))
            #score = QED.weights_max(macrocycle)
            properties = QED.properties(macrocycle)

            # Create a new row for the DataFrame
            new_row = {
                'Score': score,
                'Sequence': sequence,
                'MW': properties[0],
                'ALOGP': properties[1],
                'HBA': properties[2],
                'HBD': properties[3],
                'PSA': properties[4],
                'ROTB': properties[5],
                'AROM': properties[6],
                'ALERTS': properties[7],
            }

            # Check if the new row is not empty before adding it
            if new_row:
                new_index = len(database)
                database.loc[new_index] = new_row

    # Sort the DataFrame by the 'Score' column in descending order
    database = database.sort_values(by="Score", ascending=False)

    print(output_path)
    database.to_csv(output_path, index=False, float_format="%.20f")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and score macrocycles from input file.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("-o", "--output_file", help="Path to the output file")
    args = parser.parse_args()
    main(args.input_file, args.output_file)
    exit()


#TODO
# RING formation
# Isomerisation
# Entropy optimize weight for macrocycles
# SMILES ENTRY FORMAT
# Unnatural Aminio Acids
    
    L and D amino acids