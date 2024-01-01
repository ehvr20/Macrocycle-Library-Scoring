from rdkit import Chem
from rdkit.Chem import QED
from itertools import product
import pandas as pd
import argparse

def generate_library(peptideScaffold, buildingBlockList):
    # Returns an iterable of each possible macrocyle combination parsed as a rdkit object
    count = 0
    charList = []

    for char in peptideScaffold:
        if char == "*":
            charList.append(buildingBlockList[count])
            count += 1
        else:
            charList.append([char])
    
    for string in product(*charList):
        print(string)
        yield "".join(string)

def read_file(input_file):

    protein_scaffolds, building_blocks = [], []

    with open(input_file, 'r') as file:
        for line in file:
            input_values = line.split()
            protein_scaffolds.append(input_values[0])
            building_blocks.append(input_values[1:])
            
    return protein_scaffolds, building_blocks

def main(input_file, output_path):

    if output_path is None:
        output_path="./outs/macrocycle_scores.csv"

    protein_scaffolds, building_blocks = read_file(input_file)    

    # Define the headers for the DataFrame
    headers = ["Score", "Sequence", "MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]

    # Initialize an empty DataFrame with the defined headers
    database = pd.DataFrame(columns=headers)

    for scaffold, building_block in zip(protein_scaffolds, building_blocks):
        for sequence in generate_library(scaffold, building_block):

            macrocycle = Chem.MolFromSequence(sequence)

            if macrocycle is not None: 
                Chem.SanitizeMol(macrocycle)
            else: 
                print(f"{sequence} is invalid")
                continue

            score = QED.default(macrocycle)
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

    # Save the DataFrame to an Excel file

    print(output_path)
    database.to_csv(output_path, index=False, float_format="%.20f")

if __name__ == "__main__":

    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Example Argparse Script")

    # Define command-line arguments
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("-0,","--output_file", help="Path to the output file")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args.input_file, args.output_file)