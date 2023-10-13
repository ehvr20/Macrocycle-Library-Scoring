from Formatting import *
from Generation import *
from macrocycle import *

def main():

    peptideScaffold = "AA*DA*WY"
    buildingBlockList = ["ARN","KMF"]
    #TODO import and format both input values
    
    for sequence in generate_library(peptideScaffold=peptideScaffold,buildingBlockList=buildingBlockList):
        with macrocycle(sequence) as peptide:
            #print(peptide.vars["peptide_sequence"])
            #print(peptide.vars["smiles"])
            #print(peptide.vars["molecule"])
            peptide.calculate_score()
            print(sequence, peptide.vars["score"])


if __name__ == "__main__":
    main()