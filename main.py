from Formatting import *
from Generation import *
from Statistics import *

def main():

    peptideScaffold = "AI*DH*XL"
    buildingBlockList = ["ADHK","ZJK"]
    #TODO import and format both input values
    
    for peptide in generate_library(peptideScaffold=peptideScaffold,buildingBlockList=buildingBlockList):
        print(peptide)
        #TODO Run statistic analysis on each peptide and add it to a database with a unique identifier

if __name__ == "__main__":
    main()