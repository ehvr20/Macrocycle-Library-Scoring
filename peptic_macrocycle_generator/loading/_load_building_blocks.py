import os
from rdkit import Chem

def load_building_blocks():
    mols = []
    for file in os.listdir('./data/building_blocks'):
        file = os.path.join('./data/building_blocks/', file)
        print("Loading building blocks from", file)
        with open(file) as building_blocks:
            for building_block in building_blocks:
                molecule = Chem.MolFromSmiles(building_block)
                mols.append(molecule)
    return mols
