from rdkit import Chem

def load_macrocycles():
    macrocycle_scaffolds = []

    with open("./data/macrocycles/macrocycles.tsv") as macrocycle_file:

        for macrocycle_scaffold in macrocycle_file:

            macrocycle_scaffold = macrocycle_scaffold.strip().split("\t")

            macrocycle_scaffold = ['*' if component == '*' else [Chem.MolFromSmiles(component)] for component in macrocycle_scaffold]
            
            macrocycle_scaffolds.append(macrocycle_scaffold)

    return macrocycle_scaffolds