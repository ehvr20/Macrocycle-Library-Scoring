from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import peptic_macrocycle_generator as pmg


def scaffold_from_macrocycle(macrocycle):

    residues = pmg.chem.hydrolysis(macrocycle)

    # Draw the 2D depiction and display in the terminal
    print("Using the indexes for each residue as reference, please note down the residues that are fixed")

    # Compute 2D coordinates for the molecules
    for mol in residues:
        AllChem.Compute2DCoords(mol)

    # Create a grid image with legends and apply options
    img = Draw.MolsToGridImage(residues, molsPerRow=2, subImgSize=(400, 400),
                               legends=[f"Residue {i + 1}" for i in range(len(residues))])

    # Display the grid image
    img.show()

    fixed_residues = input("Which residues are fixed? (i.e. 1,4 )\n")

    fixed_indexes = fixed_residues.split(",")
    fixed_indexes = [int(index) for index in fixed_indexes]

    scaffold = [Chem.MolToSmiles(
        residue) if i+1 in fixed_indexes else "*" for i, residue in enumerate(residues)]
    scaffold = "\t".join(scaffold)

    print("Scaffold: ", scaffold)

    with open('./data/macrocycles/macrocycles.tsv', 'a') as file:
        file.write(scaffold + "\n")

    print("Added scaffold to data/macrocycles/macrocycles.tsv")
    return


def main():
    while True:
        macrocycle = None

        while not macrocycle:
            smiles = input("Peptic Macrocycle SMILES: ")
            macrocycle = Chem.MolFromSmiles(smiles)

        scaffold_from_macrocycle(macrocycle)

        response = None
        while not response in ("y", "n"):
            response = input("Continue?(y/n):").lower()
            if response == "y":
                continue
            elif response == "n":
                break
            else:
                print("Invalid Response")

        if response == "n":
            break


if __name__ == "__main__":
    main()
    exit()

"O=C1CCCNC([C@@H](NC(CN2CCN(C([C@@H](CNC(CC3=CC=C(S(=O)(F)=O)C=C3)=O)NC([C@@H](CC4=CN=CC=C4)N1)=O)=O)CC2)=O)C(C)(C)C)=O"
