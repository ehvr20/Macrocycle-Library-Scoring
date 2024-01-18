from rdkit import Chem
from rdkit.Chem import AllChem

def remove_Fmoc(mol):
    """Remove Fmoc group from amine

    Args:
        mol (rdkit.Mol): Molecule

    Returns:
        rdkit.Mol: Molecule with Fmoc group removed
    """    
    Fmoc = Chem.MolFromSmiles('O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2')

    clean_molecule = AllChem.DeleteSubstructs(mol, Fmoc)

    return clean_molecule