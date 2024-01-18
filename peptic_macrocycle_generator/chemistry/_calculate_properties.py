from rdkit import Chem
from rdkit.Chem import QED


def properties_summary(molecule):
    """Calculate physiochemical properties of an rdkit.Mol

    Args:
        molecule (rdkit.Mol): Molecule

    Returns:
        dict: dictionary of physiochemical properties
    """    
    score = QED.default(molecule)
    attributes = QED.properties(molecule)
    properties = {
        'Score': score,
        'MW': attributes[0],
        'ALOGP': attributes[1],
        'HBA': attributes[2],
        'HBD': attributes[3],
        'PSA': attributes[4],
        'ROTB': attributes[5],
        'AROM': attributes[6],
        'ALERTS': attributes[7],
        'SMILES': Chem.MolToSmiles(molecule)
    }
    return properties
