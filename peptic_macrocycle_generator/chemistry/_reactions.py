from rdkit import Chem
from rdkit.Chem import AllChem

from ..ui import project_2d

def peptide_synthesis(residues):
    """Function to synthesize a peptide molecule

    Args:
        residues (list(rdkit.Mol)): List of residues for condensation

    Returns:
        rdkit.Mol: Peptide Molecule
    """    

    if not len(residues) >=2:
        raise(TypeError, "At least two amino acids required")

    nh2_reaction = AllChem.ReactionFromSmarts("[C:1](=[O:2])O.[NH2:3]>>[C:1](=[O:2])[N:3]")
    nh1_reaction = AllChem.ReactionFromSmarts("[C:1](=[O:2])O.[NH1:3]>>[C:1](=[O:2])[N:3]")

    chain = residues[0]
    residues = residues[1:][::-1]

    for residue in residues:
        reactants = (chain, residue)
        try:
            chain = (nh2_reaction.RunReactants(reactants)[0][0])
        except(IndexError):
            chain = (nh1_reaction.RunReactants(reactants)[0][0])

    return chain

def peptide_cyclization(molecule):
    """Cyclization of a peptide molecule

    Args:
        molecule (rdkit.Mol): Peptide Molecule to be cyclized

    Returns:
        rdkit.Mol: Cyclic Peptide (macrocycle)
    """    
    molecule = Chem.RemoveAllHs(molecule)

    # In the case of know NH2 group - find the index of NH
    try:
        amine_pattern = Chem.MolFromSmarts('[N;H2]')
        amine_index = molecule.GetSubstructMatches(amine_pattern)[0][0]

    except(IndexError):
        amine_pattern = Chem.MolFromSmarts('[N;H1]')
        amine_index = molecule.GetSubstructMatches(amine_pattern)[0][0]

    # Hydroxyl index must lie within a known carboxyl index
    hydroxyl_pattern = Chem.MolFromSmarts('[O;H1]')
    carboxyl_pattern = Chem.MolFromSmarts("[C,c;X3](=[O,S,P])-[O;H1]")
    carboxyl_indexes = []
    for match in molecule.GetSubstructMatches(carboxyl_pattern):
        for idx in match:
            carboxyl_indexes.append(idx)

    hydroxyl_index = molecule.GetSubstructMatches(hydroxyl_pattern)
    for match in hydroxyl_index:
        if match[0] in carboxyl_indexes:
            hydroxyl_index = match[0]
            break
        else:
            continue 

    # Get the atom of known index
    carboxyl_carbon_index = molecule.GetAtomWithIdx(hydroxyl_index).GetNeighbors()[0].GetIdx()

    mw = Chem.RWMol(molecule)

    mw.RemoveAtom(hydroxyl_index)

    mw.AddBond(amine_index,carboxyl_carbon_index,Chem.BondType.SINGLE)

    macrocycle = Chem.Mol(mw)

    try:
        macrocycle.UpdatePropertyCache()
    except:
        project_2d(macrocycle)
        exit()

    return macrocycle

def hydrolysis(macrocycle):
    """Hydrolyze a macrocycle into it's individual residues.

    Args:
        macrocycle (rdkit.Mol): Macrocycle Molecule

    Returns:
        list(rdkit.Mol): List of Macrocycle Residues
    """
    peptide_bond_pattern = Chem.MolFromSmarts("NC=O")
    peptide_bonds = macrocycle.GetSubstructMatches(peptide_bond_pattern)

    molecule = Chem.RWMol(macrocycle)
    molecule.BeginBatchEdit()

    bonds_to_hydrolyze = []

    for peptide_bond in peptide_bonds:
        bond_indices = [None, None]

        for index in peptide_bond:
            atom = molecule.GetAtomWithIdx(index)

            if atom.GetAtomicNum() == 6:
                bond_indices[0] = index
            elif atom.GetAtomicNum() == 7:
                bond_indices[1] = index
            elif atom.GetAtomicNum() == 8:
                continue

        carbon = molecule.GetAtomWithIdx(bond_indices[0])

        if carbon.IsInRing():
            bonds_to_hydrolyze.append(bond_indices)

    for carbon_idx, nitrogen_idx in bonds_to_hydrolyze:
        molecule.RemoveBond(carbon_idx, nitrogen_idx)
        oxygen_idx = molecule.AddAtom(Chem.Atom(8))
        molecule.AddBond(carbon_idx, oxygen_idx, Chem.BondType.SINGLE)


    molecule.CommitBatchEdit()

    molecule.UpdatePropertyCache()

    residues = Chem.GetMolFrags(molecule, asMols=True)

    return residues
