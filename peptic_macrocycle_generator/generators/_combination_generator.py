from itertools import product, chain

def generate_residue_combinations(scaffolds, building_blocks):
    """Generator function for all possible residue combinations from a list of macrocycle scaffolds and building_blocks.

    Args:
        scaffolds (list( [scaffold] | * )): List of scaffolds
        building_blocks (list(rdkit.Mol)): List of building blocks 
    
    Returns:
        generator(--> list): A generator of all possible residues combinations
    """    
    generators = []

    for scaffold in scaffolds:

        combinations = [building_blocks if residue == "*" else residue for residue in scaffold]
        combinations = product(*combinations)
        generators.append(combinations)

    return chain(*generators)
        