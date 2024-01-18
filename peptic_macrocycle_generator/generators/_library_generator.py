from ..chemistry import peptide_cyclization, peptide_synthesis
from concurrent.futures import ProcessPoolExecutor

# def generate_library(residue_combinations):

#     for combination in residue_combinations:

#         peptide = peptide_synthesis(combination)
        
#         macrocycle =  peptide_cyclization(peptide)
        
#         yield macrocycle

def process_combination(combination):
    """A simple function to turn a residue combination into a macrocycle molecule

    Args:
        combination (list(rdkit.Mol)): A list of residues

    Returns:
        rdkit.Mol: A macrocycle molecule
    """    

    peptide = peptide_synthesis(combination)
    macrocycle = peptide_cyclization(peptide)

    return macrocycle

def generate_library(residue_combinations):
    """A Generator Function that return all possible macrocycle molecules

    Args:
        residue_combinations (list(list(rdkit.Mol))): A list of all possible residue combinations

    Returns:
        generator(--> rdkit.Mol): macrocycle molecule
    """    
    with ProcessPoolExecutor(max_workers=None) as executor:
        # Use executor.map to parallelize the generation
        library = executor.map(process_combination, residue_combinations)

        return library
