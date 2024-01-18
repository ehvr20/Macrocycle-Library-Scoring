import pandas as pd 
from concurrent.futures import ProcessPoolExecutor

from ..chemistry import properties_summary

def generate_properties_summaries(molecules, max_workers=None):
    """A generator function to calculate the physiochemical properties of a list of molecules

    Args:
        molecules (list(rdkit.Mol)): List of molecules
        max_workers (int, optional): Max #Cores. Defaults to None.

    Returns:
        generator(-->dict): A generator of dictionary values
    """    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Use executor.map to parallelize the generation
        summaries = executor.map(properties_summary, molecules)

    return summaries