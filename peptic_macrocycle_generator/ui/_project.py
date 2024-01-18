from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
import py3Dmol
from datetime import datetime

from IPython.display import display


def save_3d(molecule, width=500, height=500):
    """
    Generate an HTML for 3D molecular visualization.

    Parameters:
    - molecule: RDKit molecule object
    - width (int): Width of the viewer (default: 800)
    - height (int): Height of the viewer (default: 400)

    Returns:
    - None
    """
    molecule = Chem.Mol(molecule)

    molecule = Chem.AddHs(molecule)

    # Generate 3D coordinates
    # Use a specific seed for reproducibility
    AllChem.EmbedMolecule(molecule, randomSeed=42)
    try:
        MMFFOptimizeMolecule(molecule)
    except:
        pass

    # Convert RDKit molecule to pdb format
    mol_block = Chem.MolToPDBBlock(molecule)

    # Create Py3Dmol view
    viewer = py3Dmol.view(width=width, height=height)

    # Add the molecule to the viewer
    viewer.addModel(mol_block, format='pdb')

    # Style and zoom
    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor('black')
    viewer.zoomTo()

    filename = datetime.now().strftime("%Y-%m-%d_%H-%M")

    viewer.write_html(f"./outs/3D/{filename}.html", fullpage=True)

    print(f"HTML file generated: outs/3D/{filename}.html")


def save_2d(mol):
    """Save a PNG image of the structure of given molecule(s)

    Args:
        mol (rdkit.Mol): Molecule to visualize
    """
    if isinstance(mol, (list, tuple)):
        img = Draw.MolsToGridImage(mol)
    else:
        img = Draw.MolToImage(mol)

    filename = datetime.now().strftime("%Y-%m-%d_%H-%M")

    img.save(f"./outs/2D/{filename}.png")

    print(f"PNG file generated: outs/2D/{filename}.png")


def project_2d(mol):
    """Show simple 2D visual of molecule structure

    Args:
        mol (rdkit.Mol): Molecule to be visualized
    """
    if isinstance(mol, (list, tuple)):
        img = Draw.MolsToGridImage(mol)
    else:
        img = Draw.MolToImage(mol)
    img.show()


def project_3d(molecule, width=500, height=500):
    """
    Generate a Py3Dmol viewer for 3D molecular visualization.

    Parameters:
    - molecule: RDKit molecule object
    - width (int): Width of the viewer (default: 800)
    - height (int): Height of the viewer (default: 400)

    Returns:
    - py3dmol.view: Py3Dmol viewer object
    """
    molecule = Chem.Mol(molecule)

    molecule = Chem.AddHs(molecule)

    # Generate 3D coordinates
    # Use a specific seed for reproducibility
    AllChem.EmbedMolecule(molecule, randomSeed=42)

    MMFFOptimizeMolecule(molecule)

    # Convert RDKit molecule to pdb format
    mol_block = Chem.MolToPDBBlock(molecule)

    # Create Py3Dmol view
    viewer = py3Dmol.view(width=width, height=height)

    # Add the molecule to the viewer
    viewer.addModel(mol_block, format='pdb')

    # Style and zoom
    viewer.setStyle({'sphere': {}})
    viewer.setBackgroundColor('black')
    viewer.zoomTo()

    viewer.show()
