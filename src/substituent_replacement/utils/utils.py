from rdkit.Chem import AllChem, Draw, rdFMCS
from rdkit.Chem.rdchem import Mol
from typing import Union


def canonicalise_smiles(smiles: str, use_smarts: bool = False) -> Union[str, None]:
    """Canonicalise smiles, setting all isotopes to natural values"""

    if use_smarts:
        mol = AllChem.MolFromSmarts(smiles)  # Required for core
    else:
        mol = AllChem.MolFromSmiles(smiles)

    if mol:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
        if use_smarts:
            return AllChem.MolToSmarts(mol)
        else:
            return AllChem.MolToSmiles(mol)


def draw_molecule(mol: Mol, out_path: str, template: Mol = None, **kwargs):
    """Creates a image of molecule. Defaults to svg"""

    # Validate input
    if not isinstance(mol, Mol):
        raise TypeError(f'mol argument provided must be a RDKit Mol object. Type provided: {type(mol)}')

    # Prepare molecule 2D structure
    AllChem.Compute2DCoords(mol)

    if template:    # Align with a template if present
        if not isinstance(template, Mol):  # Validate input
            raise TypeError(f'template argument provided must be a RDKit Mol object. Type provided: {type(template)}')

        AllChem.Compute2DCoords(template)
        AllChem.GenerateDepictionMatching2DStructure(mol, template)

    # Draw
    Draw.MolToFile(mol, out_path, **kwargs)


def prep_mol(smiles: str) -> Mol:
    """Prepares and converts replacement substituent smiles into a mol object"""
    # This is pretty much verbatim from Pat Walter's example. Not sure why these exact steps are needed
    mol = AllChem.MolFromSmiles(smiles)
    rw_mol = AllChem.RWMol(mol)
    remove_idx = -1
    for atm in rw_mol.GetAtoms():
        if atm.GetAtomicNum() == 0:
            remove_idx = atm.GetIdx()
            for nbr in atm.GetNeighbors():
                nbr.SetAtomMapNum(1)
    rw_mol.RemoveAtom(remove_idx)
    AllChem.SanitizeMol(rw_mol)
    return rw_mol
