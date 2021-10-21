import os

from rdkit.Chem import AllChem
from typing import Union

import substituent_replacement.data as data_path

from substituent_replacement.models.molecule import Molecule
from substituent_replacement.models.library import ReplacementLibrary
from substituent_replacement.utils.utils import draw_molecule


IN_FILE = './source_data/top500_R_replacements.xml'
SMILES_TARGET = "c1cc(ccc1c2c(n(cn2)CC3CC3)c4ccnc(n4)N)F"  # Pat Walters example
CORE_SMARTS = "n1cncc1"  # Pat Walters example
OUT_DIR = 'results'
LEVEL = 1  # Options are 1 and 2


def substituent_replacement(target_smiles: str, core_smarts: str, out_folder: str, level: int = 1,
                            source_path: Union[str, None] = None):
    """API to perform substituent replacement on a target molecule based on a core scaffold and a library of
    transformations"""
    # Setup results folder
    os.makedirs(out_folder, exist_ok=True)

    # Load Replacement Database
    if source_path is None:
        source_path = f'{data_path.__file__}/top500_R_replacements.xml'

    replacement_library = ReplacementLibrary(xml_path=source_path)

    # Prepare Query Molecule
    molecule = Molecule(smiles=target_smiles, core=core_smarts)
    molecule.generate_derivatives(replacement_map=replacement_library, level=level)

    # Visualise mols
    svg_out = os.path.join(out_folder, 'svgs')
    os.makedirs(svg_out, exist_ok=True)

    # Draw parent molecule
    legend = {'legend': 'Parent'}
    draw_molecule(mol=molecule.mol, template=molecule.core_mol, out_path=os.path.join(svg_out, 'mol_0.svg'), **legend)

    # Output
    tsv_out = os.path.join(out_folder, 'derivatives.tsv')
    with open(tsv_out, 'w') as w:
        w.write('ID\tSMILES\tNUM_REPLACEMENTS\tSCORE\n')
        for deriv in molecule.derivatives:
            # Save data to TSV
            w.write(f'{deriv.id}\t{deriv.smiles}\t{deriv.num_replacements}\t{deriv.score}\n')

            # Draw derivatives to SVG files
            kwargs = {'legend': str(deriv.score)}
            kwargs.update(deriv.highlight)
            draw_molecule(mol=AllChem.MolFromSmiles(deriv.smiles),
                          template=AllChem.MolFromSmarts(core_smarts),
                          out_path=os.path.join(svg_out, f'mol_{deriv.id}.svg'),
                          **kwargs)


if __name__ == '__main__':
    substituent_replacement(source_path=IN_FILE, target_smiles=SMILES_TARGET,
                            core_smarts=CORE_SMARTS, out_folder=OUT_DIR, level=LEVEL)
