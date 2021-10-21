import os

from loguru import logger
from rdkit.Chem import AllChem
from typing import Union

import substituent_replacement.data as data_path

from substituent_replacement.models.molecule import Molecule
from substituent_replacement.models.library import ReplacementLibrary
from substituent_replacement.utils.utils import draw_molecule


def substituent_replacement(target_smiles: str, core_smarts: str, out_folder: str, level: int = 1,
                            source_path: Union[str, None] = None):
    """API to perform substituent replacement on a target molecule based on a core scaffold and a library of
    transformations"""
    # Setup results folder
    os.makedirs(out_folder, exist_ok=True)

    # Load Replacement Database
    if source_path is None:
        source_path = data_path.__file__.replace('__init__.py', 'top500_R_replacements.xml')

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
    logger.info(f'Saving data to: {out_folder}')
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
