import re

from collections import defaultdict
from itertools import combinations, product
from loguru import logger
from rdkit.Chem import AllChem, ReplaceCore, GetMolFrags, ReplaceSubstructs, CombineMols


from substituent_replacement.models.library import ReplacementLibrary, ReplacementMap
from substituent_replacement.utils.utils import canonicalise_smiles, prep_mol


class QueryDerivative:
    """Class to store info regarding how a derivative is formed from QuerySubstituents"""
    def __init__(self, sub_map):
        self.id = None
        self.smiles = None
        self.mol = None
        self.highlight = defaultdict(list)

        # Substituents
        self.substituents = sub_map
        self.num_replacements = len([sub for sub in self.substituents.values() if isinstance(sub, ReplacementMap)])

        # Compute score
        self.score = self._compute_score()

    def __repr__(self):
        return f"<Derivative: {self.id}, {self.smiles}>"

    def _compute_score(self) -> float:
        score = 0
        for replacement in self.substituents.values():
            if isinstance(replacement, ReplacementMap):
                score += replacement.edge_weight
        score /= self.num_replacements**2
        return round(score, 1)


class QuerySubstituent:
    """Class to store substituent info of the Query Molecule"""
    def __init__(self, mol, sub_id):
        self.id = sub_id
        self.mol = mol
        self.smiles = canonicalise_smiles(AllChem.MolToSmiles(mol))
        self.atom_indices = [x.GetIdx() for x in self.mol.GetAtoms()]
        self.attachment_index = max([x.GetIsotope() for x in self.mol.GetAtoms()])

    def __repr__(self):
        return f"<Substituent: {self.id}, {self.smiles}>"


class Molecule:
    """Class to store the Query Molecule"""
    def __init__(self, smiles: str, core: str):
        logger.info(f'Initializing Molecule...')
        self.smiles = canonicalise_smiles(smiles)
        self.core_smiles = canonicalise_smiles(core, use_smarts=True)
        self.mol = AllChem.MolFromSmiles(self.smiles)
        self.core_mol = AllChem.MolFromSmarts(self.core_smiles)
        logger.info(f'    SMILES: {self.smiles}')
        logger.info(f'    Core:   {self.core_smiles}')

        # Generate substituents
        self.substituents = self._generate_substituents(mol=self.mol, core_mol=self.core_mol)

        # Derivatives
        self.derivatives = []

        logger.info(f'    Initialization Complete.')

    def _construct_derivatives(self):
        """Connects a replacement smiles to a core scaffold to create a RDKit Mol object"""

        unique_smiles = set()
        duplicates = set()
        for deriv in self.derivatives:
            mol = AllChem.RWMol(self.core_mol)
            for query_substituent, replacement in deriv.substituents.items():

                # Get the replacement or substituent that will be attached to the core
                if isinstance(replacement, ReplacementMap):
                    rep_mol = replacement.replacement_mol
                elif isinstance(replacement, QuerySubstituent):
                    rep_mol = prep_mol(replacement.smiles)
                else:
                    raise TypeError(f'I have no idea what this is')

                # Attach replacement substituent to core
                mol = AllChem.RWMol(CombineMols(mol, rep_mol))
                end_atm = -1
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum() == 1:
                        end_atm = atom.GetIdx()
                mol.AddBond(query_substituent.attachment_index, end_atm, order=AllChem.rdchem.BondType.SINGLE)
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)

            # Create Smiles and Mol objects
            smiles = canonicalise_smiles(AllChem.MolToSmiles(mol))
            if smiles not in unique_smiles:  # Ensure only unique smiles are generated
                unique_smiles.add(smiles)
                deriv.smiles = canonicalise_smiles(AllChem.MolToSmiles(mol))
                deriv.mol = AllChem.MolFromSmiles(deriv.smiles)

                # Determine atom numbers of replaced substituents
                # 1. Get sidechain smiles
                sidechain_mol = ReplaceCore(deriv.mol, self.core_mol, labelByIndex=True)
                sidechain_smis = [canonicalise_smiles(AllChem.MolToSmiles(frag))
                                  for frag in GetMolFrags(sidechain_mol, asMols=True)]

                # 2. Create temporary mol object for derivative and label smiles with their atom numbering
                temp_mol = AllChem.MolFromSmiles(deriv.smiles)
                sidechain_atoms = []
                for atom in temp_mol.GetAtoms():
                    atom.SetAtomMapNum(atom.GetIdx() + 1)  # Needs to +1 because 0 is not labeled in the SMILES

                # 3. Find the substituents of the derivative using the core smarts. The substituent smiles returned
                #    contain the numbering of the substituent relative to the parent derivative structure. Use string
                #    matching to find the atom numbers that will be highlighted
                temp_mol = ReplaceCore(temp_mol, self.core_mol, labelByIndex=True)
                for frag in GetMolFrags(temp_mol, asMols=True):
                    smi = canonicalise_smiles(AllChem.MolToSmiles(frag))
                    atoms = re.findall(r':(\d+)]', smi)  # Find atom numbering in the parent derivative
                    atoms = [int(atom) - 1 for atom in atoms]  # Needs to -1 because above indexing is +1
                    sidechain_atoms.append(atoms)
                sidechain_map = dict(zip(sidechain_smis, sidechain_atoms))

                # 4. Using the atoms identified above, find bonds between atoms in the parent derivative and save
                #    these as dictionary entries that can be use in highlighting in SVG images
                for _, repl in deriv.substituents.items():
                    if isinstance(repl, ReplacementMap) and sidechain_map.get(repl.replacement_smiles, None):
                        atoms = sidechain_map[repl.replacement_smiles]
                        deriv.highlight['highlightAtoms'].extend(atoms)
                        for bond in deriv.mol.GetBonds():
                            if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms:
                                deriv.highlight['highlightBonds'].append(bond.GetIdx())
            else:
                duplicates.add(deriv)

        # Ensure only unique derivatives are kept
        all_derivs = set(self.derivatives.copy())
        self.derivatives = list(all_derivs - duplicates)
        logger.info(f'    {len(self.derivatives)} unique derivatives generated')

    @staticmethod
    def _generate_substituents(mol, core_mol):
        logger.info(f'    Generating Substituents...')
        substituents = []
        sidechain_mol = ReplaceCore(mol, core_mol, labelByIndex=True)
        sub_index = 0
        for sub in GetMolFrags(sidechain_mol, asMols=True):
            sub_index += 1
            substituents.append(QuerySubstituent(mol=sub, sub_id=sub_index))
        logger.info(f'        {len(substituents)} substituents found')
        return substituents

    def generate_derivatives(self, replacement_map: ReplacementLibrary, level: int = 1):
        """Generate derivatives based on substituents in the Query Molecule and the replacement options in the
        substitution_tree dictionary. Derivatives are stored in self.derivatives"""

        # Get replacements
        logger.info(f'Generating Derivatives using up to level {level} as replacements...')
        replacements = {substituent: [(substituent, rep)
                                      for rep in replacement_map.get_replacements(query_smiles=substituent.smiles,
                                                                                  level=level)]
                        for substituent in self.substituents}
        logger.info(f'    {len(replacements)}/{len(self.substituents)} substituents are replaceable')

        # Generate possible combinations of query substituents
        # E.g. replacements = {sub_a: [(sub_a, 1), (sub_a, 2), (sub_a, 3)],
        #                      sub_b: [(sub_b, 4)],
        #                      sub_c: [(sub_c, 5), (sub_c, 6)]}
        #      possibilities = [(sub_a,), (sub_b,), (sub_c,),
        #                       (sub_a, sub_b), (sub_a, sub_c), (sub_b, sub_c),
        #                       (sub_a, sub_b, sub_c)]
        possibilities = []
        for i in range(1, len(replacements.keys()) + 1, 1):
            possibilities.extend(list(combinations((list(replacements.keys())), i)))

        # Enumerate all combos for combinations of query substituents
        # all_combos = [((sub_a, 1),), ((sub_a, 2),), ((sub_a, 3),), ((sub_b, 4),), ((sub_c, 5),), ((sub_c, 6),),
        #               ((sub_a, 1), (sub_b, 4)), ((sub_a, 2), (sub_b, 4)), ((sub_a, 3), (sub_b, 4)),...]
        all_combos = []
        for p in possibilities:
            combine = [replacements[sub] for sub in p]
            all_combos.extend(list(product(*combine)))
        all_combos = set(all_combos)

        logger.info(f'    {len(all_combos)} derivatives are possible')

        # Generate derivatives
        for combo in all_combos:
            # Create substituent-derivative map
            deriv_map = {query_sub: query_sub for query_sub in self.substituents}
            for substituent, replacement in combo:
                deriv_map[substituent] = replacement

            # ToDo: Determine if combination of smiles is unique - do this better
            deriv = QueryDerivative(sub_map=deriv_map)
            self.derivatives.append(deriv)

        self.derivatives.sort(key=lambda d: d.score, reverse=True)  # Sort derivatives by score, highest returned first
        deriv_index = 1
        for deriv in self.derivatives:
            deriv.id = deriv_index
            deriv_index += 1

        self._construct_derivatives()
