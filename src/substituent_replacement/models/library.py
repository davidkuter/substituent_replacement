import xmltodict

from loguru import logger
from rdkit.Chem.rdchem import KekulizeException
from typing import List

from substituent_replacement.utils.utils import canonicalise_smiles, prep_mol


class ReplacementMap:
    """Class that maps substituents to be replaced by other substituent"""

    def __init__(self, map_id, target_smiles: str, replacement_smiles: str, edge_weight: int, level: int):
        self.id = map_id
        self.edge_weight = int(edge_weight)
        self.level = level
        self.target_smiles = canonicalise_smiles(smiles=target_smiles)
        self.replacement_smiles = canonicalise_smiles(smiles=replacement_smiles)
        self.replacement_mol = prep_mol(smiles=self.replacement_smiles)

    def __repr__(self):
        return f"<Replacement Map: ID = {self.id}, level = {self.level}, edge_weight = {self.edge_weight}, " \
               f"target_smiles = {self.target_smiles}, replacement_smiles = {self.replacement_smiles}>"


class ReplacementLibrary:
    """Class to store substituent info for moieties that will be replaced or that will replace other substituents"""

    def __init__(self, xml_path: str):
        logger.info(f'Loading Replacment Library from: {xml_path}')
        with open(xml_path) as f:
            xml = f.read()
            contents = xmltodict.parse(xml)

        self.replacement_map = []
        map_index = 0
        for center in contents['R_replacements']['center']:
            center_smiles = canonicalise_smiles(smiles=center['@SMILES'])
            # Get Substituent replacements in the first layer
            first_layer = center['first_layer']
            if isinstance(first_layer, list) is False:
                first_layer = [first_layer]

            for sub1 in first_layer:
                sub1_smiles = canonicalise_smiles(sub1['@SMILES'])
                # Check if replacement already exists
                try:
                    replacement = self._avoid_duplicates(target=center_smiles, replacement=sub1_smiles,
                                                         ew=sub1['@edge_weight'], level=1, map_index=map_index)
                except KekulizeException:
                    logger.warning(f'    Kekulization failed for {sub1_smiles}. Omitting.')
                    replacement = None

                if replacement:
                    map_index += 1
                    self.replacement_map.append(replacement)

                # Get Substituent replacements in the second layer (if any)
                if sub1.get('second_layer', None):
                    second_layer = sub1['second_layer']
                    if isinstance(second_layer, list) is False:
                        second_layer = [second_layer]

                    for sub2 in second_layer:
                        sub2_smiles = canonicalise_smiles(sub2['@SMILES'])
                        # Check if replacement already exists
                        try:
                            replacement = self._avoid_duplicates(target=sub1_smiles, replacement=sub2_smiles,
                                                                 ew=sub1['@edge_weight'], level=2, map_index=map_index)
                        except KekulizeException:
                            logger.warning(f'    Kekulization failed for {sub2_smiles}. Omitting.')
                            replacement = None

                        if replacement:  # If the is no duplicate
                            map_index += 1
                            self.replacement_map.append(replacement)

        logger.info(f'    {len(self.replacement_map)} substituents loaded comprising:')
        logger.info(f'        {len([rep for rep in self.replacement_map if rep.level == 1])} level one replacements')
        logger.info(f'        {len([rep for rep in self.replacement_map if rep.level == 2])} level two replacements')

    def _avoid_duplicates(self, target: str, replacement: str, ew: int, level: int, map_index):
        for rep in self.replacement_map:
            if rep.target_smiles == target and rep.replacement_smiles == replacement and \
                    rep.edge_weight == ew and rep.level == level:
                return None

        return ReplacementMap(target_smiles=target, replacement_smiles=replacement, edge_weight=ew,
                              level=level, map_id=map_index)

    def get_replacements(self, query_smiles: str, level: int = 1) -> List[ReplacementMap]:
        """Finds all substituents that should replace the query smiles. Level indicates how many levels you want to
        consider substituent replacements for"""

        if level < 0 or level > 2:
            raise TypeError(f'Level provided ({level}) can only be 1 or 2')

        replacements = [rep for rep in self.replacement_map if rep.target_smiles == query_smiles and rep.level == 1]
        # If the second level is required
        if level == 2:
            temp = [rep.replacement_smiles for rep in replacements]
            level_2 = [rep for rep in self.replacement_map if rep.target_smiles in temp and rep.level == 2]
            replacements.extend(level_2)

        return replacements
