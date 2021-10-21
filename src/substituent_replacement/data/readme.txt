R-group replacement database
Frequently used R-groups and replacement hierarchies are provided in an XML file. Root has 500 “center” child elements with two attributes, “degree” for node degree in the final R-group network and “SMILES” (R-group structure representation).  
In addition, “center” elements have sub-elements as s “first layer” followed by a “second layer” (if available).  Each layer has “SMILES” and “edge_weight”.
The collection of frequently used R-groups can be searched for R-groups of interest and their preferred replacements using SMILES representations. For example, to search for the hydroxyl group, search the corresponding SMILES representation in the SMILES attribute of the “center” element.  The resulting output will be five first layers and max. two associated second layers per first layer, for example:
<center degree="13627" SMILES="*O">
	<first_layer SMILES="*OC" edge_weight="2273">
		<second_layer SMILES="*Cl" edge_weight="3455"/>
		<second_layer SMILES="*OCC" edge_weight="788"/>

Final R-group network
Frequently used R-groups and preferred substitution site-specific replacements in analogue series were identified on the basis of bioactive compounds from the ChEMBL (release 26) database.  R-groups and replacements are provided in a csv file that can be used to build a network data structure.
In the csv file, “source" and “target“ represent two nodes forming edge (pairwise R-group replacement). R-group structures are represented as SMILES strings. In addition, “edge_weight" gives the frequency of occurrence of an R-groups replacement across all substitution sites.
