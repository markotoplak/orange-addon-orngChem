fragmenter.py
	Takes an input file containing SMILES codes and finds frequent fragments in them
	Arguments:

	-s file		

	specifies an input file(default "smiles.tab") containing
	SMILES codes. Can be a tab delim. file with SMILES as a second
	column or a simple txt file containg SMILES strings seperated
	with newline

	-a attrName

	the name of the attribute that contains the smiles codes in
	the input file

	-o file

	specifies an output file (default "fragments.txt") where it
	will print one fragment code per line

	-f freq

	specifies the required frequency for the fragments (default
	0.4)

	Example:
	python fragmenter.py -s smiles.tab -a SMILES -o fragmentcodes.txt -f 0.4

fragment_mapper.py
	Takes an input file containg SMILES codes of chemicals and an file contating fragment codes. Outputs an 
	tab delim. file specifying the 	number of  times a fragment is matched in each chemical.
	Arguments:

	-s file

	specifies an input file (default "smiles.tab") containing
	SMILES codes. Must be a tab delim. file with SMILES as second
	column or a column specified by -a attrName

	-a attrName

	the name of the attribute that contains the smiles codes in
	the input file

	-f file

	specifies an input file (default "fragments.txt") containing
         fragment codes seperated by newlines

	-o file

	specifies an output file (default "fragmentmap.tab")

	-b

	if present the output file will only contain values 0 or 1
	indicating the presence of a fragment in a chemical

profile_analisys.py
	Takes the sensitivity data and maps it to go terms or fragments and go terms
	Arguments:

	-s file

	specifies an input file(default "smiles.tab") containing
	SMILES codes. Must be a tab delim. file with SMILES as second
	column

	-S file		(specifies an file containing sensitivity data(default "sens.tab"))
	-b file

	specifies an file containing subset of go terms, that we are
	interested in (otherwise all terms are considered)

	-o file		(output file (default "sens_profile.tab")
	-f file		(file containing fragment codes)
	-a P|C|F	(go aspect (process|component|function) we are interested in (default P))
	-g path		(path to the local go database)
	-l 		(use slims subset only)
	-m		(output molecule based profile, otherwise fragment based profile will be generated)
	Example:
	python profile_analisys.py -s smiles.tab -S sens.tab -f fragments.txt -o ms.tab -b selected_terms.txt -m
	
profile_transpose.py
	Flips the output from profile_analisys.py
	Arguments:
	-i file		(input filename(default "sens_profile.tab")
	-o file		(output filename(default "sens_profile_t.tab")
	-n name		(meta variable name (default "name"))
	Example:
	python profile_transpose.py -i sens_profile.tab -n smiles

map_to_closest_terms.py
	Takes an input file ("Molecule-based slim profiles" output from 
	OWSensProfile) and maps terms to closest chemicals and vice versa.
	Arguments:
	-i file 	(specifies an input file)
	-k n 		(map n closest terms to each chemical, default 5)
	-j n 		(map n closest chemicals to each term, default 3)
	-o1 file	(specifies an output file for term to chemical mapping, default "term-chem.tab")
	-o2 file	(specifies an output file for chemical to term mapping, default "chem-term.tab")
	
	Example:
	python map_to_closest_terms.py -i ms.tab -k 5 -j 5 

map_genes_to_terms.py
	For each specified term get all genes from the whole
	annotation and only those from the Boone paper (reads those
	genes from genes.txt). Outputs a file named
	"terms_with_genes.txt".
	
	Example:
	python map_genes_to_terms.py

plot.py
	Draws the CA plot of any OWSensProfile or profile_analisys.py output (Molecule-based slim terms ...)
	Example:
	python plot.py ms.tab

OWSensProfile widget
	input: sensitivity data from Boone ("sens.tab")
