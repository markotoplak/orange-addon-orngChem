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
	For each specified term get all genes from the whole annotation and only those from the 
	Boone paper (reads those genes from genes.txt). Outputs a file named "terms_with_genes.txt".
	
	Example:
	python map_genes_to_terms.py

plot.py
	Draws the CA plot of any OWSensProfile output (Molecule-based slim terms ...)
	Example:
	python plot.py ms.tab

OWSensProfile widget
	input: sensitivity data from Boone ("sens.tab")
	
	

	