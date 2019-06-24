# TherapeuticABSearch


This script allows to perform sequence searches on the therapeutic antibodies curated in <a href="" target="_blank">BISDAB</a>.


## Prerequisites.

You need to install the free ANARCI software available from here.


## Usage.

Provide the amino acid sequence either the heavy chain, light chain or both and a directory where the results should be output.

A search will be performed, aligning the therapeutic antibodies to your query sequences using IMGT numbering.

The results are stored in the output directory that you provided. If you did not provide an output file, these will be stored in `search_out` in the current directory.


## Output.

The results for the heavy and light chains are given in the following files:

* `heavy/light_sequence_identities.tsv` : these are the tab-separated files giving the % sequence identities to the entire variable region CDRs and framework regions as defined by IMGT.
* `heavy/light_sequence_alignments.txt` : these are the alignments of your query sequences against the therapeutic heavy/light chains. Aligned positions are connected via `|`, misaligned via `.` and the CDRs are indicated by `^`. 
