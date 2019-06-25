# TherapeuticABSearch

Author: Konrad Krawczyk, <a href="http://naturalantibody.com" target="_blank">NaturalAntibody</a>.

This script allows to perform sequence searches on the therapeutic antibodies curated in <a href="http://opig.stats.ox.ac.uk/webapps/newsabdab/therasabdab/search/" target="_blank">BISDAB</a>.


## Prerequisites.

You need to install the free antibody numbering software <a href="http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php" target="_blank">ANARCI</a> software available from here.


## Usage.

Provide the amino acid sequence either the heavy chain, light chain or both and a directory where the results should be output.

`python Search.py  --heavy-sequence QVQLQQSGSELKKPGASVKVSCKASGYTFTNYGMNWVKQAPGQGLKWMGWINTYTGEPTYTDDFKGRFAFSLDTSVSTAYLQISSLKADDTAVYFCARGGFGSSYWYFDVWGQGSLVTVSS --light-sequence DIQMTQSPSSLSASVGDRVTITCSASQDISNYLNWYQQKPGKAPKVLIYFTSSLHSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYSTVPWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCL`

A search will be performed, aligning the therapeutic antibodies to your query sequences using IMGT numbering.

The results are stored in the output directory that you provided. If you did not provide an output file, these will be stored in `search_out` in the current directory.

The full list of options that you can specify is listed below:

* `--heavy-sequence` : sequence of the heavy chain.
* `--light-sequence` : sequence of the light chain.
* `--output` : where the output should be written to, default `search_out`
* `--cutoff` : sequence identity cutoff, specified as an integer, default `80`
* `--region` : specific region you would like to constrain on available regions are `V [default], cdr1, cdr2, cdr3, fw1, fw2, fw3, fw4`


## Output.

The results for the heavy and light chains are given in the following files:

* `heavy/light_sequence_identities.tsv` : these are the tab-separated files giving the % sequence identities to the entire variable region CDRs and framework regions as defined by IMGT.
* `heavy/light_sequence_alignments.txt` : these are the alignments of your query sequences against the therapeutic heavy/light chains. Aligned positions are connected via `|`, misaligned via `.` and the CDRs are indicated by `^`. 
* `combined.tsv` : if you specified both heavy and light chains and are not constraining to a subregion (e.g. cdr1), the names of therapeutics where both the heavy and light chains match the search criteria are given with their respective sequence identities.
