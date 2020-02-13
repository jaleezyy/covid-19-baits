# Software Requirements

* [BaitsTools](https://github.com/campanam/BaitsTools) (requires Ruby 2.4.1+)
* [OligoArrayAux](http://unafold.rna.albany.edu/?q=DINAMelt/OligoArrayAux) (requires Perl 5+ )
* BLAST+ (specifically BLASTN)
* Python 3.6.8+ (with Biopython 1.7.5)

# header_trim.py
Script trims header of a FASTA file after running BaitsTools such that each header consists of reference sequence / start co-ordinate / stop co-ordinate.

# melting_probe_filter_WIP.py
To be used with melt_parse.sh, execute melt.pl (part of OligoArrayAux) which calculates melting temperature and parses for only melting temperature data. Once melting temperature information generated, filters to keep probes within specified minimum and maximum temperature ranges.

# melt_parse.sh
Execute melt.pl (requires OligoArrayAux) and parse resulting output such that only melting temperatures remain. Order of melting temperatures matches order of FASTA file input. Will be executed as part of melting_probe_filter_WIP.py.

# off_target_filter_WIP.py
Uses BLASTN output (-outfmt '10 std staxids sscinames sblastnames sskingdoms') to removes probes that do not align to specified keep term (built corresponding to sskingdoms). 

# binb4greedy_WIP.py
Self BLASTN of candidate probes, group all hits corresponding to a query probe and remove probes covering redundant sequence space. 

