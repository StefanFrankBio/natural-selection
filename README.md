# Natural Selection
## Issues
- test_variant_records() can produce substitutions from one base to itself
- test.py main() looks messy
- Add type hinting to functions

## Next Steps
- In run.sh add ability to handle individual fastas instead of multifasta for variant sequences
    - write test variants to individual files instead of multifasta
- In run.sh move removal of previous analysis data to front of script
- Add function that writes (multi-)fasta files using Bio.SeqIO
- Reconstruct variants from variant records and compare with test variant records
    - Automate testing on github
- Build SQL database and add variant records, two approaches:
    - Either one table each per variant
        - could be integrated into variant_record.py
    - Or one table containing all positions and variations
    - Test parallelization with both approaches
- Add metadata as SQL table
- Add functionality to calculate dN/dS for whole genome in all frames
- Add functionality to calculate dN/dS for a specified sampling date or date range
- Add functionality to calculate dN/dS for stop-stop-ORFs
- Add functionality to calculate dN/dS for sliding window over genome
- Add functionality to read gene locations from annotation and calculate dN/dS for them
    - Add functionality to annotate genomes

## Research Questions
How high is the likelihood for a nonsynonymous substitution to occure at a nonsynonymous site?

How many nonsynonymous substitutions would be expected for a given sequence?

How are nonsynonymous sites spread to through a given sequence?

Do nonsynonymous substitutions spread in a similar manner to nonsynonymous sites?

Can records for multiple variants be combined to calculate dN/dS for a metadata range?

Can dN/dS of noncoding regions be used to establish baseline values for dN and dS and is this useful for error correction?