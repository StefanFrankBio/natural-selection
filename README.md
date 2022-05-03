# Natural Selection
## Issues
- Add type hinting to functions
- Testing: Instead of Gap-Match-Gap Stretcher Scoring prefers Gap-GapExtension-Mismatch. Leads to failed tests. Is that a Problem?
- In random variable records generated turing testing, a mismatch and an insert at consecutive positions might be switched around by stretcher
- If a variant record contains no synonymous substitutions for a frame, it's dN/dS will be div/0

## Next Steps
- Add BLAST as alternative for stretcher
- Add mean position of nonsynonymous sites vs mean position of nonsynonymous substitutions as additional metric
- <del>In run.sh add ability to handle individual fastas instead of multifasta for variant sequences</del>
    - <del>write test variants to individual files instead of multifasta</del>
- <del>In run.sh move removal of previous analysis data to front of script</del>
- Add function that writes (multi-)fasta files using Bio.SeqIO
- <del>Reconstruct variant sequences from variant records and compare with test variant sequence</del>
    <del>- Automate testing on github</del>
- <del>Build SQL database and add variant records, two approaches:</del>
   - <del>Either one table each per variant</del>
      - <del>could be integrated into variant_record.py</del>
    - <del>Or one table containing all positions and variations</del>
    - <del>Test parallelization with both approaches</del>
- Add metadata as SQL table
- <del>Add functionality to calculate dN/dS for whole genome in all frames</del>
- Add functionality to calculate dN/dS for a specified sampling date or date range
- Add functionality to calculate dN/dS for stop-stop-ORFs
- Add functionality to calculate dN/dS for sliding window over genome
- Add functionality to read gene locations from annotation and calculate dN/dS for them
    - Add functionality to annotate genomes
- <del>Add testing based on a set of static sequences</del>
- Clean up dNdS.py and move recurring code to nstools.py

## Research Questions
How high is the likelihood for a nonsynonymous substitution to occur at a nonsynonymous site?

How many nonsynonymous substitutions would be expected for a given sequence?

How are nonsynonymous sites spread to through a given sequence?

Do nonsynonymous substitutions spread in a similar manner to nonsynonymous sites?

Can records for multiple variants be combined to calculate dN/dS for a metadata range?

Can dN/dS of non-coding regions be used to establish baseline values for dN and dS and is this useful for error correction or significance?

Can SNP-/Variant-Calling be used to generate variant records instead of the current approach?