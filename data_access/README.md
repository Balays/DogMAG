# DogMAG data access

This directory is the stable access index for DogMAG sequence data. It
separates source reads, final dog-wise metagenome assemblies, and dereplicated
MAGs so that deposition status is not inferred from manuscript or workflow
files.

## Current status

| Data layer | Expected records | ENA status | Index |
|---|---:|---|---|
| Source and DogMAG reads | Project-level | Source records public; newly submitted DogMAG read experiments held | `reads.tsv` |
| Final dog-wise primary metagenome assemblies | 41 | Public as ERZ29880033-ERZ29880073 | `assemblies.tsv` and `../accessions/primary_metagenome_assembly_accessions.tsv` |
| 95% ANI representative MAG BioSamples | 135 | Registered and held | `mags.tsv` and `../accessions/mag_sample_accessions.tsv` |
| Multi-contig representative MAG assemblies | 122 | Accepted with ERZ accessions and held | `mags.tsv` and `../accessions/mag_assembly_accessions_pass122.tsv` |
| Single-contig representative MAG assemblies | 13 | Pending ENA-approved submission route | `mags.tsv` |

PRJEB115259 and the 41 primary assembly records are public. Public ENA reports
for the newly submitted read experiments and MAG layers may remain empty until
those held records are released; this does not invalidate the accessions
recorded in the production submission exports.

## Link policy

- ENA Browser links are stable record pages.
- ENA file-report links become populated when the corresponding records are
  publicly released.
- No sequence file is duplicated in Git.
- A future Figshare DOI may replace the repository copies of figures and
  supplementary tables, but it must not replace the primary ENA sequence links.
