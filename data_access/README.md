# DogMAG data access

This directory is the stable access index for DogMAG sequence data. It
separates source reads, final dog-wise metagenome assemblies, and dereplicated
MAGs so that deposition status is not inferred from manuscript or workflow
files.

## Current status

| Data layer | Expected records | ENA status | Index |
|---|---:|---|---|
| Source and DogMAG reads | Project-level | Linked by ENA study; PRJEB115259 currently held/private | `reads.tsv` |
| Final dog-wise primary metagenome assemblies | 41 | Accepted as ERZ29880033-ERZ29880073 | `assemblies.tsv` and `../accessions/primary_metagenome_assembly_accessions.tsv` |
| 95% ANI representative MAG BioSamples | 135 | Registered | `mags.tsv` and `../accessions/mag_sample_accessions.tsv` |
| Multi-contig representative MAG assemblies | 122 | Accepted with ERZ accessions | `mags.tsv` and `../accessions/mag_assembly_accessions_pass122.tsv` |
| Single-contig representative MAG assemblies | 13 | Pending ENA-approved submission route | `mags.tsv` |

The PRJEB115259 study is held/private. Public ENA API reports may therefore
remain empty until release; this does not invalidate the accessions recorded
in the production submission exports.

## Link policy

- ENA Browser links are stable record pages.
- ENA file-report links become populated when the corresponding records are
  publicly released.
- No sequence file is duplicated in Git.
- A future Figshare DOI may replace the repository copies of figures and
  supplementary tables, but it must not replace the primary ENA sequence links.
