# DogMAG data access

This directory is the stable access index for DogMAG sequence data. It
separates source reads, final dog-wise metagenome assemblies, and dereplicated
MAGs so that deposition status is not inferred from manuscript or workflow
files.

## Current status

| Data layer | Expected records | ENA status | Index |
|---|---:|---|---|
| DogMAG raw reads | 84 runs | Public | `reads.tsv` |
| Final dog-wise primary metagenome assemblies | 41 | Public as ERZ29880033-ERZ29880073 | `assemblies.tsv` and `../accessions/primary_metagenome_assembly_accessions.tsv` |
| 95% ANI representative MAG BioSamples | 135 | Public as SAMEA123300624-SAMEA123300758 | `mags.tsv` and `../accessions/mag_sample_accessions.tsv` |
| Multi-contig representative MAG assemblies | 122 | Public with ERZ accessions | `mags.tsv` and `../accessions/mag_assembly_accessions_pass122.tsv` |
| Single-contig representative MAG assemblies | 13 | BioSamples public; FASTAs retained in article data package | `mags.tsv` |

PRJEB115259, all 84 DogMAG raw-read runs, all 41 primary assemblies, all 135
MAG BioSamples, and all 122 submitted multi-contig MAG analyses were publicly
verified on 2026-09-01. ENA Browser may return 404 for some MAG `ERS` aliases;
the corresponding `SAMEA` BioSample records are nevertheless public.

## Link policy

- ENA Browser links are stable record pages.
- ENA file-report links provide the public downloadable read and assembly files.
- No sequence file is duplicated in Git.
- A future Figshare DOI may replace the repository copies of figures and
  supplementary tables, but it must not replace the primary ENA sequence links.
