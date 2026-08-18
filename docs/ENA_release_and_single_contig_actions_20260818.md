# DogMAG ENA release and single-contig actions

## Public visibility audit (2026-08-18)

- `PRJEB115259` is visible through the ENA Browser API as study `ERP195477`.
- The public ENA read-run report currently returns no runs for `PRJEB115259`; the submitted FASTQ runs therefore remain unreleased or otherwise unavailable through the public report.
- The public study description is an obsolete pre-finalisation version. It refers to 81 long-read assemblies, 173.66 Gb of linked reads, and bacterial and archaeal bins. These statements do not describe the final dog-first DogMAG resource.
- Supplementary Table 9 records 122 submitted multi-contig MAG analyses and 13 single-contig representatives without assembly-analysis accessions.

## Proposed corrected ENA study title

DogMAG: dog-wise canine gut metagenome assemblies and metagenome-assembled genomes

## Proposed corrected ENA study description

DogMAG is a reusable canine gut metagenome resource integrating long- and short-read sequencing data through a dog-first assembly and genome-recovery workflow. The final BASALT input layer comprises 41 dog-wise metagenome assemblies, including 30 Flye long-read-only assemblies and 11 OPERA-MS hybrid assemblies, linked to 277 FASTQ records. BASALT selected 11,276 bin/version records from 30,556 polished or reassembled candidate versions. Explicit completeness- and contamination-based re-selection retained 3,418 medium-quality-or-better MAG candidates, including 503 high-completeness/low-contamination and 2,915 medium-only candidates. External dRep dereplication yielded 792 representatives at 99% ANI and 135 species/SGB-like representatives at 95% ANI. GTDB-Tk taxonomy is provided for the 792 strain-like representatives. The resource also contains assembly metrics, provenance metadata, viral and prophage candidate sequences, reference-panel validation outputs and reproducible workflow scripts. Previously released source studies are linked to document input-read and sample provenance.

## ENA helpdesk request for the 13 single-contig MAGs

**Subject:** Request to enable submission of 13 single-contig MAG assemblies under PRJEB115259

Dear ENA Helpdesk,

We are finalising the DogMAG canine gut metagenome resource under study PRJEB115259. We registered derived MAG samples for 135 dereplicated 95% ANI representatives. Webin-CLI successfully accepted 122 multi-contig MAG assemblies, but the remaining 13 assemblies each consist of one sequence and fail CONTIG-level validation because at least two sequences are required.

These records are metagenome-assembled genomes rather than cultured isolates. They are sufficiently complete to be retained in our operational medium-quality-or-better catalogue, but we do not wish to assert that every sequence represents a finished chromosome. The ENA MAG documentation states that submitters can request single-contig assemblies to be considered through the helpdesk when chromosome-level submission would not be appropriate.

Could you please advise whether these 13 MAG samples can be enabled for single-contig MAG submission, or provide the appropriate manifest/assembly-level route? Their aliases and accessions are listed below.

| MAG alias | ENA sample | BioSample |
|---|---|---|
| DOGMAG_ANI95_MAG027 | ERS31049474 | SAMEA123300650 |
| DOGMAG_ANI95_MAG030 | ERS31049477 | SAMEA123300653 |
| DOGMAG_ANI95_MAG032 | ERS31049479 | SAMEA123300655 |
| DOGMAG_ANI95_MAG033 | ERS31049480 | SAMEA123300656 |
| DOGMAG_ANI95_MAG044 | ERS31049491 | SAMEA123300667 |
| DOGMAG_ANI95_MAG049 | ERS31049496 | SAMEA123300672 |
| DOGMAG_ANI95_MAG050 | ERS31049497 | SAMEA123300673 |
| DOGMAG_ANI95_MAG051 | ERS31049498 | SAMEA123300674 |
| DOGMAG_ANI95_MAG066 | ERS31049513 | SAMEA123300689 |
| DOGMAG_ANI95_MAG086 | ERS31049533 | SAMEA123300709 |
| DOGMAG_ANI95_MAG090 | ERS31049537 | SAMEA123300713 |
| DOGMAG_ANI95_MAG106 | ERS31049553 | SAMEA123300729 |
| DOGMAG_ANI95_MAG114 | ERS31049561 | SAMEA123300737 |

The corresponding FASTA files and complete quality metadata are available, and we can provide validation reports or manifests on request.

Kind regards,

DogMAG submitters

## Required Webin actions

1. Replace the obsolete public study title and description with the text above.
2. Set an appropriate public release date for the study, FASTQ runs, samples and analyses, or obtain a reviewer-access mechanism before manuscript submission.
3. Send the helpdesk request above and complete the 13 single-contig submissions using the route ENA approves.
4. Refresh Supplementary Table 9 after accession assignment and regenerate the manuscript and repository checksums.
