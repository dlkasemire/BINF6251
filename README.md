# BINF6251 Final Project

## Aligning Coronavirus Spike Proteins with Progressive Alignment and BLOSUM62

This project is for BINF6251. I am building a progressive multiple sequence alignment program from scratch in Python. It starts from my own Needleman-Wunsch pairwise aligner, scored with the BLOSUM62 substitution matrix, and uses it to align the spike proteins of four human coronaviruses: SARS-CoV-2, SARS-CoV-1, HCoV-OC43 and HCoV-HKU1.

The research question is whether my program can correctly align these four proteins, grouping the close relatives (the two SARS viruses, and the two common-cold viruses) together.

## Project documents
* **Part 1:** [PROPOSAL.md](PROPOSAL.md), the project proposal
* **Part 2:** [CONCEPT.md](CONCEPT.md), the conceptual progress report
* **Parts 3 and 4:** the implementation is on the `review_copy` branch. See the pull request: https://github.com/dlkasemire/BINF6251/pull/1
