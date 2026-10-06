# Finding Known Spike Mutations in SARS-CoV-2 Variants with Pairwise alignment

## Research Question
Can a pairwise alignment program I write identify known spike protein mutations in the Alpha SARS-CoV-2 variant?

## Description of the question
The spike protein is the main target of the immune response to SARS-CoV-2 and is used as the main target in many COVID-19 vaccines. Mutations in the spike protein can therefore affect how well the virus enters cells, spreads and also evades existing immunity. (Harvey et al., 2021). Each genomic variant such as Alpha or Delta, is defined by a specific combination of spike mutations (single amino acid variants and or insertions or deletions). (Mohebbi et al., 2024). 

## Why this question matters
One way to identify mutations is by comparing a viral sequence with a reference sequence and identifying the differences. My project will perform a simplified version of this process on SARS-CoV-2 spike protein sequences. In this project, I will test whether a pairwise alignment program that I write can identify these known mutations from the spike protein sequences in the Alpha SARS-CoV-2 variant.

## How the results may be innovative
I will be using Pairwise alignment, a standard algorithm as a small mutation finder. My program will turn its alignment into a list of mutations and compare that list with published mutations for each variant.

## Algorithm and Algorithm Class
The algorithm class is dynamic programming and the algorithm to be used is Needleman-Wunsch global pairwise alignment with a simple scoring scheme (+1 for a match, -1 for a mismatch and -1 for a gap).

This algorithm class is a good fit because if one is to compare two sequences letter by letter, it fails if one of them has a deletion because there is a shift in letters after that deletion making it look like a mismatch. Alignment is able to find out where the gaps belong, so the rest lines up. A variant’s spike protein is basically the same protein with just a few changes, and this is what global alignment is made for.

## Data Plan
The data reference will be the spike protein of the original SARS-CoV-2 strain (Wuhan-Hu-1), with NCBI accession number YP_009724390.1 which is 1,273 amino acids long.
I plan to select a complete Alpha variant spike protein sequence from NCBI Virus, and its accession number will be recorded in README.md file.

*The expected mutation list will be based on the published characteristic Alpha spike mutations reported by Harvey et al. (2021). The individual sequence may contain additional or missing mutations which will be reported separately.

*The data type to be used will be protein sequences in FASTA format and mutation lists as text files.

* Licensing or access considerations: The sequence data will be obtained from NCBI’s publicly accessible databases.

## Prototype Data Plan
* I plan to use very short sequences (like GAT, GT) to test the program, which I can align by hand on paper.
* I will also use a copy of the reference sequence, change a few letters, delete some that way I know exactly what the program should find.
* The prototype will use the same file format and the same reference as the real data, which means that the code that works on it should also run on the real variants.

## Success Criteria

### Expected outputs will be as follows:
* An alignment of each variant against the reference with its score
* A list of mutations found
* A table comparing the mutations the program found vs the published list

### How I will check my results
* Getting the same answer as my hand worked small examples from the program
* Finding all the changes I made in my edited copy of the reference.
* My alignment score matches the Biopython’s pairwise aligner when using the same scoring parameters.
* For the real Alpha sequence, I will also count how many published mutations my program found vs how many were missed.

## Pitfall Scan

### Data-related issues
Sequence variability.
* Why is it a concern: The published list gives the mutation that most viruses of a variant share. The single sequence I will use may have a few extra mutations or maybe missing one. A mismatch would therefore come from the data and not my code.

* What I will do: I will use a complete sequence and report both expected and any additional mutations.

### Algorithm-related issues 
Grid gets large for long sequences
* Why it’s a concern: The algorithm fills a grid with one cell for every pair of letters. Two whole virus genomes may process slowly on Python.

* What I will do: I will use only the spike protein.

### Evaluation issues 
On real data, it may be hard to tell a bug from a real difference
* Why it is a concern: If my mutation list is different from the published one, this difference could reflect an error in my code or a difference in the sequence I downloaded or in how mutations are represented. 

* What I will do: Test first on my own edited sequences where I know the right answer, compare my score with Biopython’s and then run the real variants. This will ensure that the remaining differences are from the data.

## Planned Repository Structure
```
README.md 	- what the project is and how to run it
PROPOSAL.md - this document that explains the project
aligner/ 	- my code
data/ 		- sequence files
tests/ 		- small examples for checking the code
results/ 	- output form running the code
```

## Generative AI Disclosure
Claude opus 5.5 Medium was used to understand the different algorithms and guide in selecting which one to work with. https://claude.ai/share/d67108cf-0431-4eed-9c26-bac27620ebc5

## References
Harvey, W. T., et al. (2021). SARS-CoV-2 variants, spike mutations and immune escape. Nature Reviews Microbiology, 19, 409–424. https://www.nature.com/articles/s41579-021-00573-0

Mohebbi, F., et al. (2024). Early detection of emerging viral variants through analysis of community structure of coordinated substitution networks. Nature Communications, 15, 2838. https://www.nature.com/articles/s41467-024-47304-6

