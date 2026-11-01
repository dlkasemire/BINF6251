# Aligning Coronavirus Spike Proteins with Progressive Alignment and BLOSUM62

## Research Question
Can a pairwise alignment program I write using Needleman-Wunsch and the BLOSUM62 scoring matrix, correctly align four human coronavirus spike proteins?

## Description of the question
The spike protein sits on the surface of coronaviruses and lets the virus attach to and enter human cells. It is also the main target of the immune response and of many COVID-19 vaccines (Harvey et al., 2021). Several different coronaviruses infect humans. Some, like SARS-CoV-1 and SARS-CoV-2, cause severe disease, while others, like HCoV-OC43 and HCoV-HKU1, cause the common cold.

Lining up the spike proteins of different coronaviruses shows which parts have stayed the same (conserved) and which parts have changed. Comparing more than two sequences at once requires a multiple sequence alignment. One common way to build one is progressive alignment, which aligns the most similar sequences first and then adds the others step by step (Feng & Doolittle, 1987).
 
## Why this question matters
Conserved regions of a protein are often important for how it works, because changes there are less likely to be tolerated. A multiple sequence alignment shows at a glance which parts of the spike are shared across coronaviruses and which parts vary. In this project, I will build a simplified version of the progressive alignment process used by standard tools and test whether it gives sensible results on four real spike proteins. The project is kept to four sequences so that each step can be checked before scaling up.

## How the results may be innovative
I will build progressive alignment from the ground up, starting from my own Needleman-Wunsch pairwise aligner. Instead of using a tool as a black box, I will be able to see each step; the pairwise scores, the order in which sequences are joined, and how gaps are placed.

## Algorithm and Algorithm Class
The algorithm class is dynamic programming. The project has two steps, and the second builds on the first:

### Step 1: 
Needleman-Wunsch global pairwise alignment (Needleman & Wunsch, 1970), scored with the BLOSUM62 substitution matrix (Henikoff & Henikoff, 1992). Instead of a blanket score for every mismatch, each diagonal move looks up the score for that specific pair of amino acids in the BLOSUM62 table. Because BLOSUM62 scores range from about -4 to +11, the gap penalty also has to be chosen to fit them (see Pitfall Scan).

### Step 2: 
Progressive multiple sequence alignment, built on top of the pairwise aligner:
•	Align every pair of sequences (6 pairs for 4 sequences) and calculate how similar each pair is.
•	Start with the most similar pair and align them with Needleman-Wunsch.
•	Add the next most similar sequence to the existing alignment, and repeat until all four are aligned. When a new sequence is aligned to a group, each column of the group is scored against the new letter using the average BLOSUM62 score over the letters in that column.
•	The order in which the sequences are joined forms a simple guide tree.

This algorithm class is a good fit because, as with two sequences, comparing spike proteins letter by letter fails as soon as one has an insertion or deletion, because every letter after it shifts. Alignment finds where the gaps belong so the rest lines up.

BLOSUM62 is a better fit than a simple match/mismatch score for this question. The spike proteins of different coronaviruses differ a lot, and BLOSUM62 gives higher scores to amino acid swaps that happen often in real proteins (like I and L) and lower scores to rare ones (like W and A). A simple +1/-1 scheme treats every mismatch the same, which makes distant sequences much harder to align well.


## Data Plan
The sequences come from a UniProt file of coronavirus spike glycoproteins. The file contains 2,287 spike sequences. I will use four reviewed (Swiss-Prot) entries, since these are manually checked and complete:

Virus	Disease in humans	UniProt accession	Length (amino acids)
SARS-CoV-2	|Severe (COVID-19)	P0DTC2	1,273
SARS-CoV-1	|Severe (SARS, 2003)	P59594	1,255
HCoV-OC43	|Common cold	P36334	1,353
HCoV-HKU1 |(isolate N5)	Common cold	Q0ZME7	1,351


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

