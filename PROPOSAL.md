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
SARS-CoV-2	Severe (COVID-19)	P0DTC2	1,273
SARS-CoV-1	Severe (SARS, 2003)	P59594	1,255
HCoV-OC43	Common cold	P36334	1,353
HCoV-HKU1 (isolate N5)	Common cold	Q0ZME7	1,351

## Why these four: 
They form two known pairs of close relatives, SARS-CoV-2 with SARS-CoV-1, and HCoV-OC43 with HCoV-HKU1. The two pairs are only distantly related to each other. Within each pair, roughly 65-78% of amino acids are identical, while between the pairs only about 30-35% are. This gives me an expected answer to check my guide tree against, and it tests BLOSUM62 on both close and distant sequences.

•	The data type will be protein sequences in FASTA format, the BLOSUM62 matrix as a text file, and output alignments in FASTA or Clustal format.

•	Licensing or access considerations: UniProt data is freely available under the Creative Commons Attribution (CC BY 4.0) license. I will cite UniProt and record the accession numbers and download date in README.md.

## Prototype Data Plan
•	I plan to use four very short made-up protein sequences (about 5-8 letters each) that I can align by hand on paper, including the guide tree and each joining step.
•	I will check my BLOSUM62 scoring on single pairs of letters and on short pairs of sequences that I can score by hand.
•	I will take a short segment of the SARS-CoV-2 spike, make three edited copies with a few known substitutions and deletions, and check that my program lines them up the way I expect.
•	The prototype will use the same file format and the same scoring as the real data, which means that code that works on it should also run on the real spike proteins.


## Success Criteria
### What success looks like
The project will be successful if my program correctly aligns the hand-worked and edited prototype sequences, my pairwise scores match Biopython's, and it produces a complete alignment of the four spike proteins that groups the close relatives together.

### Expected outputs will be as follows:
*	A table of the pairwise alignment scores and percent identities for all 6 pairs
*	The guide tree, showing the order in which the sequences were joined
*	The final multiple sequence alignment of all four spike proteins, showing where the gaps and mismatches are
*	A simple conservation plot showing, for each position, how many of the four sequences share the same amino acid, with a brief note on any stretches that stand out as highly conserved or highly variable

### How I will check my results
* My program gives the same answers as my hand-worked small examples, including the guide tree and final alignment.
*	My pairwise alignment scores match Biopython's pairwise aligner when using the same BLOSUM62 matrix and gap penalty.
*	My similarity scores show SARS-CoV-2 closest to SARS-CoV-1, and HCoV-OC43 closest to HCoV-HKU1, matching the known relationships between these viruses.
*	If time allows, I will also compare my final alignment with one made by Clustal Omega (Sievers et al., 2011) on the same four sequences.

## Pitfall Scan

### Data-related issues
Low similarity between the two pairs
Why it is a concern: Only about 30-35% of amino acids are identical between the SARS viruses and the common-cold viruses. At that level, some of the alignment between the two groups may line up letters that are not truly related.

What I will do: I will treat conservation results across the two groups with more caution than results within each pair.

Choosing reliable sequences
Why it is a concern: Most entries in the file are unreviewed and may be partial or mislabeled.

What I will do: I will use only reviewed (Swiss-Prot) entries, check that each is a full-length spike, and record their accession numbers.


### Algorithm-related issues 
Once a gap, always a gap
Why it is a concern: In progressive alignment, gaps placed when the first sequences are joined cannot be changed later. An early mistake is carried into the final alignment (Löytynoja & Goldman, 2005; Kruspe & Stadler, 2007; Maiolo et al., 2018).

What I will do: Join the most similar sequences first, where gaps are most reliable. Adding one sequence at a time also means the last sequence added has the least influence on the gaps, so I will note the joining order when reporting results.

Choosing a biologically relevant gap penalty, and tied alignments
Why it is a concern: BLOSUM62 scores are larger than +1/-1, so the gap penalty has to be chosen to match. If it is too small, the alignment fills with gaps. If it is too large, it avoids gaps even where they belong. With BLOSUM62, a very unlikely substitution (around -4) can cost more than a gap, so the balance between the two decides where gaps go. Also, there can be more than one equally optimal alignment, so the traceback has to pick one consistently.

What I will do: Start from the gap penalties that standard tools pair with BLOSUM62 (for example, BLAST uses 11 to open a gap and 1 to extend it), test a few values on the prototype data to see which keeps known insertions and deletions in the right place, use the same value as Biopython for the score check, and use a fixed rule for breaking ties in the traceback.

Grid size and runtime
Why it is a concern: Each pairwise alignment fills a grid of about 1,300 x 1,300 cells, and I need 6 pairwise alignments plus the joining steps. This may be slow in plain Python.

What I will do: Measure the runtime on the prototype and the real data. With only four sequences, this should still be manageable


### Evaluation issues 
There is no single “true” alignment
Why it is a concern: For real spike proteins, there is no answer key, so it may be hard to tell a bug from a reasonable alignment choice.

What I will do: Test first on my hand-worked and edited sequences, where I know the right answer, and check my pairwise scores against Biopython before running the real data.

The guide tree shows similarity, not exact evolutionary history
Why it is a concern: The guide tree is built from how similar the sequences are. It is a tool for deciding the joining order, not a true family tree.

What I will do: I will describe the guide tree as a similarity grouping and only compare it with the known relationships between these viruses as a check.


## Planned Repository Structure
```
README.md 	- what the project is and how to run it
PROPOSAL.md - this document that explains the project
aligner/ 	- my code (pairwise alignment, guide tree, progressive alignment)
data/ 		- spike sequence files and the BLOSUM62 matrix
tests/ 		- small examples for checking the code
results/ 	- alignments, score tables and conservation plot
```

## Generative AI Disclosure
Claude opus 5.5 Medium was used to understand the different algorithms and guide in selecting which one to work with. https://claude.ai/share/d67108cf-0431-4eed-9c26-bac27620ebc5

## References
Feng, D.-F., & Doolittle, R. F. (1987). Progressive sequence alignment as a prerequisite to correct phylogenetic trees. Journal of Molecular Evolution, 25, 351–360.

Harvey, W. T., et al. (2021). SARS-CoV-2 variants, spike mutations and immune escape. Nature Reviews Microbiology, 19, 409–424. https://www.nature.com/articles/s41579-021-00573-0

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. Proceedings of the National Academy of Sciences, 89(22), 10915–10919.

Kruspe, M., & Stadler, P. F. (2007). Progressive multiple sequence alignments from triplets. BMC Bioinformatics, 8, 254. https://doi.org/10.1186/1471-2105-8-254

Maiolo, M., Zhang, X., Gil, M., & Anisimova, M. (2018). Progressive multiple sequence alignment with indel evolution. BMC Bioinformatics, 19, 331. https://doi.org/10.1186/s12859-018-2357-1

Needleman, S. B., & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. Journal of Molecular Biology, 48(3), 443–453.

Sievers, F., et al. (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology, 7, 539.

The UniProt Consortium. (2023). UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523–D531




