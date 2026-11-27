# Conceptual Progress Report: Aligning Coronavirus Spike Proteins with Progressive Alignment and BLOSUM62

## Recap of Project
**Research question:** Can a progressive alignment program I write, using Needleman-Wunsch and the BLOSUM62 scoring matrix, correctly align four human coronavirus spike proteins?

**Algorithm and algorithm class:** The algorithm class is dynamic programming. The project has two steps:
1. **Needleman-Wunsch global pairwise alignment**, scored with the BLOSUM62 substitution matrix and a single (linear) gap penalty.
2. **Progressive multiple sequence alignment**, which uses the pairwise aligner to compare every pair of sequences, joins the most similar pair first, and then adds the remaining sequences to the alignment one at a time. Because sequences are added one at a time, the joining order forms a chain-shaped guide tree.

**Data:** Four reviewed spike protein sequences from a UniProt FASTA file: SARS-CoV-2 (P0DTC2), SARS-CoV-1 (P59594), HCoV-OC43 (P36334) and HCoV-HKU1 (Q0ZME7). These are about 1,255 to 1,353 amino acids long. For testing, I will also use short made-up sequences and edited copies of a short piece of the SARS-CoV-2 spike.

## Changes Since the Proposal

### Change of project
My original proposal was to use Needleman-Wunsch with the simple scoring we used in class (+1 for a match, -1 for a mismatch, -1 for a gap) to find the published spike mutations of the Alpha SARS-CoV-2 variant by aligning one Alpha sequence against the Wuhan-Hu-1 reference. This stayed very close to what we had already built in class, so I changed the project to extend it in two ways. First, I replaced the +1/-1 scoring with the BLOSUM62 substitution matrix, which scores each pair of amino acids by how often they replace each other in real proteins. Second, I built progressive multiple sequence alignment on top of my pairwise aligner, so the program can align four sequences instead of two. Needleman-Wunsch is still the core of the project.

**What this change costs:** the original project had a clear answer key (six published substitutions and two deletions), while four real spike proteins have no single correct alignment. To make up for this, I test on hand-worked examples and edited copies where the answer is known, check my scores against Biopython, and check that the known close relatives group together.

### Changes since the revised proposal
* The revised proposal planned to start from the gap penalties BLAST uses with BLOSUM62 (11 to open a gap, 1 to extend it). I have switched to a single gap penalty of -10 for every gap position. This is simpler to build and lets me check every grid by hand. The cost is that long gaps are not favored over scattered ones (see Risk 1).
* The guide tree is a chain (one sequence added at a time), not a tree that joins two pairs. See Risk 6.
* The multiple sequence alignment is saved in FASTA format only. Clustal format is dropped.

## Inputs, Outputs, and Assumptions

### Inputs
| Input | Type and format | Example |
| --- | --- | --- |
| Sequence file | FASTA text file of protein sequences | The UniProt spike file |
| Accession numbers to use | List of text strings | `["P0DTC2", "P59594", "P36334", "Q0ZME7"]` |
| BLOSUM62 matrix | Text file in the standard NCBI matrix format (letters across the top and down the side, scores in between) | `data/BLOSUM62.txt` |
| Gap penalty | One negative number | `-10` |

### Outputs
| Output | What it contains |
| --- | --- |
| Similarity table | Alignment score and percent identity for every pair of sequences (6 pairs for 4 sequences) |
| Joining order | The order in which the sequences were added to the alignment. Because sequences are added one at a time, this is a chain-shaped guide tree. |
| Multiple sequence alignment | One aligned row per sequence, all the same length, with gaps shown as `-`, saved as a FASTA file |
| Conservation plot | For each column of the final alignment, how many of the four sequences share the most common letter (1 to 4), plotted along the length of the protein |

### Assumptions
* Sequences are proteins written in one-letter amino acid codes.
* Every letter must appear in the BLOSUM62 table (the 20 standard amino acids plus B, Z, X and `*`). Lowercase letters are converted to uppercase. Any other letter stops the program with an error message.
* There are at least 2 sequences, and each sequence has a unique name.
* The alignment is **global**: every letter of every sequence is included.
* There is **one gap penalty** for every gap position (a linear gap penalty), not separate penalties for opening and extending a gap. I use -10, which is larger than the worst BLOSUM62 mismatch (-4), so the program always prefers a mismatch to a gap unless the gap lets later letters line up much better. I will also try -4, -8 and -12 on the edited copies (see Validation) to see how the choice affects gap placement.
* When two or more moves give the same best score, the program prefers **diagonal, then up, then left**.
* When a new letter is scored against a column of the existing alignment, **gaps in that column are left out of the average**.
* Percent identity is always calculated by dividing by the length of the **shorter original sequence**.

## Detailed Pseudocode

The grid follows the layout used in class: the first sequence (x) goes across the top as columns, and the second sequence (y) goes down the side as rows.

### Step 1: Pairwise alignment with BLOSUM62

#### Function 1: Read the sequences
```
FUNCTION read_fasta(file_path, wanted_accessions)

    Create an empty list called "sequences".
    Open the FASTA file.

    Go through the file line by line:
        If a line starts with ">", it is a header line.
            The accession number is the text between the first and
            second "|" characters (for example, P0DTC2).
            Start a new record with this accession number.
        Otherwise, the line is part of the current record's sequence.
            Add it to the end of that sequence.

    Keep only the records whose accession number is in wanted_accessions.

    For each wanted accession number:
        If it was not found in the file:
            Stop the program and print "Accession <number> not found in file."

    Return the list of (accession number, sequence) pairs.
```

#### Function 2: Load the BLOSUM62 table
```
FUNCTION load_blosum62(file_path)

    Create an empty lookup table called "score".
    Open the BLOSUM62 file.
    Skip any comment lines (they start with "#").

    Read the first remaining line. It lists the amino acid letters
    across the top of the table (A, R, N, D, C, ...). Save these
    letters as "column_letters".

    For every remaining line in the file:
        The first item on the line is the row letter.
        The rest of the items are the scores.
        For each column letter, store its score in the lookup table
        under the pair (row letter, column letter).

    Return the lookup table.
```

#### Function 3: Check a sequence (edge cases)
```
FUNCTION check_sequence(sequence, name, score)

    Change every letter in the sequence to uppercase.
    Remove any spaces or line breaks.

    If the sequence is now empty:
        Stop the program and print "Sequence <name> is empty."

    Make a list of every letter that appears in the BLOSUM62 table.

    Go through the sequence one letter at a time:
        If a letter is not in that list:
            Stop the program and print
            "<name> has an unexpected letter <letter> at position <number>."

    Return the cleaned sequence.
```

#### Function 4: Needleman-Wunsch
```
FUNCTION needleman_wunsch(x, y, score, gap_penalty)

    x is the sequence across the top (columns).
    y is the sequence down the side (rows).
    gap_penalty is a negative number, for example -10.

    # Set up two grids
    Make a score grid with one more row than the length of y,
    and one more column than the length of x.
    Make a second grid of the same size, called the trace grid.
    The trace grid remembers which move gave each cell its score.

    # Fill the edges
    Put 0 in the top-left cell.
    Fill the top row: each cell is the cell to its left plus the gap penalty.
        Mark each of these in the trace grid as "left".
    Fill the left column: each cell is the cell above it plus the gap penalty.
        Mark each of these in the trace grid as "up".

    # Fill the rest of the grid, row by row, left to right
    For each cell:
        Diagonal score = the cell up and to the left,
                         plus the BLOSUM62 score for this row's letter
                         and this column's letter.
        Up score       = the cell above, plus the gap penalty.
        Left score     = the cell to the left, plus the gap penalty.

        Write the largest of the three scores in this cell.

        In the trace grid, record which move gave that score.
        If two or more moves tie, prefer diagonal first, then up, then left.

    # Trace back
    Start at the bottom-right cell.
    Make two empty lists for the aligned sequences.

    Repeat until you reach the top-left cell:
        Look at the trace grid for the current cell.
        If it says "diagonal":
            Add this column's letter of x and this row's letter of y.
            Move up one row and left one column.
        If it says "up":
            Add a gap "-" for x and this row's letter of y.
            Move up one row.
        If it says "left":
            Add this column's letter of x and a gap "-" for y.
            Move left one column.

    The lists were built from the end backwards, so reverse both lists.

    Return the two aligned sequences and the score in the bottom-right cell.
```

#### Function 5: Percent identity
```
FUNCTION percent_identity(aligned_x, aligned_y, x, y)

    Start a counter at 0.

    Go through the two aligned sequences one column at a time:
        If both letters are the same, and neither is a gap:
            Add 1 to the counter.

    Divide the counter by the length of the shorter original sequence.
    Multiply by 100 to get a percentage.
    Return the percentage.
```

#### Worked example for Step 1
x = `HEAG` (top), y = `HAG` (side), gap penalty = -10.
BLOSUM62 scores used: H/H = 8, A/A = 4, G/G = 6, E/A = -1, A/G = 0, H/E = 0, H/A = -2, H/G = -2, G/E = -2.

|  | start | H | E | A | G |
| --- | --- | --- | --- | --- | --- |
| **start** | 0 | -10 | -20 | -30 | -40 |
| **H** | -10 | 8 | -2 | -12 | -22 |
| **A** | -20 | -2 | 7 | 2 | -8 |
| **G** | -30 | -12 | -3 | 7 | **8** |

Tracing back from the bottom-right cell (diagonal, diagonal, left, diagonal) gives:
```
HEAG
H-AG
```
Score: 8 (H/H) - 10 (gap) + 4 (A/A) + 6 (G/G) = **8**. Biopython's pairwise aligner gives the same alignment and score with the same settings (see Validation).

### Step 2: Progressive alignment

#### Function 6: Compare every pair
```
FUNCTION compare_all_pairs(sequences, score, gap_penalty)

    sequences is a list of (name, sequence) pairs.

    Create an empty table called "similarity".
    Create an empty table called "pair_alignments".

    For every pair of different sequences (A and B):
        Align A and B using needleman_wunsch.
        Calculate their percent identity.
        Store (alignment score, percent identity) in the similarity table
        under the pair (A, B), and also under (B, A), so it can be
        looked up in either order.
        Store the two aligned rows in pair_alignments under (A, B),
        so the alignment can be reused later without aligning again.

    Return the similarity table and pair_alignments.
```

#### Function 7: Decide the joining order
```
FUNCTION choose_join_order(names, similarity)

    Find the pair with the highest percent identity.
    Start the order list with those two names.
    Put every other name in a "not yet added" list.

    While the "not yet added" list is not empty:
        For each name in the "not yet added" list:
            Find its highest similarity to ANY name already in the order.
        Pick the name with the highest of these values.
        Add it to the end of the order list.
        Remove it from the "not yet added" list.

    Return the order list.
```

Because this adds one name at a time, the result is a chain. For my data the expected order is SARS-CoV-2, SARS-CoV-1, then one of the common-cold viruses, then the other. The two common-cold viruses are never aligned with each other first (see Risk 6).

#### Function 8: Score a letter against a column of the group
```
FUNCTION column_score(column, letter, score, gap_penalty)

    column is a list of the letters in one column of the group,
    for example ["E", "-"].

    Make a list of the letters in the column, leaving out any gaps.

    If that list is empty:
        Return the gap penalty.
        (Safety rule: this should not normally happen.)

    For each letter in the list:
        Look up its BLOSUM62 score with the new letter.
    Return the average of these scores.
```

#### Function 9: Add one sequence to the group
```
FUNCTION add_to_group(group, new_sequence, score, gap_penalty)

    group is a list of aligned rows, all the same length,
    for example ["HEAG", "H-AG"].

    Split the group into columns. Column 1 is the first letter of
    every row, column 2 is the second letter of every row, and so on.

    # Set up the grids
    The columns go across the top, and the new sequence goes down the side.
    Make a score grid and a trace grid, each with one more row than the
    new sequence's length and one more column than the number of columns.

    # Fill the edges
    Same as needleman_wunsch: 0 in the top-left, then add the gap penalty
    along the top row ("left") and down the left column ("up").

    # Fill the rest of the grid
    For each cell:
        Diagonal score = the cell up and to the left, plus
                         column_score(this column, this row's letter).
        Up score       = the cell above, plus the gap penalty.
        Left score     = the cell to the left, plus the gap penalty.
        Write the largest in the cell, and record the move in the trace grid.
        Ties: prefer diagonal, then up, then left.

    # Trace back
    Start at the bottom-right cell.
    Make an empty list for each existing row, and one for the new sequence.

    Repeat until you reach the top-left cell:
        If the trace says "diagonal":
            Add this column's letters to the existing rows.
            Add this row's letter to the new sequence.
            Move up and left.
        If the trace says "up":
            Add a gap "-" to EVERY existing row.
            Add this row's letter to the new sequence.
            Move up.
        If the trace says "left":
            Add this column's letters to the existing rows.
            Add a gap "-" to the new sequence.
            Move left.

    Reverse all the lists, since they were built backwards.

    Return the existing rows plus the new row as the updated group,
    and the score in the bottom-right cell.
```

Existing gaps are never removed, only added to. This is where the "once a gap, always a gap" behavior of progressive alignment comes from.

#### Function 10: Count conservation
```
FUNCTION conservation(group)

    Make an empty list called "counts".

    For each column of the group:
        Make a list of the letters in the column, leaving out any gaps.
        If that list is empty:
            Add 0 to counts.
        Otherwise:
            Find the letter that appears most often.
            Add the number of times it appears to counts.

    Return counts.
    (Each value is from 1 to 4 for four sequences. 4 means every
    sequence has the same letter in that column.)
```

#### Function 11: Put it all together
```
FUNCTION progressive_alignment(fasta_file, wanted_accessions, blosum_file, gap_penalty)

    If fewer than 2 accession numbers are given:
        Stop and print "Need at least 2 sequences to align."

    Read the sequences with read_fasta.
    Load the BLOSUM62 table with load_blosum62.
    Check every sequence with check_sequence.

    Compare all pairs with compare_all_pairs.
    Decide the order with choose_join_order.

    Look up the first two sequences in the order in pair_alignments,
    and start the group with those two aligned rows.
    (They were already aligned by compare_all_pairs, so they are not aligned again.)

    For each remaining sequence in the order:
        Add it to the group with add_to_group.

    # Sanity checks
    Check that every row in the group has the same length.
    For each row, remove the gaps and check that it matches
    the original sequence exactly.
    If either check fails, stop and print which row is wrong.

    Count conservation with conservation(group).

    Save the group as a FASTA file.
    Save the similarity table and the joining order.
    Plot the conservation counts along the length of the alignment.
    Return the group, the joining order, the similarity table and the conservation counts.
```

#### Worked example for Step 2
Start with the group from the Step 1 example and add `HEG` (gap penalty -10). The group's columns are `[H, H]`, `[E, -]`, `[A, A]` and `[G, G]`.

Column scores used (gaps left out of the average):

| New letter | `[H,H]` | `[E,-]` | `[A,A]` | `[G,G]` |
| --- | --- | --- | --- | --- |
| **H** | 8 | 0 | -2 | -2 |
| **E** | 0 | 5 | -1 | -2 |
| **G** | -2 | -2 | 0 | 6 |

|  | start | `[H,H]` | `[E,-]` | `[A,A]` | `[G,G]` |
| --- | --- | --- | --- | --- | --- |
| **start** | 0 | -10 | -20 | -30 | -40 |
| **H** | -10 | 8 | -2 | -12 | -22 |
| **E** | -20 | -2 | 13 | 3 | -7 |
| **G** | -30 | -12 | 3 | 13 | **9** |

Tracing back from the bottom-right cell (diagonal, left, diagonal, diagonal) gives:
```
HEAG
H-AG
HE-G
```
The new E lines up under column 2 (which holds E and a gap, so it scores BLOSUM62 E/E = 5), and the new sequence gets a gap under the A column. Score: 8 + 5 - 10 + 6 = **9**.

Conservation counts for this alignment: 3, 2, 2, 3.

## Complexity and Bottlenecks

### How long it takes (time)
Most of the work is filling in the grid. Each cell takes the same small amount of work: one BLOSUM62 lookup, three additions, and picking the largest score. So the time depends on how many cells there are. For two sequences of length n and m, the grid has about n × m cells, so the time is **O(n × m)**.

For my spike proteins, which are about 1,300 amino acids long, that is about 1,300 × 1,300 = 1.7 million cells for each pairwise alignment.

For the whole project:
* Comparing every pair of my 4 sequences means 6 pairwise alignments.
* The first pair in the joining order reuses its alignment from that comparison, so it is not aligned again.
* Adding the third and fourth sequences to the group means 2 more grids of about the same size. These are a little slower, because each cell averages the new letter's score against every letter already in that column.
* So the whole program fills about 8 grids of about 1.7 million cells each.

### How much memory it uses (space)
Each alignment keeps two grids (the score grid and the trace grid), each with about n × m cells, so the memory is also **O(n × m)**. Once an alignment is finished, its grids are no longer needed, so only one alignment's grids have to be in memory at a time. The aligned rows from each pair are kept, but these are only about 1,300 letters each.

### Bottlenecks
* Filling the grids is the slowest part of the program.
* Memory grows quickly with sequence length. If both sequences are twice as long, the grid has four times as many cells. Whole coronavirus genomes are about 30,000 letters long, which would make the grid far too large. This is why I am using only the spike protein.
* More sequences means many more pairs to compare. With 4 sequences there are 6 pairs, but with 40 sequences there would be 780.

### How I will keep it manageable
* Use only the spike protein and only 4 sequences.
* Align each pair only once, save the results, and reuse the first pair's alignment to start the group.
* Delete each alignment's grids once the traceback is done.
* Store the trace grid as single characters (for example `D`, `U`, `L`) to save memory.
* Before filling a group grid, work out each column's average score against all 20 amino acids once, so each cell needs only one lookup instead of an average.
* Measure how long one real pairwise alignment takes and how much memory it uses before running the full program.

## Validation and Testing Plan

### 1. Small hand-worked examples (answer known)
| Test | Expected result |
| --- | --- |
| Align `HEAG` and `HAG` (gap penalty -10) | Alignment `HEAG` / `H-AG`, score 8 (worked out by hand above) |
| Add `HEG` to that group | Alignment `HEAG` / `H-AG` / `HE-G`, score 9 (worked out by hand above) |
| Align a sequence with an exact copy of itself | No gaps, 100% identity |
| A sequence that is empty | The program stops with the message "Sequence [name] is empty." |
| A sequence containing a letter not in BLOSUM62 (for example J) | The program stops and names the letter and its position |
| An accession number that is not in the file | The program stops with "Accession [number] not found in file." |

### 2. Edited copies of a real sequence (answer known)
I will take a short piece (about 50 letters) of the SARS-CoV-2 spike and make three copies by hand:

| Copy | Change |
| --- | --- |
| Copy 1 | Two letters swapped for different amino acids |
| Copy 2 | Three letters deleted in a row |
| Copy 3 | Both changes |

**Expected result:** in the final alignment, the swapped letters line up as mismatches, and the deleted letters appear as three gaps together in the right place. Since I made the changes myself, I know exactly where they should be.

I will run this test with gap penalties of -4, -8, -10 and -12 and record which ones keep the deletion together and in the right place.

### 3. The real data (four spike proteins)
Expected results:
* My pairwise alignment scores match Biopython's `PairwiseAligner` for all 6 pairs, using these settings: `mode="global"`, `substitution_matrix = substitution_matrices.load("BLOSUM62")`, `open_gap_score = -10` and `extend_gap_score = -10` (equal open and extend scores give a linear gap penalty).
* In the similarity table, SARS-CoV-2 is most similar to SARS-CoV-1, and HCoV-OC43 is most similar to HCoV-HKU1. Because the joining order is a chain, the OC43–HKU1 grouping is checked with the similarity table, not with the joining order.
* So the joining order starts with SARS-CoV-2 and SARS-CoV-1.
* The program finishes in a reasonable time. This is also my stress test, since these are the largest sequences the program will handle.

### Signs that something is wrong
* A pairwise score does not match Biopython's.
* The rows of the final alignment are not all the same length.
* Removing the gaps from a row does not give back the original sequence.
* A sequence aligned with an exact copy of itself gets gaps.
* SARS-CoV-2 is not grouped with SARS-CoV-1 first.
* In the similarity table, OC43 and HKU1 are not each other's closest match.
* In the edited copies, the known changes are not where I put them.

### Automated tests
Unit tests for each function, for example:
* `load_blosum62` gives -1 for (E, A) and 11 for (W, W).
* `check_sequence` stops on an empty sequence and on an unexpected letter.
* `needleman_wunsch` gives the `HEAG` / `H-AG` example.
* `column_score` gives 5 for the letter E against the column `["E", "-"]`.
* `add_to_group` gives the `HEAG` / `H-AG` / `HE-G` example.
* `conservation` gives 3, 2, 2, 3 for that example.

**End-to-end test:** run the whole program on four short made-up sequences and compare the result with my hand-worked answer.

**Checks that should always be true** (run after every alignment):
* All rows have the same length.
* Removing the gaps from each row gives back the original sequence.
* Adding up the score of each column of the returned alignment gives the same total as the returned score. For a pairwise alignment, each column is a BLOSUM62 score or a gap penalty. For a group alignment, each column of the new row is scored against the existing rows only (gaps left out of the average), or is a gap penalty.

## Updated Pitfall and Risk Log

### 1. One gap penalty may scatter gaps
**Why it matters:** With one gap penalty, a gap of 3 letters in a row costs the same as 3 separate one-letter gaps. Real insertions and deletions usually happen as a chunk, so the program may break up a deletion that should stay together. This is a result of switching from the BLAST-style open and extend penalties in the proposal to a single penalty.

**How I will handle it:** The edited copy with 3 letters deleted in a row will show whether the gap stays together, and I will test several gap penalties on it. If the gap does not stay together, I will note it as a limitation. Using separate penalties for opening and extending a gap would fix this, but it is outside the scope of this project.

### 2. How to score a column that contains a gap
**Why it matters:** When a new letter is scored against a column of the group, I had to decide what to do with gaps in that column. In the Step 2 worked example, counting the gap as -10 gives E against `[E, -]` an average of (5 - 10) / 2 = -2.5, which pushed the new E under the A column instead of under the matching E.

**How I will handle it:** `column_score` leaves gaps out of the average. I will state this choice clearly, because a different choice can give a different alignment.

### 3. Equally good alignments can look different
**Why it matters:** Sometimes two alignments have exactly the same score, especially where the same letter repeats (for example, a deletion in "YY" could go in either spot). My tie rule picks one, but Biopython might pick the other.

**How I will handle it:** When checking against Biopython, I will compare the scores, not the exact positions of the gaps.

### 4. Off-by-one mistakes in the grid
**Why it matters:** The grid has one extra row and column for the starting edge, so row i of the grid matches letter i - 1 of the sequence. Mixing these up is an easy mistake that would shift the whole alignment.

**How I will handle it:** The small hand-worked examples and the "remove the gaps and get back the original sequence" check would catch this.

### 5. Different ways to calculate percent identity
**Why it matters:** Dividing the number of matches by the shorter sequence's length gives a different number than dividing by the alignment length. For example, `HEAG` / `H-AG` has 3 matching columns. Dividing by the shorter sequence (3 letters) gives 100%, while dividing by the alignment length (4 columns) gives 75%. Using different methods would make the similarity table inconsistent.

**How I will handle it:** I will always divide by the length of the shorter original sequence, and state this in the README. This method does not count extra letters in the longer sequence against it, so I will keep that in mind when reading the table.

### 6. The joining order is a chain, not a tree of pairs
**Why it matters:** My program adds one sequence at a time. For my data, this means one of the common-cold viruses (OC43 or HKU1) is aligned against the SARS pair, where only about 30-35% of letters are identical, before it is ever aligned with its close relative. This is the least reliable step in the program, and because of "once a gap, always a gap," any gaps placed badly there stay in the final alignment. Standard tools such as Clustal Omega instead align each close pair first and then join the two groups.

**How I will handle it:** I will describe the joining order as a chain, check the OC43–HKU1 relationship with the similarity table, and treat the alignment between the SARS group and the common-cold group with extra caution when reading the conservation plot. Joining two groups together would need a group-to-group alignment function, which is outside the scope of this project. If time allows, I will add it as an extension and compare the two results.

## Appendix: Generative AI Usage
Claude Opus 5.5 (Anthropic) was used to restructure both my proposal and concept from finding known mutations in the Alpha SARS-CoV-2 variant with pairwise alignment to aligning four coronavirus spike proteins with progressive alignment and the BLOSUM62 scoring matrix. BLOSUM62 and progressive alignment were new to me, so I used Claude to help me understand them and plan how to build them as a beginner.

Prompts used include:
* I need to understand progressive alignment and how to use BLOSUM62 to align a few protein sequences. What is BLOSUM62, and how is it different from scoring +1 for a match and −1 for a mismatch?
* How does Needleman-Wunsch change when I use BLOSUM62?
* What is progressive alignment, and how does it build on pairwise alignment?
* What is a guide tree, explain to a beginner?
* How does it decide which sequences to join first?
* How do I align a new sequence to a group of sequences that are already aligned?
* How do I choose which sequences to use from a large file?
* What gap penalty should I use with BLOSUM62
* Why does it matter?
* How big can my sequences be before the program gets too slow?
* How will I know my alignment is correct?
* What mistakes am I likely to make, remember I am a beginner?
* Can we work through the pseudocode? I need to understand every part?
