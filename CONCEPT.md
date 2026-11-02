# Conceptual Progress Report: Aligning Coronavirus Spike Proteins with Progressive Alignment and BLOSUM62

## Recap of Project
**Research question:** Can a progressive alignment program I write, using Needleman-Wunsch and the BLOSUM62 scoring matrix, correctly align four human coronavirus spike proteins?

**Algorithm and algorithm class:** The algorithm class is dynamic programming. The project has two steps:
1. **Needleman-Wunsch global pairwise alignment**, scored with the BLOSUM62 substitution matrix and a single (linear) gap penalty.
2. **Progressive multiple sequence alignment**, which uses the pairwise aligner to compare every pair of sequences, joins the most similar pair first, and then adds the remaining sequences to the alignment one at a time.

**Data:** Four reviewed spike protein sequences from a UniProt FASTA file provided by my professor: SARS-CoV-2 (P0DTC2), SARS-CoV-1 (P59594), HCoV-OC43 (P36334) and HCoV-HKU1 (Q0ZME7). These are about 1,255 to 1,353 amino acids long. For testing, I will also use short made-up sequences and edited copies of a short piece of the SARS-CoV-2 spike.

## Inputs, Outputs, and Assumptions

### Inputs
| Input | Type and format | Example |
| --- | --- | --- |
| Sequence file | FASTA text file of protein sequences | The UniProt spike file from my professor |
| Accession numbers to use | List of text strings | `["P0DTC2", "P59594", "P36334", "Q0ZME7"]` |
| BLOSUM62 matrix | Text file in the standard NCBI matrix format (letters across the top and down the side, scores in between) | `data/BLOSUM62.txt` |
| Gap penalty | One negative number | `-10` |

### Outputs
| Output | What it contains |
| --- | --- |
| Similarity table | Percent identity and alignment score for every pair of sequences (6 pairs for 4 sequences) |
| Joining order | The order in which the sequences were added to the alignment (the simple guide tree) |
| Multiple sequence alignment | One aligned row per sequence, all the same length, with gaps shown as `-`, saved as a FASTA file |

### Assumptions
* Sequences are proteins written in one-letter amino acid codes.
* Every letter must appear in the BLOSUM62 table (the 20 standard amino acids plus B, Z, X and `*`). Lowercase letters are converted to uppercase. Any other letter stops the program with an error message.
* There are at least 2 sequences, and each sequence has a unique name.
* The alignment is **global**: every letter of every sequence is included.
* There is **one gap penalty** for every gap position (a linear gap penalty), not separate penalties for opening and extending a gap.
* When two or more moves give the same best score, the program prefers **diagonal, then up, then left**.
* When a new letter is scored against a column of the existing alignment, **gaps in that column are left out of the average**.

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
BLOSUM62 scores used: H/H = 8, A/A = 4, G/G = 6, E/A = -1, A/G = 0.

|  | start | H | E | A | G |
| --- | --- | --- | --- | --- | --- |
| **start** | 0 | -10 | -20 | -30 | -40 |
| **H** | -10 | 8 | -2 | -12 | -22 |
| **A** | -20 | -2 | 7 | 2 | -8 |
| **G** | -30 | -12 | -3 | 7 | **8** |

Tracing back from the bottom-right cell gives:
```
HEAG
H-AG
```
Score: 8 (H/H) - 10 (gap) + 4 (A/A) + 6 (G/G) = **8**. Biopython's pairwise aligner gives the same alignment and score with the same settings.

### Step 2: Progressive alignment

#### Function 6: Compare every pair
```
FUNCTION compare_all_pairs(sequences, score, gap_penalty)

    sequences is a list of (name, sequence) pairs.

    Create an empty table called "similarity".

    For every pair of different sequences (A and B):
        Align A and B using needleman_wunsch.
        Calculate their percent identity.
        Store the percent identity in the table under the pair (A, B),
        and also under (B, A), so it can be looked up in either order.

    Return the similarity table.
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

    Return the existing rows plus the new row as the updated group.
```

Existing gaps are never removed, only added to. This is where the "once a gap, always a gap" behavior of progressive alignment comes from.

#### Function 10: Put it all together
```
FUNCTION progressive_alignment(fasta_file, wanted_accessions, blosum_file, gap_penalty)

    If fewer than 2 accession numbers are given:
        Stop and print "Need at least 2 sequences to align."

    Read the sequences with read_fasta.
    Load the BLOSUM62 table with load_blosum62.
    Check every sequence with check_sequence.

    Compare all pairs with compare_all_pairs.
    Decide the order with choose_join_order.

    Align the first two sequences in the order with needleman_wunsch.
    Start the group with those two aligned rows.

    For each remaining sequence in the order:
        Add it to the group with add_to_group.

    # Sanity checks
    Check that every row in the group has the same length.
    For each row, remove the gaps and check that it matches
    the original sequence exactly.
    If either check fails, stop and print which row is wrong.

    Save the group as a FASTA file.
    Return the group, the joining order, and the similarity table.
```

#### Worked example for Step 2
Start with the group from the Step 1 example and add `HEG` (gap penalty -10):
```
HEAG
H-AG
HE-G
```
The new E lines up under column 2 (which holds E and a gap, so it scores BLOSUM62 E/E = 5), and the new sequence gets a gap under the A column. Score: 8 + 5 - 10 + 6 = **9**.

## Complexity and Bottlenecks
*To be completed.*

## Validation and Testing Plan
*To be completed.*

## Updated Pitfall and Risk Log
*To be completed.*

## Appendix: Generative AI Usage
*To be completed.*
