import numpy as np

def needleman_wunsch(seq1, seq2, match_score=1, gap_penalty=-1, mismatch_penalty=-1):
    # Create the scoring matrix
    n = len(seq1) + 1
    m = len(seq2) + 1
    score_matrix = np.zeros((n, m), dtype=int)
    
    # Initialize boundary conditions (gap penalties)
    for i in range(n):
        score_matrix[i][0] = i * gap_penalty
    for j in range(m):
        score_matrix[0][j] = j * gap_penalty
    
    # Fill the scoring matrix
    for i in range(1, n):
        for j in range(1, m):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
    
    # Recursively traceback all possible optimal paths
    def traceback(i, j, aligned_seq1, aligned_seq2):
        if i == 0 and j == 0:
            # When reaching the top-left corner of the matrix, return the current path
            return [(aligned_seq1, aligned_seq2)]
        
        alignments = []
        current_score = score_matrix[i][j]
        
        # Check the current cell to decide the traceback path
        if i > 0 and current_score == score_matrix[i-1][j] + gap_penalty:
            alignments += traceback(i-1, j, seq1[i-1] + aligned_seq1, '-' + aligned_seq2)
        
        if j > 0 and current_score == score_matrix[i][j-1] + gap_penalty:
            alignments += traceback(i, j-1, '-' + aligned_seq1, seq2[j-1] + aligned_seq2)
        
        if i > 0 and j > 0:
            if current_score == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
                alignments += traceback(i-1, j-1, seq1[i-1] + aligned_seq1, seq2[j-1] + aligned_seq2)
        
        return alignments

    # Get all alignment results
    alignments = traceback(n-1, m-1, "", "")
    
    # Return the optimal alignment results and all paths
    return alignments, score_matrix[n-1][m-1]

# Input sequences dynamically
seq1 = input("Enter the first sequence: ").strip().upper()
seq2 = input("Enter the second sequence: ").strip().upper()

# Perform alignment
alignments, score = needleman_wunsch(seq1, seq2)

# Print results
print(f"Alignment Score: {score}")
for idx, (aligned_seq1, aligned_seq2) in enumerate(alignments, 1):
    print(f"Alignment {idx}:")
    print("Aligned Sequence 1:", aligned_seq1)
    print("Aligned Sequence 2:", aligned_seq2)