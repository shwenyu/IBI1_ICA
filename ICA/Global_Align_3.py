import numpy as np

def needleman_wunsch(seq1, seq2, match_score=1, gap_penalty=-1, mismatch_penalty=-1, gap_extension=-0.5):
    # 创建评分矩阵
    n = len(seq1) + 1
    m = len(seq2) + 1
    score_matrix = np.zeros((n, m), dtype=int)
    
    # 初始化边界条件（gap penalties）
    score_matrix[1][0]=gap_penalty
    score_matrix[0][1]=gap_penalty
    for i in range(2,n):
        score_matrix[i][0] = score_matrix[i-1][0]+gap_extension
    for j in range(2,m):
        score_matrix[0][j] = score_matrix[0][j-1]+gap_extension
    
    # 填充评分矩阵
    for i in range(1, n):
        for j in range(1, m):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            if (score_matrix[i-1][j-1] == score_matrix[i-2][j-2]+match_score or score_matrix[i-1][j-1] == score_matrix[i-2][j-2]+mismatch_penalty) and i>2 and j>2:
                delete = score_matrix[i-1][j] + gap_penalty
                insert = score_matrix[i][j-1] + gap_penalty
            else:
                delete = score_matrix[i-1][j] + gap_extension
                insert = score_matrix[i][j-1] + gap_extension
            score_matrix[i][j] = max(match, delete, insert)
    
    # 递归回溯所有可能的最优路径
    def traceback(i, j, aligned_seq1, aligned_seq2):
        if i == 0 and j == 0:
            # 当到达矩阵的左上角时，返回当前路径
            return [(aligned_seq1, aligned_seq2)]
        
        alignments = []
        current_score = score_matrix[i][j]
        
        # 检查当前格子来决定回溯的路径
        if i > 0 and (current_score == score_matrix[i-1][j] + gap_penalty or current_score == score_matrix[i-1][j] + gap_extension ):
            alignments += traceback(i-1, j, seq1[i-1] + aligned_seq1, '-' + aligned_seq2)
        
        if j > 0 and (current_score == score_matrix[i][j-1] + gap_penalty or current_score == score_matrix[i-1][j] + gap_extension ):
            alignments += traceback(i, j-1, '-' + aligned_seq1, seq2[j-1] + aligned_seq2)
        
        if i > 0 and j > 0:
            if current_score == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
                alignments += traceback(i-1, j-1, seq1[i-1] + aligned_seq1, seq2[j-1] + aligned_seq2)
        
        return alignments

    # 获取所有比对结果
    alignments = traceback(n-1, m-1, "", "")
    
    # 返回最优比对结果及所有路径
    return alignments, score_matrix[n-1][m-1]

# 示例
seq1 = input("seq1:")
seq2 = input("seq2:")
alignments, score = needleman_wunsch(seq1, seq2)

# 打印结果
print(f"Alignment Score: {score}")
for idx, (aligned_seq1, aligned_seq2) in enumerate(alignments, 1):
    print(f"Alignment {idx}:")
    print("Aligned Sequence 1:", aligned_seq1)
    print("Aligned Sequence 2:", aligned_seq2)