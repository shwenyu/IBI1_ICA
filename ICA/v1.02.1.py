import matplotlib.pyplot as plt
import numpy as np
def check_begin(seq):  
    begin = ["AUG"]  
    leng = len(seq)  
    for i in range(leng-2):  
        if seq[i:i+3] in begin:  
            return (i, leng)  
    return (None,leng)  

def counter(seq):  
    end = ["UGA", "UAA", "UAG"]  
    count = {}  
    pos, leng = check_begin(seq)  
    if pos is None:  
        return {}  
    i = pos + 3  
    while i + 3 <= leng and seq[i:i+3] not in end:  
        codon = seq[i:i+3]  
        count[codon] = count.get(codon, 0) + 1  
        i += 3  
    return count  

def max_count(seq):  
    count = counter(seq)  
    if not count:  
        return {}  
    max_codon = max(count.items(), key=lambda item: item[1])  
    return {max_codon[0]: max_codon[1]}

def max_ami(i):
    codon_amino_acid = [("UUU", "Phe"),("UUC", "Phe"),("UUA", "Leu"),("UUG", "Leu"),("CUU", "Leu"),("CUC", "Leu"),("CUA", "Leu"),("CUG", "Leu"),("AUU", "Ile"),  
    ("AUC", "Ile"),("AUA", "Ile"),("AUG", "Met"),("GUU", "Val"),("GUC", "Val"),("GUA", "Val"),("GUG", "Val"),("UCU", "Ser"),("UCC", "Ser"),("UCA", "Ser"),("UCG", "Ser"),  
    ("CCU", "Pro"),("CCC", "Pro"),("CCA", "Pro"),("CCG", "Pro"),("ACU", "Thr"),("ACC", "Thr"),("ACA", "Thr"),("ACG", "Thr"),("GCU", "Ala"),("GCC", "Ala"),("GCA", "Ala"),
    ("GCG", "Ala"),("UAU", "Tyr"),("UAC", "Tyr"),("UAA", "Stop"),("UAG", "Stop"),("CAU", "His"),("CAC", "His"),("CAA", "Gln"),("CAG", "Gln"),("AAU", "Asn"),("AAC", "Asn"),  
    ("AAA", "Lys"),("AAG", "Lys"),("GAU", "Asp"),("GAC", "Asp"),("GAA", "Glu"),("GAG", "Glu"),("UGU", "Cys"),("UGC", "Cys"),("UGA", "Stop"),("UGG", "Trp"),("CGU", "Arg"),  
    ("CGC", "Arg"),("CGA", "Arg"),("CGG", "Arg"),("AGU", "Ser"),("AGC", "Ser"),("AGA", "Arg"),("AGG", "Arg"),("GGU", "Gly"),("GGC", "Gly"),("GGA", "Gly"),("GGG", "Gly")]
    for codon, amino_acid in codon_amino_acid:
        if codon == i :
            print (f"the amino acid that coden {codon} translates is {amino_acid}.")
            return 

def plot(seq):
    fig, axs = plt.subplots(1,2,figsize = (10,5))
    frequences=list(counter(seq).values())
    codon_index=list(counter(seq).keys())
    
    axs[0].pie(frequences, labels=codon_index, autopct='%.1f%%' ,startangle=140)
    axs[0].set_title("the frequence of each codon")
    
    
    
    axs[1].bar(codon_index, frequences)
    axs[1].set_title("the frequence of each codon ")
    plt.xticks(codon_index, rotation = 90)
    plt.xlabel("codens")
    plt.ylabel("frequences(%)")
    plt.tight_layout()
    plt.show()

def align(seq1, seq2, match_score=1, gap_penalty=-1, mismatch_penalty=-1):
    # 创建评分矩阵
    n = len(seq1) + 1
    m = len(seq2) + 1
    score_matrix = np.zeros((n, m), dtype=int)
    
    # 初始化边界条件（gap penalties）
    for i in range(n):
        score_matrix[i][0] = i * gap_penalty
    for j in range(m):
        score_matrix[0][j] = j * gap_penalty
    
    # 填充评分矩阵
    for i in range(1, n):
        for j in range(1, m):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
    
    # 回溯找最优比对路径
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = n - 1, m - 1
    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        if current_score == score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif current_score == score_matrix[i-1][j] + gap_penalty:
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            j -= 1
    
    # 处理剩余部分
    while i > 0:
        aligned_seq1.append(seq1[i-1])
        aligned_seq2.append('-')
        i -= 1
    while j > 0:
        aligned_seq1.append('-')
        aligned_seq2.append(seq2[j-1])
        j -= 1
    
    # 返回比对结果
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), score_matrix[n-1][m-1]



def mode_choose(stri):
    if "1" in stri:
        print("Welcome to codon counter!")
        seq = input("Inputting your sequence: ")  
        max_count_res = max_count(seq)  
        for i in max_count_res.keys():  
            print(i, "is the maximum, with", max_count_res[i], "of its kind.")
        max_ami(i)
        #print(counter(seq))
        plot(seq)
        return()
    elif "2" in stri:
        print("Welcome to Sequence Alignment(test version)!")
        seq1 = input("The first sequence is:")
        seq2 = input("The second sequence is:")
        algin(seq1, seq2)
        al_seq1, al_seq2, score = align(seq1, seq2)
        print("Aligned Sequence 1:", al_seq1)
        print("Aligned Sequence 2:", al_seq2)
        print("Alignment Score:", score)
        return()

print("Please choose your mode: 1, Codon counter(default); 2, Sequence alignment(test). What's the mode you choose:", end="")
stri = input()
mode_choose(stri)


'''e.g to input:AUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUGGUGGCCGAGCUGCGGAAGUACUUCAAGAUCGUCGGCAUCAUGCGGCGCAAGAAGGCACCCGGGUG

output "AAG is the maximum, with 20 of its kind.
the amino acid that coden AAG translates is Lys.
                pie plot                       
                bar plot         "



'''