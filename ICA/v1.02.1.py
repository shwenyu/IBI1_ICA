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
    #print(count) 
    return count  

def trans_ami(seq):
    count = counter(seq)
    trans_ami = {}
    codon_amino_acid = [("UUU", "Phe"),("UUC", "Phe"),("UUA", "Leu"),("UUG", "Leu"),("CUU", "Leu"),("CUC", "Leu"),("CUA", "Leu"),("CUG", "Leu"),("AUU", "Ile"),  
    ("AUC", "Ile"),("AUA", "Ile"),("AUG", "Met"),("GUU", "Val"),("GUC", "Val"),("GUA", "Val"),("GUG", "Val"),("UCU", "Ser"),("UCC", "Ser"),("UCA", "Ser"),("UCG", "Ser"),  
    ("CCU", "Pro"),("CCC", "Pro"),("CCA", "Pro"),("CCG", "Pro"),("ACU", "Thr"),("ACC", "Thr"),("ACA", "Thr"),("ACG", "Thr"),("GCU", "Ala"),("GCC", "Ala"),("GCA", "Ala"),
    ("GCG", "Ala"),("UAU", "Tyr"),("UAC", "Tyr"),("UAA", "Stop"),("UAG", "Stop"),("CAU", "His"),("CAC", "His"),("CAA", "Gln"),("CAG", "Gln"),("AAU", "Asn"),("AAC", "Asn"),  
    ("AAA", "Lys"),("AAG", "Lys"),("GAU", "Asp"),("GAC", "Asp"),("GAA", "Glu"),("GAG", "Glu"),("UGU", "Cys"),("UGC", "Cys"),("UGA", "Stop"),("UGG", "Trp"),("CGU", "Arg"),  
    ("CGC", "Arg"),("CGA", "Arg"),("CGG", "Arg"),("AGU", "Ser"),("AGC", "Ser"),("AGA", "Arg"),("AGG", "Arg"),("GGU", "Gly"),("GGC", "Gly"),("GGA", "Gly"),("GGG", "Gly")]
    for codon in count.keys():
        for codon2, amino_acid in codon_amino_acid:
            if codon == codon2 :
                if not amino_acid in trans_ami.keys():
                    trans_ami[amino_acid] = count[codon]
                else:
                    trans_ami[amino_acid] += count[codon]
    #print(trans_ami)
    return(trans_ami)

def max_count(seq):  
    count = counter(seq)
     
 
    max_codon = max(count.items(), key=lambda item: item[1])
      
    return {max_codon[0]: max_codon[1]}

def max_count_amino(seq):  
    
    ami = trans_ami(seq)  
     
    
    max_amino = max(ami.items(), key=lambda item: item[1])   
    return {max_amino[0]: max_amino[1]}

def max_ami(i):
    count_amino = {}
    codon_amino_acid = [("UUU", "Phe"),("UUC", "Phe"),("UUA", "Leu"),("UUG", "Leu"),("CUU", "Leu"),("CUC", "Leu"),("CUA", "Leu"),("CUG", "Leu"),("AUU", "Ile"),  
    ("AUC", "Ile"),("AUA", "Ile"),("AUG", "Met"),("GUU", "Val"),("GUC", "Val"),("GUA", "Val"),("GUG", "Val"),("UCU", "Ser"),("UCC", "Ser"),("UCA", "Ser"),("UCG", "Ser"),  
    ("CCU", "Pro"),("CCC", "Pro"),("CCA", "Pro"),("CCG", "Pro"),("ACU", "Thr"),("ACC", "Thr"),("ACA", "Thr"),("ACG", "Thr"),("GCU", "Ala"),("GCC", "Ala"),("GCA", "Ala"),
    ("GCG", "Ala"),("UAU", "Tyr"),("UAC", "Tyr"),("UAA", "Stop"),("UAG", "Stop"),("CAU", "His"),("CAC", "His"),("CAA", "Gln"),("CAG", "Gln"),("AAU", "Asn"),("AAC", "Asn"),  
    ("AAA", "Lys"),("AAG", "Lys"),("GAU", "Asp"),("GAC", "Asp"),("GAA", "Glu"),("GAG", "Glu"),("UGU", "Cys"),("UGC", "Cys"),("UGA", "Stop"),("UGG", "Trp"),("CGU", "Arg"),  
    ("CGC", "Arg"),("CGA", "Arg"),("CGG", "Arg"),("AGU", "Ser"),("AGC", "Ser"),("AGA", "Arg"),("AGG", "Arg"),("GGU", "Gly"),("GGC", "Gly"),("GGA", "Gly"),("GGG", "Gly")]
    for codon, amino_acid in codon_amino_acid:
        if codon == i :
            return(f"The amino acid that coden {codon} translates is {amino_acid}.")
            

def plot(seq):
    fig, axs = plt.subplots(2,1,figsize = (8,8),gridspec_kw={'height_ratios': [6, 1]})
    frequences=list(trans_ami(seq).values())
    codon_index=list(trans_ami(seq).keys())
    
    wedges, texts, autotexts = axs[0].pie(frequences, labels=codon_index, autopct='%.1f%%' ,
               startangle=160,
               textprops={'fontsize': 12, 'weight': 'bold'},
               pctdistance=0.8
            )  # Adjust distance of percentage text toward the center)
    axs[0].set_title("the frequence of each codon")
 # Adjust rotation for each percentage label
    for wedge, autotext in zip(wedges, autotexts):
    # Compute the angle of the wedge's center
        angle = (wedge.theta1 + wedge.theta2) / 2

    # Make sure the labels on the left half are rotated correctly
        if  angle > 180 :
            rotation = angle-180  # For angles > 180, rotate back to keep text upright
        else:
            rotation = 0  
        autotext.set_rotation(rotation)
    axs[1].axis('off')
    max_count_res = max_count(seq)
    max_count_ami = max_count_amino(seq)
    res = []  
    for i in max_count_res.keys():  
        #print(i, "is the maximum, with", max_count_res[i], "of its kind.")
        res.append(str(f"Codon {i} is the maximum, with {max_count_res[i]} of its kind."))
        res.append(max_ami(i))
    for i in max_count_ami.keys():  
        #print(i, "is the maximum, with", max_count_res[i], "of its kind.")
        res.append(str(f"Amino acid {i} is the maximum, with {max_count_ami[i]} of its kind."))
    

    # Add text in the lower subplot
    commentary = "\n".join(res)
    axs[1].text(0.5, 0.5, commentary, ha='center', va='center', fontsize=12, wrap=True)

    plt.get_current_fig_manager().full_screen_toggle()
    plt.tight_layout()
    plt.show()
    return()

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
        #max_count_res = max_count(seq)  
        #for i in max_count_res.keys():  
            #print(i, "is the maximum, with", max_count_res[i], "of its kind.")
        #max_ami(i)
        #print(counter(seq))
        plot(seq)
        return()
    elif "2" in stri:
        print("Welcome to Sequence Alignment(test version)!")
        seq1 = input("The first sequence is:")
        seq2 = input("The second sequence is:")
        align(seq1, seq2)
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