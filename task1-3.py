
import matplotlib.pyplot as plt

codon_to_amino_acid = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

mrna_sequence=input("enter mrna sequence:")

#task 1
def most_frequent_codon(mrna_sequence):
    
    stop_codons = {"UAA", "UAG", "UGA"}
    
    codon_counts = {}
    
    for i in range(0, len(mrna_sequence) - 2, 3):
        codon = mrna_sequence[i:i+3]
        
        if codon in stop_codons:
            break
        
        if codon in codon_counts:
            codon_counts[codon] += 1
        else:
            codon_counts[codon] = 1
    
    if not codon_counts:
        return None

    most_frequent = max(codon_counts, key=codon_counts.get)
    count = codon_counts[most_frequent]
    
    return most_frequent, count


result=most_frequent_codon(mrna_sequence)

if result:
    print(f"most frequent code is {result[0]}, time is {result[1]}")
else:
    print("none")

#task 2
def most_frequent_amino_acid(mrna_sequence):
    
    stop_codons = {"UAA", "UAG", "UGA"}
    
    amino_acid_counts = {}
    
    for i in range(0, len(mrna_sequence) - 2, 3):
        codon = mrna_sequence[i:i+3]
        
        if codon in stop_codons:
            break
        
        amino_acid = codon_to_amino_acid.get(codon)
         
        if amino_acid and amino_acid != "Stop":
           
            if amino_acid in amino_acid_counts:
                amino_acid_counts[amino_acid] += 1
            else:
                amino_acid_counts[amino_acid] = 1
    
    if not amino_acid_counts:
        return None
    
    most_frequent_amino_acid = max(amino_acid_counts, key=amino_acid_counts.get)
    count = amino_acid_counts[most_frequent_amino_acid]
    
    return most_frequent_amino_acid, count


result = most_frequent_amino_acid(mrna_sequence)

if result:
    print(f"most frequent amino acid: {result[0]}, occur {result[1]} times")
else:
    print("none")
    
#task3
def plot_amino_acid_frequencies(mrna_sequence):
    
    stop_codons = {"UAA", "UAG", "UGA"}
    
    amino_acid_counts = {}
    
    for i in range(0, len(mrna_sequence) - 2, 3):
        codon = mrna_sequence[i:i+3]
        
        if codon in stop_codons:
            break
        
        amino_acid = codon_to_amino_acid.get(codon)
        
        if amino_acid and amino_acid != "Stop":

            if amino_acid in amino_acid_counts:
                amino_acid_counts[amino_acid] += 1
            else:
                amino_acid_counts[amino_acid] = 1
    
    if not amino_acid_counts:
        print("none")
        return

    labels = amino_acid_counts.keys()
    sizes = amino_acid_counts.values()
    
    plt.figure(figsize=(8, 8))
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', pctdistance=0.8, startangle=160, colors=plt.cm.tab20.colors)
    plt.title('Frequency Distribution of Encoded Amino Acids')
    plt.axis('equal')  
    plt.show()

plot_amino_acid_frequencies(mrna_sequence)