
import matplotlib.pyplot as plt

# Dictionary mapping codons to their corresponding amino acids
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

# Input an mRNA sequence
mrna_sequence = input("enter mrna sequence:")

# Helper function to process the mRNA sequence and count codons or amino acids
def process_mrna_sequence(mrna_sequence, count_type="codon"):

    stop_codons = {"UAA", "UAG", "UGA"}  # Define stop codons
    counts = {}  # Dictionary to store counts
    
    # Iterate through the mRNA sequence in steps of 3 (codon length)
    for i in range(0, len(mrna_sequence) - 2, 3):
        codon = mrna_sequence[i:i+3]  # Extract a codon (3 bases)
        
        # Stop processing if a stop codon is encountered
        if codon in stop_codons:
            break
        
        # Determine what to count (codons or amino acids)
        if count_type == "amino_acid":
            item = codon_to_amino_acid.get(codon)  # Map codon to amino acid
            if not item or item == "Stop":  # Skip stop codons
                continue
        else:
            item = codon  # Count codons directly
        
        # Increment the count for the item or initialize it
        if item in counts:
            counts[item] += 1
        else:
            counts[item] = 1
    
    return counts

# Task 1: Find the most frequent codon
def most_frequent_codon(mrna_sequence):
    codon_counts = process_mrna_sequence(mrna_sequence, count_type="codon")
    if not codon_counts:
        return None
    most_frequent = max(codon_counts, key=codon_counts.get)
    return most_frequent, codon_counts[most_frequent]

# Task 2: Find the most frequent amino acid
def most_frequent_amino_acid(mrna_sequence):
    amino_acid_counts = process_mrna_sequence(mrna_sequence, count_type="amino_acid")
    if not amino_acid_counts:
        return None
    most_frequent = max(amino_acid_counts, key=amino_acid_counts.get)
    return most_frequent, amino_acid_counts[most_frequent]

# Task 3: Plot the frequency distribution of amino acids
def plot_amino_acid_frequencies(mrna_sequence, result1, result2):
    amino_acid_counts = process_mrna_sequence(mrna_sequence, count_type="amino_acid")
    if not amino_acid_counts:
        print("none")
        return

    # Prepare data for the pie chart
    labels = amino_acid_counts.keys()  # Amino acid labels
    sizes = amino_acid_counts.values()  # Corresponding frequencies
    
    # Plot the pie chart
    plt.figure(figsize=(12, 8), dpi=100)  # Set figure size and resolution
    plt.pie(
        sizes, labels=labels, autopct='%1.1f%%', pctdistance=0.8, 
        startangle=160, colors=plt.cm.tab20.colors, textprops={"weight": "bold"}
    )
    plt.suptitle('Frequency Distribution of Encoded Amino Acids', fontsize=16, fontweight='bold')
    plt.title(
        f'The most frequent code is {result1[0]}, occurs {result1[1]} times.\n'
        f'The most frequent amino acid: {result2[0]}, occurs {result2[1]} times.'
    )
    plt.axis('equal')  # Ensure the pie chart is circular
    plt.tight_layout()  # Adjust layout to prevent overlap
    plt.show()

# Call the functions and display results
result1 = most_frequent_codon(mrna_sequence)
if result1:
    print(f"The most frequent codon is {result1[0]}, occurs {result1[1]} times.")
else:
    print("none")

result2 = most_frequent_amino_acid(mrna_sequence)
if result2:
    print(f"The most frequent amino acid: {result2[0]}, occurs {result2[1]} times.")
else:
    print("none")

# Plot the frequency distribution of amino acids
plot_amino_acid_frequencies(mrna_sequence, result1, result2)