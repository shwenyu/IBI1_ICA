begin = ["AUG"]
end = ["UGA","UAA","UAG"]
seq = input("Inputing your sequence:")
count = {}
lens = len(seq)
i = 0
while i+3 <= lens:
    codon = seq[i:i+3]
    if codon in end:
        break
    if codon in begin:
        i += 3
        continue
    count[codon] += 1
    i += 3
maxco = max(count.values())
print("the most frequent cdon is", maxco)