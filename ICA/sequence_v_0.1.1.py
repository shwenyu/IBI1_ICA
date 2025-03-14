def check_begin (seq):
    begin = ["AUG"]
    leng = len(seq)
    pos = 0
    for i in range (leng-2):
        if seq[i:i+3] in begin:
            pos = i
            return (pos,leng)

def counter (seq):
    end = ["UGA","UAA","UAG"]
    count = {}
    pos, leng = check_begin(seq)
    i = pos+3
    while i+3 <= leng and not(seq[i:i+3] in end) :
        codon = seq[i:i+3]
        if not codon in count:
            count[codon] = 1
        else:
            count[codon] += 1
        i += 3
    return (count)
'''
def max_count (seq):
    count = counter(seq)
    a = count.keys(0)
    b = count[a]
    max_count_res = {a: b}
    for i in count.keys():
        j = 0
        max_val = max_count_res[max_count_res.keys(j)]
        while j+1 < len(max_count_res):
            val = max_count_res[max_count_res.keys(j+1)]
            if val > max_val:
                max_val = val
                max_pos = max_count_res.keys(j+1)
            j += 1
        if count[i] >= max_val:
            max_val = count[i]
            max_pos = i
        if not max_pos in max_count_res.keys():
            max_count_res[max_pos] = 1
        else:
            max_count_res[max_pos] += 1
    return (max_count_res) #This part has problem, the dictionary can't work, maybe change to two-dimensional list can fix!
'''
def max_count_2 (seq):
    count = counter(seq)
    key = list(count.keys())
    val = list(count.values())
    max_val = max(val)
    max_key =[]
    for i in range (len(key)):
        if val[i] == max_val:
            max_key.append(key[i])
    return(max_key, max_val)

seq = input("Inputing your sequence:")
max_key_list, max_val = max_count_2 (seq)
for i in range (len(max_key_list)):
    print(max_key_list[i], "is the maximum, with", max_val, "of its kind")
