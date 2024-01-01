import sys
# This script was created by Tucker J Lancaster
# This script was modified by Lauren G Sabo & Karthik Krishnan (09/11-20/23)

def valid_DNA_sequence(DNA):
    DNA = DNA.upper()
    for base in DNA:
        if base != "G" and base != "C" and base != "A" and base != "T":
            return False
        
    if len(DNA)%3 != 0:
        return False
    return True

def print_DNA_sequence(DNA, mode):

    if mode == "Compact" or mode == "Verbose":
        for i in range(0,3):
            print("5' to 3' Frame: ", i)
            print(translate(DNA[i:len(DNA)], mode))
            
        rev_DNA_sequence = reverse_complement(DNA)

        for i in range(0,3):
            print("3' to 5' Frame: ", i)
            print(translate(rev_DNA_sequence[i:len(DNA)], mode))

    elif mode == "DNA":
        for i in range(0,3):
            print("5' to 3' Frame: ", i)

            broken_DNA = DNA_split(DNA[i:len(DNA)])
            broken_Codons = translate(DNA[i:len(DNA)], "DNA")
            
            for j in range(0,len(broken_DNA)):
                if j == len(broken_DNA) - 1:
                    cutoff = len(broken_DNA[j]) % 3
                    if cutoff > 0:
                        broken_DNA[j] = broken_DNA[j][:-1*cutoff]

                print(broken_DNA[j])
                print(broken_Codons[j])
            
            print("\n")
            
        rev_DNA_sequence = reverse_complement(DNA)

        for i in range(0,3):
            print("3' to 5' Frame: ", i)

            rev_broken_DNA = DNA_split(rev_DNA_sequence[i:len(rev_DNA_sequence)])
            rev_broken_Codons = translate(rev_DNA_sequence[i:len(rev_DNA_sequence)], "DNA")
            
            for j in range(0,len(rev_broken_DNA)):
                if j == len(rev_broken_DNA) - 1:
                    cutoff = len(rev_broken_DNA[j]) % 3
                    if cutoff > 0:
                        rev_broken_DNA[j] = rev_broken_DNA[j][:-1*cutoff]

                print(rev_broken_DNA[j])
                print(rev_broken_Codons[j])

            print("\n")

def DNA_split(DNA_sequence):
    i = 0
    j = 0
    if len(DNA_sequence) >= 60:
        k = len(DNA_sequence) - (len(DNA_sequence) % 60) #beginning of the last
        splitout_seq = []
        newLine = []
        while j < k:
            newLine = str(DNA_sequence[j:j+60])
            newLine = "".join(newLine)
            splitout_seq += [newLine]
            j += 60 
        
        newLine = str(DNA_sequence[k:len(DNA_sequence)])
        newLine = "".join(newLine)
        splitout_seq += [newLine]
        return splitout_seq
    
    else:
        return DNA_sequence[j:len(DNA_sequence)]
        
def translate(DNA_sequence, mode):
    out_seq = ""

    compCodonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
    'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
    }

    verbCodonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'Met',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'Stop', 'TAG':'Stop',
    'TGC':'C', 'TGT':'C', 'TGA':'Stop', 'TGG':'W',
    }

    if mode == "Compact":
        i = 0
        while i <= (len(DNA_sequence)-3):
            currCodon = str(DNA_sequence[i] + DNA_sequence[i+1] + DNA_sequence[i+2])
            out_seq += compCodonTable[currCodon]
            i += 3

    elif mode == "Verbose":
        i = 0
        while i <= (len(DNA_sequence)-3):
            currCodon = str(DNA_sequence[i] + DNA_sequence[i+1] + DNA_sequence[i+2])
            out_seq += verbCodonTable[currCodon] + " "
            i += 3

    if mode == "DNA":
        i = 0
        while i <= (len(DNA_sequence)-3):
            currCodon = str(DNA_sequence[i] + DNA_sequence[i+1] + DNA_sequence[i+2])
            out_seq += " " + compCodonTable[currCodon] + " "
            i += 3
        
        j = 0
        if len(out_seq) >= 60:
            k = len(out_seq) - (len(out_seq) % 60)
            splitout_seq = []
            newLine = []
            while j < k:
                newLine = str(out_seq[j:j+60])
                newLine = "".join(newLine)
                splitout_seq += [newLine]
                j += 60 
            
            newLine = str(out_seq[k:len(out_seq)])
            newLine = "".join(newLine)
            splitout_seq += [newLine]
            return splitout_seq
        
        else:
            return out_seq[j:len(out_seq)]

    return out_seq + "\n"

def reverse_complement(DNA_sequence):
    revCompDNASequence = list(DNA_sequence[::-1])
    for i in range(0, len(revCompDNASequence)):
        if revCompDNASequence[i] == "G":
            revCompDNASequence[i] = "C"
            continue
        if revCompDNASequence[i] == "C":
            revCompDNASequence[i] = "G"
            continue
        if revCompDNASequence[i] == "A":
            revCompDNASequence[i] = "T"
            continue
        if revCompDNASequence[i] == "T":
            revCompDNASequence[i] = "A"
            continue
    revCompDNASequence = "".join(revCompDNASequence)
    return revCompDNASequence
    

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Invalid number of inputs \nUsage: python3 Assignment2_solution.py <mode>")
        print("Mode can be one of the following options: \n COMPACT \n VERBOSE \n DNA")
        exit()
    if sys.argv[1].upper() != "COMPACT" and sys.argv[1].upper() != "VERBOSE" and sys.argv[1].upper() != "DNA":
        print(sys.argv[1].upper(), "is not a valid option \nUsage: python3 Assignment2_solution.py <mode>")
        print("Mode can be one of the following options: \n COMPACT \n VERBOSE \n DNA")
        exit()

    DNA_seq = (input("Enter DNA sequence (or Exit to quit the program): "))
    if DNA_seq.upper() == "EXIT":
        exit()
    while (valid_DNA_sequence(DNA_seq) == False):
        print("Invalid DNA sequence. Characters must be one of A, a, C, c, G, g, T, or t")
        DNA_seq = (input("Enter DNA sequence (or Exit to quit the program): "))
        if DNA_seq.upper() == "EXIT":
            exit()
    
    print_DNA_sequence(DNA_seq, sys.argv[1])
