##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties.


### DNA Class: INHERITS Seq class

### RNA Class:  INHERITS DNA class

### Protein Class: INHERITS Seq class


import re

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

#Make a dictionary called kyte_doolittle which has the hydrophobicity values for
# each amino acid from the Kyte-Doolittle scale: https://kelleybioinfo.org/algorithms/default.php?o=2
# For unknown ('X') amino acids, the value should be zero.

kyte_doolittle = {'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,
                  'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,
                  'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y' :-1.3}  

class Seq:

    def __init__(self, sequence, gene,species,kmers = []):
        self.sequence=str(sequence.upper().strip())
        self.gene=gene
        self.species=species
        self.kmers = kmers

    def __str__(self):
        return self.sequence

    def print_record(self):
        print(self.species + " " + self.gene + ": " + self.sequence)

    def make_kmers(self, k=3):
        for i in range(0, len(self.sequence)-(k-1)):
            self.kmers.append(self.sequence[i: i+k])

    def fasta(self):
        info = ">{} {} \n{}".format(self.species,self.gene,self.sequence)
        return info
    
    def count_nucleotides(self):
        totalA = totalT = totalC = totalG = 0
        for char in self.sequence:
            if char == 'A': totalA+=1
            if char == 'T': totalT+=1
            if char == 'C': totalC+=1
            if char == 'G': totalG+=1
        print('A:' + str(totalA) + ' T:' + str(totalT) + ' C:' + str(totalC) + 
              ' G:' + str(totalG))
            

    
    
class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence= re.sub('[^ATCGU]','N',self.sequence)
        self.geneid=geneid
 
    def analysis(self):
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

    def print_info(self):
        print(self.geneid + " ", end = '')
        super().print_record()
    
    def transcribe(self):
        return self.sequence.replace('T', 'U')
    
    def complement(self):        
        return self.sequence.translate(str.maketrans('TAGCtagc', 'ATCGATCG'))[::-1]
    

class RNA(DNA):

    def __init__(self, sequence, gene, species, geneid, codons = [], **kwargs):
        super().__init__(sequence, gene, species, geneid)
        self.sequence = self.sequence.replace('T', 'U')
        self.codons = codons
        
    def make_codons(self):
        for i in range(0, len(self.sequence), 3):
            if len(self.sequence[i:i+3]) == 3:
                self.codons.append(self.sequence[i: i+3])
 
    def translate(self):
        protein = ""
        for i in range(0, len(self.codons)):
            if self.codons[i] in standard_code:
                protein += standard_code.get(self.codons[i])
            else:
                protein += 'X'
        return protein
    
    def count_nucleotides(self):
        totalA = totalU = totalC = totalG = 0
        for char in self.sequence:
            if char == 'A': totalA+=1
            if char == 'U': totalU+=1
            if char == 'C': totalC+=1
            if char == 'G': totalG+=1
        print('A:' + str(totalA) + ' U:' + str(totalU) + ' C:' + str(totalC) + 
              ' G:' + str(totalG))
                

class Protein(Seq):

    def __init__(self, sequence, gene, species, kmers = [], aa_counts = {}, **kwargs):
        super().__init__(sequence, gene, species)
        self.sequence = re.sub('[^A-Z]','X', self.sequence)
        self.kmers = []
        self.aa_counts={'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'U':0,'V':0,'W':0,'X':0,'Y' :0}

    def tabulate_amino_acids(self):
        for i in range(0, len(self.sequence)):
            if self.sequence[i] in self.aa_counts:
                self.aa_counts[self.sequence[i]] += 1

    def total_hydro(self):
        sum = 0
        for i in range(0, len(self.sequence)):
            if self.sequence[i] in kyte_doolittle:
                sum += kyte_doolittle.get(self.sequence[i])
        return sum

    def hydro_scan(self):
        if len(self.kmers) == 0:
            self.make_kmers(5)
        hydro = []
        for i in range(0, len(self.kmers)):
            sum = 0
            for j in self.kmers[i]:
                if j in kyte_doolittle:
                    sum += kyte_doolittle.get(j)
            hydro.append(sum/len(self.kmers[i]))
        return hydro      
        







