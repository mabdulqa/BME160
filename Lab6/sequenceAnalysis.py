#!/usr/bin/env python3
# Name: Mohammad Abdulqader (mabdulqa)
# Group Members: None

class NucParams:
    '''
    The class is focused on using rnaCodonTable; dnaCodonTable
    The class will give you the Amino Acid Composition, the
    codon composition, the nucleotide composition, and the count
    of each nucleotide, codon, and AA in the given sequence.
    '''
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
    dnaCodonTable = {key.replace('U', 'T'): value for key, value in rnaCodonTable.items()}

    def __init__(self, inString=''):  # optional 2nd parameter
        self.nucComp = {nuc: 0 for nuc in 'ACGTUN'}
        self.aaComp = {aa: 0 for aa in NucParams.rnaCodonTable.values()}
        self.codonComp = {codon: 0 for codon in NucParams.rnaCodonTable.keys()}
        self.addSequence(inString)


    def addSequence(self, inSeq):
        '''
        addSequence creates the dictionaries that store the count of
        codons, nucleotides, and amino acids.
        '''
        for nuc in inSeq:
            if nuc in self.nucComp:
                self.nucComp[nuc] += 1
        rnaSeq = inSeq.replace('T', 'U')  # changes the DNA sequences into RNA
        for nuc in range(0, len(rnaSeq), 3):  # assigns the codons in the sequence
            codon = rnaSeq[nuc:nuc + 3]
            if codon in self.codonComp:  # adds the codon and its corresponding aa into dictionary if real.
                self.codonComp[codon] += 1
                amino = NucParams.rnaCodonTable[codon]
                self.aaComp[amino] += 1


    def aaComposition(self):
        ''' Gives the composition of the Amino acids in the sequence. '''
        return self.aaComp

    def nucComposition(self):
        ''' nucCompostion gives the nucleotide composition of the sequence. '''
        return self.nucComp

    def codonComposition(self):
        ''' Function gives the count of the codons in the sequence. '''
        return self.codonComp

    def nucCount(self):
        ''' The count of the entire sequence is given. '''
        return sum(self.nucComp.values())


'''
The program overview: ProteinParam.py is a class file with an attached main file that runs 
on command line that should be able to take any protein sequence you give it and will give the
sequence length, molar and mass extinction, pI (isoelectric point), molar mass, and amino acid 
composition. 
'''
class ProteinParam:
    '''
    The following class ProteinParam is a class that gives the following
    paramerters for any given protein sequence:
    - amino acid count
    - amino acid composition
    - pI
    - molecular weight
    - molar and mass extinction
    The class is run by running python ProtienParam.py and then you input any given sequence.
    '''
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):
        '''
        This block creates the dictionary which contains all the aa in the sequence and how many times it is in
        the sequence.
        '''
        myProtien = str.upper(protein)  # name the input
        real = self.aa2mw.keys()  # real is a list of the keys of aa
        data = []  # empty list that will only add if the aa is in the aa2mw key list.
        for amino in myProtien:
            if amino in real:
                data.append(amino)
        space = ''.join(data).split() #seperates the sequence with spaces
        self.newData = ''.join(space).upper() #joins the sequence togehter with the spaces and caps the sequence
        self.AAtotal = 0  # initial count for total of aa, set to 0 so stuff doesnt break.
        self.compostion = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}  # dictionary which will house aa: #of aa in seq
        for aacid in real:
            self.compostion[aacid] = data.count(aacid)

    def aaCount(self):
        '''
        Gives how long the sequence is.
        '''
        for count in self.compostion.keys():
            self.AAtotal += self.compostion[count]
        return self.AAtotal

    def pI(self):
        '''
        Gives the stable pH of the sequence for it to be neutral.
        '''
        high = 14.0
        low = 0.0
        while ((high - low) > 0.01):
            mid = (high + low) / 2
            thisCharge = self._charge_(mid)
            if thisCharge < 0:
                high = mid
            else:
                low = mid
        return mid

    def aaComposition(self):
        '''
        Gives the composition of the sequence.
        '''
        return self.compostion

    def _charge_(self, pH):
        '''
        Gives the pH of the sequence to be used to find pI.
        '''
        posCharge = 0
        negCharge = 0
        nterm = (10 ** self.aaNterm) / ((10 ** self.aaNterm) + (10 ** pH))
        cterm = (10 ** self.aaCterm) / ((10 ** self.aaCterm) + (10 ** pH))
        posTerm = 0
        negTerm = 0
        for aa in self.newData:
            if aa in self.aa2chargePos.keys():  # so its done by finding the numerator and denominator seperately
                topPos = 10 ** (self.aa2chargePos[aa])
                bottomPos = (10 ** self.aa2chargePos[aa]) + (10 ** pH)
                posCharge += topPos / bottomPos  # then dividing them and adding them to posCharge
            if aa in self.aa2chargeNeg.keys():  # then repeat for negative charge
                topNeg = 10 ** pH
                bottomNeg = (10 ** self.aa2chargeNeg[aa]) + (10 ** pH)
                negCharge += topNeg / bottomNeg
            posTerm = posCharge + nterm
            negTerm = negCharge + cterm
        return posTerm - negTerm

    def molarExtinction(self):
        '''
        Gives the molar extinction of the sequence.
        '''
        molarExt = 0
        for aa in self.newData:
            if aa in self.aa2abs280.keys():
                molarExt += self.aa2abs280[aa]
        return molarExt

    def massExtinction(self):
        '''
        Gives the mass extinction of the sequence.
        '''
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        '''
        Gives molecular mass of the sequence.
        '''
        mass = 0
        if self.AAtotal != 0:
            dehydrate = (len(self.newData) - 1)
            for amino in self.newData:
                mass += self.aa2mw.get(amino)
            return mass - (dehydrate * self.mwH2O)
        else:
            return 0.0


import sys

class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class OrfFinder:
    '''
    This class will use FastAreader class to take the input fa file. The commandLine class will take the
    parameters for the OrfFinder class, from minimum gene length, what start and stop codons, or whether you
    want the longest gene in the open frame or all the possible options. It will first read through the orignal
    sequences by frame and then produce the complementary strand and then read through the inverse strand.
    The output will be sent to a text file using a list of lists that contain the frame, start position,
    end position, and the sequences length. The name of the output file is in the discretion of the user.
    '''

    def __init__(self, seq, header, minGene, largeGene, start, stop):
        '''
        __init__ defines the parameters of the OrfFinder and includes start codons, stop codons, minimum length,
        and whether the longest gene is desired or not.
        '''
        self.seq = seq
        self.header = header
        self.minGene = minGene
        if stop:
            self.stop = stop
        else:
            self.stop = ['TAA, TGA, TAG']  # sets a default stop parameter.
        if start:
            self.start = start
        else:
            self.start = ['ATG']  # sets a default start parameter.
        if minGene:
            self.minGene = minGene
        else:
            self.minGene = 100    # sets a default minimum gene length.
        if largeGene:
            self.largeGene = True  # boolean for whether only largest ORF gene wanted or not.
        else:
            self.largeGene = False
        self.dataList = []

    def ORF(self, complementary=False):
        '''
        The block of code below will go through a given sequence
        The first for loop goes through the frame and within it, it will go by codon
        to find either start or stop or codons. The first position is counted as a start, if another is found
        it gets added to the startList until a stop is found and it is all sent to dataFile function to be formatted
        for the output file.
        '''

        startList = [0]  # first element here for the dangling start case
        if complementary:
            sequence = self.reverseStrand()  # creates the reverse strand.
            orientation = -1
        else:
            sequence = self.seq
            orientation = 1
        for frame in range(3):  # goes through each frame
            oFrame = orientation * (frame + 1)  # sets the frame as pos or neg based on sequence given.
            for i in range(frame, len(sequence), 3):
                codon = ''.join(sequence[i: i+3])
                if codon in self.start:
                    startList.append(i)
                if codon in self.stop:
                    stopPos = i + 2
                    if self.largeGene and len(startList) > 0:
                        startList = [startList[0]]  # only runs largest ORF in gene
                    for start in startList:
                        self.dataFile(start, stopPos, oFrame)  # sends data to dataFile function
                    startList = []
                if i is len(sequence):  # handles the dangling end sequence
                    stopPos = i
                    for start in startList:
                        self.dataFile(start, stopPos, oFrame)
                    startList = []

    def dataFile(self, begin, end, oFrame):
        '''
        This function adds the inputs and appends it to the dataList list.
        '''
        seqLen = end - begin + 1  # finds length.
        if oFrame < 0:
            temporary = end
            end = (len(self.seq) - begin - 1)  # finds position relative to original string
            begin = (len(self.seq) - temporary - 1)
        if seqLen >= self.minGene:  # ensures that the sequence is >= minimum length.
            data = [oFrame, begin + 1, end + 1, seqLen]
            self.dataList.append(data)  # saves to main dataList

    def compStrand(self):
        '''
        This function reads the reverse strand of the sequence.
        '''
        self.ORF(complementary=True)  # runs the inverse sequence.

    def outputFile(self):
        '''
        The function outputFile sorts dataList by order of largest length and it is added into
        an output file
        '''
        sys.stdout.write('{0}\n'.format(self.header))
        for data in sorted(self.dataList, key=lambda x: (x[3], x[1]), reverse=True):   # sorts dataList
            output = '{:+d} {:>5d}..{:>5d} {:>5d}\n'.format(data[0], data[1], data[2], data[3])
            sys.stdout.write(output)

    def clearData(self):
        '''
        The function clears the dataList once the output is already made.
        '''
        self.dataList = []  # clears dataList if needed

    def reverseStrand(self):
        '''
        Creates the reverse strand of the sequence given.
        '''
        reverse = self.seq.translate(str.maketrans('AGCTUN', 'TCGAAN'))[::-1]   # creates reverse sequence.
        return reverse
