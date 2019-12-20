#!/usr/bin/env python3
# Name: Mohammad Abdulqader (mabdulqa)
# Group Members: Melody Azimi (msazimi)

from sequenceAnalysis import NucParams
from sequenceAnalysis import FastAreader

def main():
    myReader = FastAreader()  # formatted for stdin
    '''
    The variables here are renamed for more efficient calling.
    '''
    myNuc = NucParams()
    codonComp = myNuc.codonComposition()
    nCount = myNuc.nucComposition()
    aaComp = myNuc.aaComposition()
    '''
    These empty values are used as place holders if an empty fa file is added as an input.
    '''
    # sequence length place holder and a list place holder.
    fullData = []
    long = 0
    '''
    The following block takes in the input fa file, runs addSequence from 
    '''
    for head, seq in myReader.readFasta():
        myNuc.addSequence(seq)
        long += len(seq)
    for codon, aa in myNuc.rnaCodonTable.items():
        data = [codon, aa, (codonComp[codon]/aaComp[aa]), codonComp[codon]]
        fullData.append(data)
    dataList = sorted(fullData, key=lambda x: (x[1], x[0]))
    '''
    This block is the finding the seq length in Mega Base Pairs, or Mb.
    '''
    length = long
    mb = length / 1000000
    '''
    This block finds the GC content of the entire sequence. 
    '''
    gBases = nCount['G']
    cBases = nCount['C']
    total = length
    gc = ((gBases + cBases) / total)*100
    '''
    This is where the data is printed. Includes the sequence length, the GC content, 
    and the data of each relative codon in order of its amino acid single letter code.
    '''
    # These print statements print length in Mb and GC content in %.
    print('sequence length = {0:0.2f} Mb'.format(mb))
    print('\nGC content = {0:0.1f}%\n'.format(gc))
    # For loop runs the data in dataList, the sorted list of the lists called data
    for data in dataList:
        print('{:s} : {:s} {:5.1f} ({:6d})'.format(data[0], data[1], data[2] * 100, data[3]))


if __name__ == "__main__":
    main()
