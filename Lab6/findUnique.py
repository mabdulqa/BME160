#!/usr/bin/env python3
# Name: Mohammad Abdulqader(mabdulqa)
# Group Members: None

from sequenceAnalysis import FastAreader


class FindUnique:
    '''
    The following class has the following goals: to first to order and save the original sequences
    with thier original headers, to then create the powerset of each sequence, and to then merge and
    find the sequences that are unique to each tRna and produce them in a header and then a dotted +
    unique sequence fashion.

    The code contains:
    - an initializer where the tRNA data is stored and the original powersets are made.
    - a cleanCharacter function that cleans certain characters from the sequence.
    - a functinon to produce the powersets for each tRNA
    - a function focused on finding the unique sequences in each tRNA when compared to
      other tRNA powersets
    - and an output function that produces a stdout output.
    '''
    def __init__(self):
        '''
        __init__ produces the tRNA dictionary and the dictionary containing
        the powerset of each tRNA sequence.
        '''
        self.powerSet = []
        self.tRnaDic = {}
        self.uniqueList = []
        count = 0
        for head, seq in FastAreader().readFasta():
            self.tRnaDic[count] = [head, self.cleanChar(seq)]  # creates the dictionary with header and cleaned seq
            self.powerSet.append(self.thePowerSet(self.cleanChar(seq)))  # list containing the powerset of each seq
            count += 1

    def thePowerSet(self, sequence):
        '''
        The function produces a powerset of the given tRNA string.
        '''
        x = len(sequence)
        tRna = sequence
        thePower = set()
        for size in range(x):  # produces the powerset
            sizeOfSeq = len(sequence)
            while sizeOfSeq > size:
                thePower.add(tRna[size: sizeOfSeq])
                sizeOfSeq -= 1
        return thePower

    def cleanChar(self, sequence):
        '''
        cleanChar is a function that cleans certain characters from the tRna sequence.
        '''
        return sequence.replace('_', '').replace('-', '').replace('.', '')

    def unique(self):
        '''
        The function below makes copies of the powerSets in order to use one of the copies to put the
        sets in union in order to eliminate duplicates. It is then compared by thier original sequences
        to find the differneces in the unionized set, those differenes are then compared in order to
        remove the larger strings that contain smaller sequences that are in the unique set.
        '''
        for pset in self.powerSet:
            alike = set()  # alike is an empty set that will merge all the sequences but the interested one
            copiedPowerSet = self.powerSet.copy()
            theCopy = pset.copy()
            copiedPowerSet.remove(theCopy)  # removes the powerset of itself to prevent a crash.
            for pS in copiedPowerSet:
                alike = alike.union(pS)
            theCopy.difference_update(alike)  # alike will then be compared to the powerset of interest for uniques.

            anotherSet = theCopy.copy()  # now we need to find the essential sequences
            for tRna in theCopy:
                uniqueSet = theCopy.copy()
                uniqueSet.remove(tRna)  # pulls tRNA so that it doesnt get compared with itself.
                for tRna2 in uniqueSet:
                    if tRna in tRna2 and len(tRna) < len(tRna2):
                            anotherSet.discard(tRna2)  # Gets rid of the larger substrings.
            self.uniqueList.append(anotherSet)

    def output(self):
        '''
        The following function prints out the STD out of the STDin input.
        '''
        for sequence in self.tRnaDic:
            print(self.tRnaDic[sequence][0])  # prints the header
            print(self.tRnaDic[sequence][1])  # prints the sequence
            for position in range(len(self.tRnaDic[sequence][1])):  # for loop print the powerset elements in order
                for string in self.uniqueList[sequence]:
                    if string == self.tRnaDic[sequence][1][position:position + len(string)]:
                        print(('.' * position) + string)

######################################################################
#
# Main
# Here is the main program
#
######################################################################


def main(inCL=None):
    '''
    The main function runs FindUnique class and runs the unique function to find
    the unique sequences and then output function from the same class to produce
    the STDout output.
    '''
    tRnaClass = FindUnique()  # FindUnique cleans the sequence and produces the powerset list
    tRnaClass.unique()  # finds the uniques and essential sequences in the powerset list.
    tRnaClass.output()  # produces the output

if __name__ == "__main__":
    main()