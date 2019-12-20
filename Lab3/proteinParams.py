#!/usr/bin/env python3
# Name: Mohammad Abdulqader (mabdulqa)
# Group Members: None

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
        nterm = (10 ** self.aaNterm) / ((10 ** self.aaNterm) + (10 ** pH))  # finds the nTerminus value
        cterm = (10 ** self.aaCterm) / ((10 ** self.aaCterm) + (10 ** pH))  # finds the cTerminus value
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


# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys


def main():
    '''
    The main function will produce the desired output; printing out all of the parameters
    mentioned earlier. It is made so the standard output can be parsed.
    '''
    inString = input('protein sequence?')
    while inString:
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print("Number of Amino Acids: {aaNum}".format(aaNum=myAAnumber))
        print("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0: myAAnumber = 1  # handles the case where no AA are present
        for key in keys:
            print("\t{} = {:.2%}".format(key, myAAcomposition[key] / myAAnumber))

        inString = input('protein sequence?')


if __name__ == "__main__":
    main()
