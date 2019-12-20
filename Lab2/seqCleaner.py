#!/usr/bin/env python3
# Name: Mohammad Abdulqader (mabdulqa)
# Group Members: None

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example:
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''


class DNAstring(str):
    def length(self):
        return (len(self))

    def purify(self):  # returns an upcased and collapsed version of a sequence.
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        seq = self.upper()
        cleaner = seq.split('N')
        nNumber = seq.count("N")
        newNumber = str(nNumber)
        clean = cleaner[0] + "{" + newNumber + "}" + cleaner[-1]
        return clean


def main():  # asks for the sequence that needs to be cleaned and prints the result.
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring(data)
    pureData = thisDNA.purify()
    print(pureData)


main()