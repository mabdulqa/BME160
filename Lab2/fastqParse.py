#!/usr/bin/env python3
# Name: Mohammad Abdulqader (mabdulqa)
# Group Members: None

'''
Program docstring goes here
'''


class FastqString(str):
    ''' Class docstring goes here'''

    def parse(self):
        ''' Method docstring goes here'''
        editSeq = self.strip('@')
        data = editSeq.split(':')
        return data


def main(): #takes in the data input and breaks up the data into specs.
    ''' Function docstring goes here'''
    seq = input("What is the sequence?  ")
    thisSeq = FastqString(seq)
    newData = thisSeq.parse()
    print("\nInstrument = {0}\nRun ID = {1}\nFlow Cell ID = {2}".format(newData[0], newData[1], newData[2]))
    print("Flow Cell Lane = {0} \nTile Number = {1}".format(newData[3], newData[4]))
    print("X-coord = {0} \nY-coord = {1}".format(newData[5], newData[6]))


main()