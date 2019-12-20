#!/usr/bin/env python3
# Names: Mohammad Abdulqader (mabdulqa)
# Group Members:


import sys
import math

'''
The following code below contains 3 classes, all with the goal of converting 
the PHRED sequences of 4 main formats in fastQ files in PHRED+33 and PHRED+64
formats. The following classes in file are:
- FastQreader: which reads the FastQ file and strips assigns each line in the file
- FastQconverter: which cleans the DNA\RNA sequence and converts the PHRED quality 
    score into the desired output PHRED format. 
- CommandLine: parses the commands of the command line from STDIN and STDOUT 
    which sets the parameters for the FastQconverter class.
'''


class FastQreader:
    '''
    Reads fastq file and modifies to requested stdin or new file output
    '''

    def __init__(self, Qname=''):
        '''contructor: saves attribute fname '''
        self.Qname = Qname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.Qname is '':
            return sys.stdin
        else:
            return open(self.Qname)

    def readFastq(self):
        ''' Read an entire FastQ record and return the header-sequence'''

        with self.doOpen() as fileH:

            header = ''
            headerTwo = ''
            sequence = ''
            PHRED = ''

            # skip to first fastq header/stores as header
            line = fileH.readline()
            while not line.startswith('@'):
                line = fileH.readline()
            header = line[:].rstrip()

            # stores header/ sequence / header2 / PHRED sequence
            for line in fileH:
                if line.startswith('+'):
                    headerTwo = line[:].rstrip()
                    # Skip this line, and save the next line as score
                    line = fileH.readline()
                    PHRED += ''.join(line.rstrip().split())
                elif line.startswith('@'):
                    yield header, sequence, headerTwo, PHRED
                    header = line[:].rstrip()
                    sequence = ''
                    PHRED = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence, headerTwo, PHRED


class FastQconverter:
    '''
    The following code will have the following goals:
    - to change the PHRED sequence of the file to intended output.
    - clean the original DNA/RNA sequence given.
    - if there is a B offset, to turn all portions marked B to N in the DNA/RNA sequence.
    '''

    def __init__(self, phredIn, phredOut, sequence, Pseq, solexa, Boffset):
        '''
        Initializes class and has optional parameter Qname if one were to give a
        file in an IDE rather than the command line.
        '''
        self.Boffset = Boffset
        self.solexa = solexa
        self.inValue = phredIn
        self.outValue = phredOut
        self.sequence = sequence
        self.Pseq = Pseq
        self.differential = self.outValue - self.inValue  # differential is needed for converting PHRED seq.

    def cleanSeq(self):
        '''
        Cleans the sequence given to the function making it all AGCTUN.
        '''
        self.sequence = self.sequence.replace('*', 'N').replace('.', 'N')

    def convert(self):
        '''
        Converts the PHRED sequences to the desired Pout.
        '''
        self.newPseq = ''
        Phredseq = ''
        if self.solexa is True:
            for symbol in self.Pseq:
                Qsol = ord(symbol) - 64  # finds Qsolexa index (-5 to 62)
                Qphred = 10 * math.log10((10 ** (Qsol/10)) + 1)  # converts from Qsolexa to Qphred
                Qphred += self.differential + 64
                Pout = int(round(Qphred))
                Phredseq += chr(Pout)
        else:
            if self.Boffset is True:  # handles Boffset case.
                seqList = list(self.sequence)
                letterCount = 0
                for letter in self.Pseq:
                    if letter == 'B':
                        seqList[letterCount] = 'N'  # replaces all nucs with B phred score with N
                    letterCount += 1
                self.sequence = ''.join(seqList)
                self.Pseq.replace('B', '@')
            for symbol in self.Pseq:
                Qphred = ord(symbol) + self.differential
                Phredseq += chr(Qphred)
        self.newPseq += Phredseq

    def returnValue(self):
        '''
        Returns the cleaned sequence and the new PHRED seq.
        '''
        return self.sequence, self.newPseq


class CommandLine:
    '''
    This command line class is tailored for the FastQreader class with intention of taking in
    input PHRED quality scale and taking in an output PHRED quality scale.
    '''
    def __init__(self, inOpts=None):
        '''
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - The following program converts PHRED sequences of any given fastQ file.',
            epilog='Program epilog - Program takes in either PHRED 33, 64, Solexa, and 64B offset  fastQ files',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )

        self.parser.add_argument(
            '-P33in', '--PHRED33input', action='store', nargs='?', const=True,
            default=False, help='input is in PHRED33 format')
        self.parser.add_argument(
            '-P64in', '--PHRED64input', action='store', nargs='?', const=True,
            default=False, help='input is in PHRED64 format')
        self.parser.add_argument(
            '-P64Bin', action='store', nargs='?', const=True,
            default=False, help='input is in PHRED64B format')
        self.parser.add_argument(
            '-P64SOLin', action='store', nargs='?', const=True,
            default=False, help='input is in PHRED64SOL format')
        self.parser.add_argument(
            '-P33out', '--PHRED33output', action='store', nargs='?', const=True,
            default=False, help='output should be in PHRED33 format')
        self.parser.add_argument(
            '-P64out', '--PHRED64output', action='store', nargs='?', const=True,
            default=False, help='output should be in PHRED64 format')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


def main(inCL = None):
    '''
    Lets convert some fastQ files!
    '''
    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    inputCount = 0
    solexa = False
    Boffset = False
    if myCommandLine.args.P64SOLin:  # handles input PHRED quality
        inValue = 64
        solexa = True
        inputCount += 1
    if myCommandLine.args.P64Bin:
        inValue = 64
        Boffset = True
        inputCount += 1
    if myCommandLine.args.PHRED64input:
        inValue = 64
        inputCount += 1
    if myCommandLine.args.PHRED33input:
        inValue = 33
        inputCount += 1

    if inputCount is 0:
        print('Error: Please add input format.')
        sys.exit(1)
    if inputCount > 1:
        print('Error: You cannot have more than one input format.')
        sys.exit(1)

    if myCommandLine.args.PHRED33output:  # handles output PHRED quality
        outValue = 33
    elif myCommandLine.args.PHRED64output:
        outValue = 33
    else:
        outValue = 64

    reader = FastQreader()

    for head, seq, head2, PHRED in reader.readFastq():
        if len(seq) != len(PHRED):  # if length of sequence and PHRED different Error returns
            print('Error : The PHRED score length and the sequence length for {0} are not equal.'.format(head))
            sys.exit()
        else:   # Runs FastQconverter's functions to produce output
            converter = FastQconverter(inValue, outValue, seq, PHRED, solexa, Boffset)
            converter.cleanSeq()
            converter.convert()
            new_seq, new_PHRED = converter.returnValue()
            print('{0}\n{1}\n{2}\n{3}'.format(head, new_seq, head2, new_PHRED))


if __name__ == "__main__":
    main()

