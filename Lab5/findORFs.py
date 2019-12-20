#!/usr/bin/env python3
# Name: Mohammad Abdulqader (mabdulqa)
# Group Members: None

from sequenceAnalysis import FastAreader, OrfFinder

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, default=100,
                                 action='store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', default=['ATG'], nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-t', '--stop', action='append', default=['TAG', 'TGA', 'TAA'], nargs='?',
                                 help='stop Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


########################################################################
# Main
# Here is the main program
#
#
########################################################################
def main(inCL=None):
    '''
    Find some genes.
    '''
    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    ###### replace the code between comments.
    lg = myCommandLine.args.longestGene  # is True if only the longest Gene is desired
    start = myCommandLine.args.start  # is a list of start codons
    minGene = myCommandLine.args.minGene  # is the minimum Gene length to include
    stop = myCommandLine.args.stop  # takes stop arguments

    orfReader = FastAreader()  # runs FastAreader
    #open(output, 'w').close()  # output file is created
    #outfile = open(output, 'a')

    for head, seq in orfReader.readFasta():  # for loop runs OrfFinder and finds the orfs.
        orfFunc = OrfFinder(seq, head, minGene, lg, start, stop)
        orfFunc.ORF(False)
        orfFunc.compStrand()
        orfFunc.outputFile()



#######


if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN


