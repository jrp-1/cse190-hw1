# Author: Janne Rapakko (A08240805)
# denovo.py
# CSE 190 (Bandeira) - Homework 1 
import sys

print(sys.argv[1:])


def getvalue(amino):
    print(amino)
    values = {
            "A": 71,
            "R": 156,
            "N": 114,
            "D": 115,
            "C": 103,
            "E": 129,
            "Q": 128,
            "G": 57,
            "H": 137,
            "I": 113,
            "L": 113,
            "K": 128,
            "M": 131,
            "F": 147,
            "P": 97,
            "S": 87,
            "T": 101,
            "W": 186,
            "Y": 163,
            "V": 99
        }
    assert(str(amino).length == 1)
    return values[str(amino)]

# q1a & q2b BOTH DOn'T INCLUDE FULL-LENGTH/UNFRAGMENTED PEPTIDE IONS
#
# q1a Peptide match score using intensity of only main b-ions

# string input filename (\t separated), string sequence
def q1a(spectrumFile, sequence):

#####
#####   SPECTRUM READ LOGIC

#####
#####
    with open(spectrumFile) as f:
        spectrumList = f.readlines()

    spectra = {}
    for line in spectrumList:
        separator = line.find('\t')
        if (separator == -1):
            spectra[int(line[0:-1])] = 0
        else:
            spectra[int(line[0:separator])] = int(line[separator+1:-1])

    print(spectra)

#####
#####   END SPECTRUM READ LOGIC

#####
#####


# q1b peptide match score using intensity of main b -ins and y-ions


# string input filename, string sequence
def q1b(spectrumFile, sequence):
    file = open(r""+spectrumFile, "r")

# TODO: argv validations

# TODO: split spectra read in from logic

# dynamic fn call
if(sys.argv[1] == "q1a"):
    q1a(sys.argv[2], sys.argv[3])

