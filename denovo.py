"""Author: Janne Rapakko (A08240805)
denovo.py
CSE 190 (Bandeira) - Homework 1"""
import sys # cmd line args

### ION OFFSETS
ION_OFFSETS = {
    "B": 1,
    "B13C": 2,
    "BNH3": -16,
    "BH20": -17,
    "A": -27,
    "ANH3": -44,
    "AH20": -45,
    "Y": 19,
    "Y13C": 20,
    "YNH3": 2,
    "YH20": 1
}
ION_LOG_GAIN = {
    "B": 2.33,
    "B13C": 0.14,
    "BNH3": 0.58,
    "BH20": 0.26,
    "A": 0.77,
    "ANH3": 0.14,
    "AH20": 0.58,
    "Y": 2.38,
    "Y13C": 0.68,
    "YNH3": 0.14,
    "YH20": 0.49
}
ION_LOG_LOSS = {  # non-negative (subtract values)
    "B": 0.85,
    "B13C": 0.02,
    "BNH3": 0.08,
    "BH20": 0.03,
    "A": 0.12,
    "ANH3": 0.02,
    "AH20": 0.08,
    "Y": 0.91,
    "Y13C": 0.1,
    "YNH3": 0.02,
    "YH20": 0.07
}
PREFIX_IONS = ["B", "B13C", "BNH3", "BH20", "A", "ANH3", "AH20"]
SUFFIX_IONS = ["Y", "Y13C", "YNH3", "YH20"]

PEPTIDE_VALUES = {
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

# amino acid weight in Da
def getvalue(amino):
    """getvalue(amino) -> int"""
    assert len(amino) == 1
    return PEPTIDE_VALUES[str(amino)]

# return dict of spectrum
def get_spectrum(spectrum_file):
    """get_spectrum(spectrum_file) -> dict """
    with open(spectrum_file, mode="r", encoding="utf-8") as file:
        spectrum_list = file.readlines()

    spectra = {}
    for line in spectrum_list:
        separator = line.find('\t')
        if separator == -1:
            spectra[int(line[0:-1])] = None
        else:
            spectra[int(line[0:separator])] = int(line[separator+1:-1])

    return spectra

def get_seq_dict_b(sequence):
    """get_seq_dict_b(sequence) -> dict
    fragment scores b-ions"""
    frag_scores = {}
    assert isinstance(sequence, str)
    for i in range(1, len(sequence) + 1):
        if frag_scores.get(sequence[:i-1]) is None: # we don't a previous value
            frag_scores[sequence[:i]] = getvalue(sequence[i-1]) + ION_OFFSETS["B"]
        else: # base case (we have no keys)
            frag_scores[sequence[:i]] = frag_scores[sequence[:i-1]] + getvalue(sequence[i-1])
    return frag_scores

def get_seq_dict_y(sequence):
    """get_seq_dict_b(sequence) -> dict
    fragment scores y-ions"""
    frag_scores = {}
    assert isinstance(sequence, str)
    y_seq = sequence[::-1]
    ## NOTE: KEYS ARE NOT IN ORDER (REVERSE OF ACTUAL y-ion fragments)
    for i in range(1, len(y_seq)):
        if frag_scores.get(y_seq[:i-1]) is None: # we don't a previous value
            frag_scores[y_seq[:i]] = getvalue(y_seq[i-1]) + ION_OFFSETS["Y"]
        else: # base case (we have no keys)
            frag_scores[y_seq[:i]] = frag_scores[y_seq[:i-1]] + getvalue(y_seq[i-1])
    return frag_scores

# q1a & q2b BOTH DOn'T INCLUDE FULL-LENGTH/UNFRAGMENTED PEPTIDE IONS

def q1a(spectrum_file, sequence):
    """q1a(spectrum_file, sequence) -> print int
    str spectrum_file, str sequence
    q1a Peptide match score using intensity of only main b-ions"""
    score = 0
    spectrum = get_spectrum(spectrum_file)
    frag_scores = get_seq_dict_b(sequence)
    for b_val in frag_scores.values():
        if spectrum.get(b_val) is None:
            continue
        score += spectrum[b_val]
    print(sequence, score)
    sys.exit()

def q1b(spectrum_file, sequence):
    """q1b(spectrum_file, sequence) -> print int
    str spectrum_file, str sequence
    q1b peptide match score using intensity of main b -ins and y-ions"""
    score = 0
    spectrum = get_spectrum(spectrum_file)
    frag_scores_b = get_seq_dict_b(sequence)
    frag_scores_y = get_seq_dict_y(sequence)
    for b_val in frag_scores_b.values():
        if spectrum.get(b_val) is None:
            continue
        score += spectrum[b_val]
    for y_val in frag_scores_y.values():
        if spectrum.get(y_val) is None:
            continue
        score += spectrum[y_val]
    print(sequence,score)
    sys.exit()

class Fragment:
    def __init__(self, p, s):
        self.prefix = p
        self.suffix = s
        self.prefix_score = 0
        self.suffix_score = 0
        for i in range(0, len(self.prefix)):
            self.prefix_score += getvalue(self.prefix[i])
        for i in range(0, len(self.suffix)):
            self.suffix_score += getvalue(self.suffix[i])

    def __str__(self):
        return "{0}.{1}:{2},{3}".format(self.prefix, self.suffix, self.prefix_score, self.suffix_score)

def get_fragments(sequence):
    """get_fragments(sequence) -> list[Fragment]"""
    fragments = []
    assert isinstance(sequence, str)
    seq_len = len(sequence)
    for i in range(1, seq_len):
        fragments.append(Fragment(sequence[:i], sequence[i:seq_len]))
    return fragments

def q2(spectrum_file, sequence):
    """q2(spectrum_file, sequence) -> print int
    str spectrum_file, sequence
    q2 peptide match score using log-likelihood scores for all a/b/y-ions"""
    log_score = 0.0
    spectrum = get_spectrum(spectrum_file)
    fragments = get_fragments(sequence)
    for fragment in fragments: # iterate through our prefixes & suffixes
        for p_ion in PREFIX_IONS: # iterate through b/a-ions
            if fragment.prefix_score + ION_OFFSETS[p_ion] in spectrum.keys():
                log_score += ION_LOG_GAIN[p_ion]
            else:
                log_score -= ION_LOG_LOSS[p_ion]
        for s_ion in SUFFIX_IONS:
            if fragment.suffix_score + ION_OFFSETS[s_ion] in spectrum.keys():
                log_score += ION_LOG_GAIN[s_ion]
            else:
                log_score -= ION_LOG_LOSS[s_ion]
    print(sequence, log_score)
    sys.exit()

def q3a(spectrum_file):
    """denovo sequence using b-ions -> prints int of sequence score"""
    sys.exit()

def q3b(spectrum_file):
    """denovo sequence using b & y-ions -> prints int of sequence score"""
    sys.exit()

def q3c(spectrum_file):
    """denovo sequence using log-likelihoods -> prints int of sequence score"""
    sys.exit()

# TODO: argv validations & better err messages
if len(sys.argv) == 1:
    print("No inputs", file=sys.stderr)
    sys.exit(1)

# dynamic fn call
if sys.argv[1] == "q1a":
    q1a(sys.argv[2], sys.argv[3])
elif sys.argv[1] == "q1b":
    q1b(sys.argv[2], sys.argv[3])
elif sys.argv[1] == "q2":
    q2(sys.argv[2], sys.argv[3])
elif sys.argv[1] == "q3a":
    q3a(sys.argv[2])
elif sys.argv[1] == "q3b":
    q3b(sys.argv[2])
elif sys.argv[1] == "q3c":
    q3c(sys.argv[2])
else:
    print("incorrect input parameters", file=sys.stderr)
    sys.exit(2)
