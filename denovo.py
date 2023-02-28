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

AA_VALUES = {
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
    return AA_VALUES[str(amino)]

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
    for i in range(1, len(sequence)):
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

class Node_b:
    def __init__(self, l, a, v, i):
        self.level = l
        self.aa = a
        self.value = v
        self.intensity = i

    def __str__(self):
        return "{0}:{1}:{2}:level:{3}".format(self.aa, self.value, self.intensity, self.level)

def node_cmp_b(n):
    return n.intensity

def q3a(spectrum_file):
    """denovo sequence using b-ions -> prints int of sequence score"""
    spectrum = get_spectrum(spectrum_file)
    ### We only need the keys for b-ion based denovo
    par_mass = next(iter(spectrum))
    last_b_ion = par_mass - 18
    intensity = 0
    level = 0
    routes = []

    # STEP 1 BUILD A GRAPH while updating intensities
    # STEP 2 PICK ROUTE WITH HIHGEST INTENSITY
    if last_b_ion in spectrum.keys():
        intensity = spectrum[last_b_ion]
    if 1 not in spectrum.keys():
        spectrum[1] = 0         # WE END at H-ion

    for k,v in AA_VALUES.items(): # start our graph, only edges which are possible from end
        if last_b_ion - v in spectrum.keys():
            routes.append(Node_b(0, k, last_b_ion - v, spectrum[last_b_ion-v]))
            # find route to start
    for route in routes:
        for k,v in AA_VALUES.items():
            if (route.value - v) in spectrum.keys():
                routes.append(Node_b(route.level + 1, k + route.aa, route.value - v, route.intensity + spectrum[route.value - v]))
                if route.level + 1 > level:
                    level = route.level + 1

    final_routes = []
    for route in routes:
        if route.level == level:
            final_routes.append(route)

    final_routes.sort(key=node_cmp_b)

    final_route = final_routes.pop()

    print(final_route.aa, final_route.intensity)

    sys.exit()

class Node_by:
    def __init__(self, l, a, vb, vy, i):
        self.level = l
        self.aa = [] + a
        self.val_b = vb
        self.val_y = vy
        self.intensity = i

    def __str__(self):
        return "{0}::b:{1}::y:{2}:{3}:level:{4}".format(self.aa, self.val_b, self.val_y, self.intensity, self.level)

def q3b(spectrum_file):
    """denovo sequence using b & y-ions -> prints int of sequence score"""
    spectrum = get_spectrum(spectrum_file)
    parent_mass = next(iter(spectrum))
    last_b_ion = parent_mass - 18
    end = parent_mass - ION_OFFSETS["Y"]
    intensity = 0
    level = 0
    routes = []
    if last_b_ion in spectrum.keys():
        intensity = spectrum[last_b_ion]
    if 1 not in spectrum.keys():    ## might not need this one here
        spectrum[1] = 0         # WE END at H-ion

    spectrum[end] = 0              # OTHER END at parent mass (y-ions)

    for k,v in AA_VALUES.items(): # only possible, check both b and y-ions
        if last_b_ion - v in spectrum.keys():
            # TODO: reverse
            # initialize y-ion to 0
            routes.append(Node_by(0, [k], last_b_ion - v, 0, spectrum[last_b_ion-v]))
            # print(k, last_b_ion - v, spectrum[last_b_ion - v], parent_mass - v, spectrum[parent_mass - v])

    for route in routes:
        for k,v in AA_VALUES.items():
            if (route.val_b - v) in spectrum.keys():
                routes.append(Node_by(route.level + 1, [k] + route.aa, route.val_b - v, 0, route.intensity + spectrum[route.val_b - v]))
                if route.level + 1 > level:
                    level = route.level + 1

    final_routes = []
    for route in routes:
        if route.level == level:
            final_routes.append(route)
            print(route)

    # TODO; reverse -- y to parent_max


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
