from itertools import product


def rev_comp(string):
    complement = string.translate(str.maketrans('ACGT', 'TGCA'))
    return complement[::-1]


def cycle_motif(motif):
    return motif[1:] + motif[0]



def num_cycles(motif, cycle):
    for i in range(1, len(motif)):
        motif = motif[1:] + motif[0]
        if motif == cycle: break
    j = len(motif) - i
    if i > j: return j
    return i


def expand_repeat(string, size):
    return_string = ''
    i = 0
    while len(return_string) < size:
        return_string += string[i]
        i += 1
        if i >= len(string):
            i = 0
    return return_string


def get_cycles(string):
    cycles = []
    for i in range(len(string)):
        cycles.append(string[i:] + string[:i])
    return cycles


def generate_repeats(motif_size):
    alphabet = ['A', 'C', 'G', 'T']
    repeat_classes = {}
    for combination in product(alphabet, repeat=motif_size):
        repeat = ''.join(combination)
        if repeat not in repeat_classes:
            cycles = get_cycles(repeat)
            rev_cycles = get_cycles(rev_comp(repeat))
            repeat_class = sorted(cycles + rev_cycles)[0]
            for motif in (cycles+rev_cycles):
                repeat_classes[motif] = repeat_class
    return repeat_classes
