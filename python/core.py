import re
from Levenshtein import distance
from Bio import SeqIO
from tqdm import tqdm
from datetime import datetime
import pprint

from utils import cycle_motif,expand_repeat

import sys

REPEAT_CLASSES = {}
repeat_track = {}
debug = False

motif_size = 6
percent_mutations = 0.2

pp = pprint.PrettyPrinter(depth=2)

start_time = datetime.now()
if debug: log = open('chr1_MUT10.log', 'w')
out = open('chr1_MUT20.out', 'w')

def repeat_class(motif):
    try: return REPEAT_CLASSES[motif]
    except KeyError:
        cycles = []
        for i in range(len(motif)):
            motif = cycle_motif(motif)
            cycles.append(motif)
        cycles.sort()
        for cycle in cycles: REPEAT_CLASSES[cycle] = cycles[0]
        return cycles[0]


# awdmuts - Allowed mutations
# wins - With insertion
# woins - Without insertion
def insertion_mutations(repeat_track):
    terminate = False

    start = repeat_track['positions'][0]
    end = repeat_track['positions'][1]
    rlen = end - start + motif_size
    mutations = repeat_track['mutations']
    
    insert = repeat_track['ins_seq'][-1]
    motif = insert.split('|')[0]
    m = motif_size
    insert_seq = insert.split('|')[1]
    s = len(insert_seq)

    awdmuts_wins = int(percent_mutations * (rlen + s))
    remain_muts = awdmuts_wins - mutations
    refrep_ul = int((s + remain_muts)//m )+ 1
    refrep_ll = int((s - remain_muts)//m)
    if refrep_ll < 0: refrep_ll = 0
    min_d = float('inf')
    for u in range(refrep_ll, refrep_ul + 1):
        repeat_seq = expand_repeat(motif, u * m)
        d = distance(repeat_seq, insert_seq)
        if (debug): print(f"Edit distance between {repeat_seq} and {insert_seq} is {d}", file=log)
        if d < min_d: min_d = d
    
    # Treating the insertion as a normal insertion if the
    # distance is greater than insertion length
    if min_d > s: min_d = s

    if min_d > remain_muts:
        if (debug): print(f"\n*** Partial extension: {remain_muts} ***")
        terminate = True
        min_d = float('inf')
        if remain_muts == 0:
            s = 0
            min_d = 0
        else:
            while min_d > remain_muts:
                insert_seq = insert_seq[:-1]
                s = len(insert_seq)
                if s == 0:
                    min_d = 0
                    break
                for u in range(s-remain_muts, s + remain_muts + 1):
                    repeat_seq = expand_repeat(motif, u)
                    d = distance(repeat_seq, insert_seq)
                    if (debug): pprint(f"Edit distance between {repeat_seq} and {insert_seq} is {d}", file=log)
                    if d < min_d: min_d = d
    
    if (debug): print(f'\nFinal edit distance is {min_d}', file=log)
    if (debug): print("", file=log)
    return [min_d, s, terminate]

def extension_mutations(repeat_track):
    terminate = False
    
    insert = repeat_track['ins_seq'][-1]
    motif = insert.split('|')[0]
    m = motif_size
    insert_seq = insert.split('|')[1]
    s = len(insert_seq)

    min_d = float('inf')
    min_d_rep = ""
    if s == m:
        if (debug): print("*** Extension Check ***", file=log)
        for u in range(1, 2*m + 1):
            repeat_seq = expand_repeat(motif, u)
            d = distance(repeat_seq, insert_seq)
            if (debug): print(f"Edit distance between {repeat_seq} and {insert_seq} is {d}", file=log)
            if d < min_d:
                min_d = d
                min_d_rep = repeat_seq

        start = repeat_track['positions'][0]
        end = repeat_track['positions'][1]
        rlen = end - start + motif_size
        mutations = repeat_track['mutations']

        awdmuts_wins = int(percent_mutations * (rlen + s))
        remain_muts = awdmuts_wins - mutations

        if min_d > remain_muts:
            if (debug): print(f"*** Partial extension: {remain_muts} ***", file=log)
            terminate = True
            min_d = float('inf')
            if remain_muts == 0:
                s = 0
                min_d = 0
            else:
                while min_d > remain_muts:
                    insert_seq = insert_seq[:-1]
                    s = len(insert_seq)
                    if s == 0:
                        min_d = 0
                        break
                    for u in range(s-remain_muts, s + remain_muts + 1):
                        repeat_seq = expand_repeat(motif, u)
                        d = distance(repeat_seq, insert_seq)
                        if (debug): print(f"Edit distance between {repeat_seq} and {insert_seq} is {d}", file=log)
                        if d < min_d: min_d = d
    
        if (debug): 
            print(f"\nFinal edit distance is {min_d}", file=log)
            print("", file=log)
    return [min_d, s, terminate, min_d_rep]


fh = open(sys.argv[1])
# out = open(sys.argv[1] + '_imperf.out', 'w')

fasta = SeqIO.parse(fh, 'fasta')

for record in fasta:
    chrom = record.id
    sequence = str(record.seq)
    iterations = len(sequence) - (motif_size-1)
    for i in tqdm(range(iterations), total=iterations):
        curr_motif = sequence[i:i+motif_size].upper()
        curr_nuc = curr_motif[-1]
        curr_rclass = repeat_class(curr_motif)
        curr_rclass_first_check = False

        if (debug): 
            print(f"\n\n ****************************** After iteration number {i+1} ****************************** \n", file=log)
            print("", file=log)
            print("************************", file=log)
            print(f"*  Nucleotide: {curr_nuc}       *", file=log)
            print(f"*  Motif     : {curr_motif}  *", file=log)
            print("************************", file=log)
            print("", file=log)
        
        if curr_rclass not in repeat_track:
            curr_rclass_first_check = True
            repeat_track[curr_rclass] = { 
                'positions': [i, i], 
                'next': curr_motif[0], 
                'prev': curr_motif, 
                'ins': 0, 
                'continue' : 1,    # checked as continuing repeat  
                'ins_seq': [], 
                'mutations': 0 
            }

        drop_rclasses = []

        for rclass in repeat_track:
            drop_rclass = False

            if debug: print(f"\n=============== {rclass} ===============", file=log)
            if (rclass == curr_rclass) and curr_rclass_first_check:
                pass
            else:
                rclass_continue = repeat_track[rclass]['continue']
                prev_motif = repeat_track[rclass]['prev']
                if repeat_track[rclass]['next'] == curr_nuc:
                    if rclass_continue:
                        if debug: print('NN: True\tRC: True', file=log)
                        repeat_track[rclass]['prev'] = cycle_motif(prev_motif)
                        repeat_track[rclass]['next'] = repeat_track[rclass]['prev'][0]
                        repeat_track[rclass]['positions'][1] = i
                    
                    else:
                        # The repeat is being continued after a previous insertion
                        # - Prior insertion repeat stayed qualified
                        # - Check for qualification of repeat after insertion
                        if debug: print('NN: True\tRC: False', file=log)

                        # Calculating the mutations in the insertion
                        insert = repeat_track[rclass]['ins_seq'][-1]
                        if debug: print("\n**** Insertion Check ****\n", file=log)
                        d = insertion_mutations(repeat_track[rclass])
                        insertion_len = d[1]
                        terminate = d[2]
                        d = d[0]

                        if terminate:
                            if debug: print("*** Repeat Termination ***", file=log)
                            start = repeat_track[rclass]['positions'][0]
                            end = repeat_track[rclass]['positions'][1] + motif_size + insertion_len
                            rlen = end - start
                            mutations = repeat_track[rclass]["mutations"] + d
                            if rlen > 12: # and (rlen * percent_mutations) >= mutations:
                                if debug: print("*** Valid Repeat ***", file=log)
                                print(chrom, start, end, rclass, rlen, mutations, sep='\t', file=out)
                            # If the repeat is stopped due to an insertion and it has to be 
                            # continued again from here
                            if rclass == curr_rclass:
                                repeat_track[rclass] = {
                                    'positions': [i, i],
                                    'next': curr_motif[0],
                                    'prev': curr_motif,
                                    'ins': 0,
                                    'continue' : 1,
                                    'ins_seq': [],
                                    'mutations': 0 
                                }
                            else:
                                drop_rclass = True
                        else:
                            repeat_track[rclass]['prev'] = cycle_motif(prev_motif)
                            repeat_track[rclass]['next'] = repeat_track[rclass]['prev'][0]
                            repeat_track[rclass]['positions'][1] = i
                            repeat_track[rclass]['mutations'] += d
                            repeat_track[rclass]['continue'] = 1

                else:
                    if rclass_continue:
                        if debug: print('NN: False\tRC: True', file=log)
                        repeat_track[rclass]['ins_seq'].append(f'{prev_motif}|{curr_nuc}')
                        repeat_track[rclass]['continue'] = 0
                    
                    else:
                        if debug: print('NN: False\tRC: False', file=log)
                        repeat_track[rclass]['ins_seq'][-1] += curr_nuc

                        d = extension_mutations(repeat_track[rclass])
                        insertion_len = d[1]
                        terminate = d[2]
                        min_d_repeat = d[3]
                        d = d[0]

                        if terminate:
                            if debug: print("*** Repeat Termination ***", file=log)
                            start = repeat_track[rclass]['positions'][0]
                            end = repeat_track[rclass]['positions'][1] + motif_size + insertion_len
                            rlen = end - start + 1
                            mutations = repeat_track[rclass]["mutations"] + d
                            if rlen > 12: # and (rlen * percent_mutations) >= mutations:
                                if debug: print("*** Valid Repeat ***", file=log)
                                print(chrom, start, end, rclass, rlen, mutations, sep='\t', file=out)
                            # If the repeat is stopped due to an insertion and it has to be 
                            # continued again from here
                            if rclass == curr_rclass:
                                repeat_track[rclass] = {
                                    'positions': [i, i],
                                    'next': curr_motif[0],
                                    'prev': curr_motif,
                                    'ins': 0,
                                    'continue' : 1,
                                    'ins_seq': [],
                                    'mutations': 0 
                                }
                            else:
                                drop_rclass = True
                        else:
                            if d != float('inf'):
                                a = len(min_d_repeat)
                                if a < motif_size: repeat_track[rclass]['prev'] = repeat_track[rclass]['prev'][a:] + repeat_track[rclass]['prev'][:a]
                                else: repeat_track[rclass]['prev'] = min_d_repeat[-1*motif_size:]
                                repeat_track[rclass]['next'] = repeat_track[rclass]['prev'][0]
                                repeat_track[rclass]['positions'][1] = i
                                repeat_track[rclass]['mutations'] += d
                                repeat_track[rclass]['continue'] = 1
            if (drop_rclass): drop_rclasses.append(rclass)
            else:
                if debug:
                    pp.pprint(repeat_track[rclass], file=log)
                    print("", file=log)
        
        for rclass in drop_rclasses:
            del repeat_track[rclass]
fh.close()
if debug: log.close()
out.close()