#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "utils.h"

#include "levenshtein.h"
using namespace std;
using namespace utils;

// global variable
unordered_map<string, string> rClassMap;
unordered_map<string, utils::RepeatTracker> globalRepeatTracker;
const uint motif_size = 6;
const float fraction_mutations = 0.1;
uint debug = 1;

vector<uint> insertion_mutations(utils::RepeatTracker rtracker) {

    uint terminate = 0;

    uint start = rtracker.start;
    uint end   = rtracker.end;
    uint rlen  = end - start;
    uint muts  = rtracker.mutations;

    string valid_motif = rtracker.valid_motif;
    string insert = rtracker.insert;
    uint s = insert.length();
    uint m = motif_size;
    uint min_d = -1;

    uint threshold_muts = fraction_mutations * (rlen + s);
    uint remain_muts  = threshold_muts - muts;

    if (debug) {
        cout << "\n*** Insertion " << insert << " ***\n";
        cout << "*** Remaining mutations: " << remain_muts << " ***\n";
    }

    if (remain_muts == 0) {
        if (debug) { cout << "*** Zero Remaining mutations! ***\n*** Termination ***\n\n" ; }
        s = 0; min_d = 0; terminate = 1;
        return vector<uint> { min_d, s, terminate };
    }

    uint refrep_ul = ((s + remain_muts)/m) + 1;
    uint refrep_ll = 0;
    if (remain_muts < s) { refrep_ll = ((s - remain_muts)/m); }
    if (debug) { cout << "\n*** Unit(s) Insertions ***\n"; }
    for (uint u = refrep_ll; u <= refrep_ul; u++) {
        string repeat_seq = utils::expand_repeat(valid_motif, u*m);
        int x = repeat_seq.length();
        int y = insert.length();
        int* d_matrix = levenshtein_matrix(repeat_seq.length(), repeat_seq, insert.length(), insert);
        uint d = (&d_matrix[0])[((x+1)*(y+1))-1];
        if (debug) { cout << "Comparing:  " << insert << "\t" << repeat_seq << "\tMutations: " << d << "\n"; }
        if (debug) {
            cout << "*** Levenshtein Matrix ***\n";
            for (int j=0; j<((x+1)*(y+1)); j++) {
                if (repeat_seq.length()>0 && j > 0) {
                    if (j % (repeat_seq.length()+1) == 0) { cout << "\n"; }
                }
                cout << (&d_matrix[0])[j] << " ";
            }
            cout << "\n";
        }
        if (d < min_d) { min_d = d; }
        int nexpmuts = ((u+1)*m) - s;
        if (min_d < nexpmuts) {
            if (debug) { cout << "*** Next expected mutations exceed minimum ***\n"; }
            break;
        }
    }
    // If mutations greater than size of insert
    // Treating the sequence as actual insert
    if (min_d > s) { min_d = s; }

    // If minimum mutations are greater than allowed remaining mutations
    if (min_d > remain_muts) {
        if (debug) { cout << "\n*** Least mutations greater than allowed ***\n" ; }
        terminate = 1;
        uint max_rlen = end - start + s;
        // cout << "Maximum Repeat length: " << max_rlen << "\n";
        if (max_rlen < 12) {
            s = 0; min_d = 0;
            if (debug) { cout << "*** Short repeat length ***" ; }
        }
        else {
            min_d = -1;
            insert += 'N';
            while (min_d > remain_muts) {
                insert = insert.substr(0, insert.length()-1);
                s = insert.length();
                if (s == 0) { min_d = 0; break; }
                for (uint u = s-remain_muts; u <= s+remain_muts; u++) {
                    string repeat_seq = utils::expand_repeat(valid_motif, u);
                    uint d = levenshtein_distance(repeat_seq.length(), repeat_seq, insert.length(), insert);
                    if (debug) { cout << "Comparing:  " << insert << "\t" << repeat_seq << "\tMutations: " << d << "\n"; }
                    if (d < min_d) { min_d = d; }
                }
            }
        }
    }

    if (debug) {
        cout << "\nLeast Mutations: " << min_d << "\nInsert Length: " << s << "\n\n";
        if (terminate) { cout << "*** Termination ***\n" ; }
    }
    return vector<uint> { min_d, s, terminate };
}

vector<uint> extension_mutations(utils::RepeatTracker rtracker) {
    if (debug) { cout << "\n*** Extension ***\n"; }
    uint terminate = 0;

    string valid_motif = rtracker.valid_motif;
    string insert = rtracker.insert;
    uint s = insert.length();
    uint m = motif_size;

    uint min_d = -1;
    uint plen = 0;
    for (uint u = 1; u <= 2*m; u++) {
        string repeat_seq = utils::expand_repeat(valid_motif, u);
        uint d = levenshtein_distance(repeat_seq.length(), repeat_seq, insert.length(), insert);
        if (d < min_d) { min_d = d; plen = u; }
    }

    uint start = rtracker.start;
    uint end   = rtracker.end;
    uint rlen  = end - start;
    uint muts  = rtracker.mutations;

    uint awdmuts_wins = fraction_mutations * (rlen + s);
    uint remain_muts  = awdmuts_wins - muts;

    if (remain_muts == 0) {
        if (debug) { cout << "Remaining mutations: 0\n" << "*** Termination ***\n\n" ; }
        s = 0; min_d = 0; terminate = 1;
        return vector<uint> { min_d, s, terminate };
    }

    if (min_d > remain_muts) {
        terminate = 1;
        min_d = -1;
        if (remain_muts == 0) { s = 0; min_d = 0; }
        else {
            while (min_d > remain_muts) {
                insert = insert.substr(0, insert.length()-1);
                s = insert.length();
                if (s == 0) { min_d = 0; break; }
                for (uint u = s-remain_muts; u <= s+remain_muts; u++) {
                    string repeat_seq = utils::expand_repeat(valid_motif, u);
                    uint d = levenshtein_distance(repeat_seq.length(), repeat_seq, insert.length(), insert);
                    if (d < min_d) { min_d = d; }
                }
            }
        }
    }

    if (debug) { cout << "Least Mutations: " << min_d << "\nInsert Length: " << s << "\n\n"; }
    return vector<uint> { min_d, s , terminate, plen };
}

int main(int argc, char* argv[]) {
    string fin = argv[1];
    string line;
    ifstream ins(fin);
    string seq_name;
    utils::SequenceWindow window;
    while(getline(ins, line)) {
        if (line[0] == '>') {
            seq_name = line.substr(1, line.find(' ')-1);
            window.reset();
        }
        else {
            for (const auto c: line) {
                char curr_nuc = toupper(c);
                window.update(curr_nuc);

                if (window.count >= 2*motif_size) {
                    uint position = window.count - 2*motif_size;
                    curr_nuc = window.sequence[motif_size-1];
                    if (debug) {
                        cout << "\n\n******************************  Position: " << position << "  ******************************\n";
                        cout << "Nucleotide: " << curr_nuc << "\n" << "Motif: " << window.motif << "\n";
                    }
                    string curr_rclass = utils::get_repeat_class(window.motif, rClassMap);
                    uint curr_rclass_first_check = 0;
                    vector<string> drop_rclasses;

                    if (globalRepeatTracker.find(curr_rclass) == globalRepeatTracker.end()) {
                        curr_rclass_first_check = 1;
                        globalRepeatTracker[curr_rclass] = utils::RepeatTracker();
                        globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count - motif_size, window.next);
                    }

                    std::unordered_map<string, utils::RepeatTracker>::iterator iter = globalRepeatTracker.begin();
                    std::unordered_map<string, utils::RepeatTracker>::iterator end_iter = globalRepeatTracker.end();
                    for(; iter != end_iter; ++iter) {
                        uint drop_rclass = 0;
                        string rclass = iter->first;
                        globalRepeatTracker[rclass].next_motif = window.next;
                        if (debug) { cout << "\n==========  " << rclass << "  ==========\n"; }
                        if (rclass == curr_rclass) {
                            if (curr_rclass_first_check == 0) {
                                uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                                string valid_motif = globalRepeatTracker[rclass].valid_motif;
                                char valid_nuc = globalRepeatTracker[rclass].valid_nuc;

                                if (window.motif == valid_motif) {
                                    if (rclass_continue) {
                                        if (debug) { cout << "NN: True\tRC: True" << "\n"; }
                                        // repeat continuation
                                        globalRepeatTracker[rclass].valid_motif = valid_motif.substr(1) + valid_motif[0];
                                        globalRepeatTracker[rclass].valid_nuc = globalRepeatTracker[rclass].valid_motif[0];
                                        globalRepeatTracker[rclass].repeat += curr_nuc;
                                        globalRepeatTracker[rclass].end = window.count - motif_size;
                                    }
                                    else {
                                        if (debug) { cout << "NN: True\tRC: False" << "\n"; }
                                        // repeat continuation after insertion
                                        vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass]);
                                        uint d = insertion_result[0];
                                        uint insert_len = insertion_result[1];
                                        uint terminate = insertion_result[2];

                                        if (terminate) {
                                            uint start = globalRepeatTracker[rclass].start;
                                            globalRepeatTracker[rclass].end += insert_len;
                                            uint end = globalRepeatTracker[rclass].end;
                                            uint rlen = end - start;
                                            globalRepeatTracker[rclass].mutations += d;
                                            uint muts = globalRepeatTracker[rclass].mutations;
                                            globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, insert_len);
                                            if (rlen >= 12) {
                                                if (debug) {
                                                    cout << "*** Valid repeat ***\n";
                                                    globalRepeatTracker[rclass].print();
                                                }
                                                else {
                                                    cout << seq_name << "\t" << start << "\t" << end << "\t" << rclass << "\t" \
                                                    << rlen << "\t" << muts << "\n";
                                                }
                                            }

                                            if (rclass == curr_rclass) {
                                                globalRepeatTracker[rclass].initialise(window.motif, position, window.count - motif_size, window.next);
                                            }
                                            else { drop_rclass = 1; }
                                        }

                                        else {
                                            globalRepeatTracker[rclass].valid_motif = valid_motif.substr(1) + valid_motif[0];
                                            globalRepeatTracker[rclass].valid_nuc = globalRepeatTracker[rclass].valid_motif[0];
                                            globalRepeatTracker[rclass].end = window.count - motif_size;
                                            globalRepeatTracker[rclass].repeat += (globalRepeatTracker[rclass].insert + curr_nuc);
                                            globalRepeatTracker[rclass].mutations += d;
                                            globalRepeatTracker[rclass].interrupt = 0;
                                        }
                                    }
                                }
                                else {
                                    if (position < globalRepeatTracker[rclass].end) {
                                        if (debug) {
                                            cout << "*** Cycle found mid motif ***\n";
                                            cout << "*** " << valid_motif << " encountered as " << window.motif << " ***\n";
                                            globalRepeatTracker[rclass].print();
                                            cout << "\n";
                                        }
                                        globalRepeatTracker[rclass].interrupt = 0;
                                        globalRepeatTracker[rclass].end = window.count - motif_size;
                                        globalRepeatTracker[rclass].valid_motif = window.motif.substr(1) + window.motif[0];
                                    }
                                    else if (position == globalRepeatTracker[rclass].end) {
                                        if (debug) {
                                            cout << "*** Cycle found book ended ***\n";
                                            cout << "*** " << valid_motif << " encountered as " << window.motif << " ***\n";
                                            globalRepeatTracker[rclass].print();
                                            cout << "\n";
                                        }
                                        globalRepeatTracker[rclass].interrupt = 0;
                                        globalRepeatTracker[rclass].end = window.count - motif_size;
                                        globalRepeatTracker[rclass].valid_motif = window.motif.substr(1) + window.motif[0];
                                    }
                                    else {
                                        if (debug) {
                                            cout << "*** Cycle found after insertion ***\n";
                                            cout << "*** " << valid_motif << " encountered as " << window.motif << " ***\n";
                                            globalRepeatTracker[rclass].print();
                                            cout << "\n";
                                        }
                                        string cycle = "";
                                        uint c = 1;
                                        for (; c < motif_size; c++) {
                                            cycle = valid_motif.substr(c) + valid_motif.substr(0, c);
                                            if (cycle == window.motif) { break; }
                                        }
                                        string insert = globalRepeatTracker[rclass].insert;
                                        uint s = insert.length();
                                        cout << "Insert length: " << c << "\n"; 
                                        if (debug) {
                                            cout << "*** Match " << insert.substr(0, s-motif_size) << " with " << valid_motif+cycle.substr(motif_size-c) << " ***\n";
                                        }
                                    }
                                }
                            }
                        }
                        // if ((rclass == curr_rclass) && curr_rclass_first_check) {
                        //     // Current repeat class first occurrence
                        // }
                        else {
                            uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                            string valid_motif = globalRepeatTracker[rclass].valid_motif;
                            char valid_nuc = globalRepeatTracker[rclass].valid_nuc;
                            globalRepeatTracker[rclass].curr_motif = window.motif;

                            if (rclass_continue) {
                                if (debug) { cout << "NN: False\tRC: True" << "\n"; }
                                globalRepeatTracker[rclass].interrupt = 1;
                                globalRepeatTracker[rclass].insert = curr_nuc;
                                globalRepeatTracker[rclass].valid_motif = valid_motif[motif_size-1] + valid_motif.substr(0, motif_size-1);
                            }

                            else {
                                if (debug) { cout << "NN: False\tRC: False" << "\n"; }
                                globalRepeatTracker[rclass].insert += curr_nuc;
                                // if (globalRepeatTracker[rclass].insert.length() == motif_size) {
                                //     if (debug) { cout << "\n+++ Unit extension +++\n"; }
                                //     vector<uint> extension_result = extension_mutations(globalRepeatTracker[rclass]);
                                //     uint d = extension_result[0];
                                //     uint insert_len = extension_result[1];
                                //     uint terminate = extension_result[2];
                                //     uint plen = extension_result[3];

                                //     if (terminate) {
                                //         uint start = globalRepeatTracker[rclass].start;
                                //         globalRepeatTracker[rclass].end += insert_len;
                                //         uint end = globalRepeatTracker[rclass].end;
                                //         uint rlen = end - start;
                                //         globalRepeatTracker[rclass].mutations += d;
                                //         uint muts = globalRepeatTracker[rclass].mutations;
                                //         globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, insert_len);
                                //         if (rlen >= 12) {
                                //             if (debug) {
                                //                 cout << "*** Valid repeat ***\n";
                                //                 globalRepeatTracker[rclass].print();
                                //             }
                                //             else {
                                //                 cout << seq_name << "\t" << start << "\t" << end << "\t" << rclass << "\t" \
                                //                 << rlen << "\t" << muts << "\n";
                                //             }
                                //         }

                                //         if (rclass == curr_rclass) {
                                //             globalRepeatTracker[rclass].initialise(window.motif, position, window.count - motif_size, window.next);
                                //         }
                                //         else { drop_rclass = 1; }
                                //     }
                                //     else {
                                //         uint check = -1;

                                //         if (d != check) {
                                //             if (plen < motif_size) {
                                //                 globalRepeatTracker[rclass].valid_motif = valid_motif.substr(plen, valid_motif.length()) + valid_motif.substr(0, plen);
                                //             }
                                //             else {
                                //                 string imotif = globalRepeatTracker[rclass].valid_motif;
                                //                 globalRepeatTracker[rclass].valid_motif = utils::expand_repeat(imotif, plen).substr(plen-motif_size, plen);
                                //             }
                                //             globalRepeatTracker[rclass].valid_nuc = globalRepeatTracker[rclass].valid_motif[0];
                                //             globalRepeatTracker[rclass].end = window.count - motif_size;
                                //             globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert;
                                //             globalRepeatTracker[rclass].mutations += d;
                                //             globalRepeatTracker[rclass].interrupt = 0;
                                //         }
                                //     }
                                // }
                            }
                        }
                        if (drop_rclass) { drop_rclasses.push_back(rclass); }
                        else if (debug) {
                            globalRepeatTracker[rclass].print();
                        }
                    }
                    for (uint j=0; j<drop_rclasses.size(); j++) {
                        globalRepeatTracker.erase(drop_rclasses[j]);
                    }
                }
            }
        }
    }
    ins.close();
}