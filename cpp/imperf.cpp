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
unordered_map<string, utils::repeat_tracker> globalRepeatTracker;
const uint motif_size = 6;
const float fraction_mutations = 0.1;
uint debug = 0;

vector<uint> insertion_mutations(utils::repeat_tracker rtracker) {
    uint terminate = 0;

    uint start = rtracker.start;
    uint end   = rtracker.end;
    uint rlen  = end - start + motif_size;
    uint muts  = rtracker.mutations;

    string imotif = rtracker.inter_motif;
    string insert = rtracker.insertion;
    uint s = insert.length();
    uint m = motif_size;

    uint awdmuts_wins = fraction_mutations * (rlen + s);
    uint remain_muts  = awdmuts_wins - muts;
    uint refrep_ul = ((s + remain_muts)/m);
    uint refrep_ll = 0;
    if (remain_muts < s) { refrep_ll = ((s - remain_muts)/m); }
    uint min_d = -1;
    for (uint u = refrep_ll; u <= refrep_ul; u++) {
        string repeat_seq = utils::expand_repeat(imotif, u*m);
        uint d = levenshtein_distance(repeat_seq.length(), repeat_seq, insert.length(), insert);
        if (d < min_d) { min_d = d; }
    }

    if (min_d > s) { min_d = s; }
    if (min_d > remain_muts) {
        terminate = 1;
        min_d = -1;
        if (remain_muts == 0) { s = 0; min_d = 0; }
        else {
            while (min_d > remain_muts) {
                insert = insert.substr(0,insert.length()-1);
                s = insert.length();
                if (s == 0) { min_d = 0; break; }
                for (uint u = s-remain_muts; u <= s+remain_muts; u++) {
                    string repeat_seq = utils::expand_repeat(imotif, u);
                    uint d = levenshtein_distance(repeat_seq.length(), repeat_seq, insert.length(), insert);
                    if (d < min_d) { min_d = d; }
                }
            }
        }
    }

    return vector<uint> { min_d, s , terminate };
}

vector<uint> extension_mutations(utils::repeat_tracker rtracker) {
    uint terminate = 0;

    string imotif = rtracker.inter_motif;
    string insert = rtracker.insertion;
    uint s = insert.length();
    uint m = motif_size;

    uint min_d = -1;
    uint plen = 0;
    for (uint u = 1; u <= 2*m; u++) {
        string repeat_seq = utils::expand_repeat(imotif, u);
        uint d = levenshtein_distance(repeat_seq.length(), repeat_seq, insert.length(), insert);
        if (d < min_d) { min_d = d; plen = u; }
    }

    uint start = rtracker.start;
    uint end   = rtracker.end;
    uint rlen  = end - start + motif_size;
    uint muts  = rtracker.mutations;

    uint awdmuts_wins = fraction_mutations * (rlen + s);
    uint remain_muts  = awdmuts_wins - muts;

    if (min_d > remain_muts) {
        terminate = 1;
        min_d = -1;
        if (remain_muts == 0) { s = 0; min_d = 0; }
        else {
            while (min_d > remain_muts) {
                insert = insert.substr(0,insert.length()-1);
                s = insert.length();
                if (s == 0) { min_d = 0; break; }
                for (uint u = s-remain_muts; u <= s+remain_muts; u++) {
                    string repeat_seq = utils::expand_repeat(imotif, u);
                    uint d = levenshtein_distance(repeat_seq.length(), repeat_seq, insert.length(), insert);
                    if (d < min_d) { min_d = d; }
                }
            }
        }
    }

    return vector<uint> { min_d, s , terminate, plen };
}

int main(int argc, char* argv[]) {
    string fin = argv[1];
    string line;
    ifstream ins(fin);
    string seq_name;
    utils::seqWindow window;
    while(getline(ins, line)) {
        if (line[0] == '>') {
            seq_name = line.substr(1, line.find(' ')-1);
            window.reset();
        }
        else {
            for (const auto c: line) {
                char curr_nuc = toupper(c);
                window.update(curr_nuc);

                if (window.count >= motif_size) {
                    if (debug) {
                        cout << "\n\n******************************  ";
                        cout << "Position: " << window.count - motif_size;
                        cout << "  ******************************\n";
                        cout << "Nucleotide: " << curr_nuc << "\n" << "Motif: " << window.motif << "\n";
                    }
                    string curr_rclass = utils::get_repeat_class(window.motif, rClassMap);
                    uint curr_rclass_first_check = 0;
                    vector<string> drop_rclasses;

                    if (globalRepeatTracker.find(curr_rclass) == globalRepeatTracker.end()) {
                        curr_rclass_first_check = 1;
                        globalRepeatTracker[curr_rclass] = utils::repeat_tracker();
                        globalRepeatTracker[curr_rclass].reset(window.motif, window.count);
                    }

                    std::unordered_map<string, utils::repeat_tracker>::iterator iter = globalRepeatTracker.begin();
                    std::unordered_map<string, utils::repeat_tracker>::iterator end_iter = globalRepeatTracker.end();
                    for(; iter != end_iter; ++iter) {
                        uint drop_rclass = 0;
                        string rclass = iter->first;
                        if (debug) {
                            cout << "\n==========  " << rclass << "  ==========\n";
                        }
                        if ((rclass == curr_rclass) && curr_rclass_first_check) {
                            // pass
                        }
                        else {
                            uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                            string prev_motif = globalRepeatTracker[rclass].prev_motif;
                            char next_nuc = globalRepeatTracker[rclass].next_nuc;

                            if (next_nuc == curr_nuc) {
                                if (rclass_continue) {
                                    if (debug) { cout << "NN: True\tRC: True" << "\n"; }
                                    // repeat continuation
                                    globalRepeatTracker[rclass].prev_motif = prev_motif.substr(1) + prev_motif[0];
                                    globalRepeatTracker[rclass].next_nuc = globalRepeatTracker[rclass].prev_motif[0];
                                    globalRepeatTracker[rclass].end = window.count;
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
                                        uint end   = globalRepeatTracker[rclass].end + motif_size + insert_len;
                                        uint rlen = end - start;
                                        uint muts = globalRepeatTracker[rclass].mutations + d;
                                        if (rlen >= 12) {
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rclass << "\t" \
                                            << rlen << "\t" << muts << "\n";
                                        }

                                        if (rclass == curr_rclass) {
                                            globalRepeatTracker[rclass].reset(window.motif, window.count);
                                        }
                                        else { drop_rclass = 1; }
                                    }

                                    else {
                                        globalRepeatTracker[rclass].prev_motif = prev_motif.substr(1) + prev_motif[0];
                                        globalRepeatTracker[rclass].next_nuc = globalRepeatTracker[rclass].prev_motif[0];
                                        globalRepeatTracker[rclass].end = window.count;
                                        globalRepeatTracker[rclass].mutations += d;
                                        globalRepeatTracker[rclass].interrupt = 0;
                                    }
                                }
                            }

                            else {

                                if (rclass_continue) {
                                    if (debug) { cout << "NN: False\tRC: True" << "\n"; }
                                    globalRepeatTracker[rclass].inter_motif = globalRepeatTracker[rclass].prev_motif;
                                    globalRepeatTracker[rclass].insertion = curr_nuc;
                                    globalRepeatTracker[rclass].interrupt = 1;
                                }

                                else {
                                    if (debug) { cout << "NN: False\tRC: False" << "\n"; }
                                    globalRepeatTracker[rclass].insertion += curr_nuc;
                                    if (globalRepeatTracker[rclass].insertion.length() == motif_size) {
                                        if (debug) { cout << "\n+++ Unit extension +++\n"; }
                                        vector<uint> extension_result = extension_mutations(globalRepeatTracker[rclass]);
                                        uint d = extension_result[0];
                                        uint insert_len = extension_result[1];
                                        uint terminate = extension_result[2];
                                        uint plen = extension_result[3];

                                        if (terminate) {
                                            uint start = globalRepeatTracker[rclass].start;
                                            uint end   = globalRepeatTracker[rclass].end + motif_size + insert_len;
                                            uint rlen = end - start;
                                            uint muts = globalRepeatTracker[rclass].mutations + d;
                                            if (rlen >= 12) {
                                                cout << seq_name << "\t" << start << "\t" << end << "\t" << rclass << "\t" \
                                                << rlen << "\t" << muts << "\n";
                                            }

                                            if (rclass == curr_rclass) {
                                                globalRepeatTracker[rclass].reset(window.motif, window.count);
                                            }
                                            else { drop_rclass = 1; }
                                        }
                                        else {
                                            uint check = -1;

                                            if (d != check) {
                                                if (plen < motif_size) {
                                                    globalRepeatTracker[rclass].prev_motif = prev_motif.substr(plen, prev_motif.length()) + prev_motif.substr(0, plen);
                                                }
                                                else {
                                                    string imotif = globalRepeatTracker[rclass].inter_motif;
                                                    globalRepeatTracker[rclass].prev_motif = utils::expand_repeat(imotif, plen).substr(plen-motif_size, plen);
                                                }
                                                globalRepeatTracker[rclass].next_nuc = globalRepeatTracker[rclass].prev_motif[0];
                                                globalRepeatTracker[rclass].end = window.count;
                                                globalRepeatTracker[rclass].mutations += d;
                                                globalRepeatTracker[rclass].interrupt = 0;
                                            }
                                        }
                                    }
                                }

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