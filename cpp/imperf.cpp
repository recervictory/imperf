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
const uint motif_size = 4;
const float fraction_mutations = 0.1;
uint debug = 0;

vector<uint> insertion_mutations(utils::RepeatTracker rtracker, string valid_insert, uint terminate) {

    string valid_motif = rtracker.valid_motif;    
    uint start = rtracker.start;
    uint end   = rtracker.end;
    uint rlen  = end - start;
    uint runits = rlen / motif_size;
    uint muts  = rtracker.mutations;

    string insert = rtracker.insert;
    uint ilen = insert.length();
    uint iunits = ilen / motif_size;
    uint m = motif_size;
    uint min_d = -1;
    uint s = 0;

    // if repeat already reached threshold mutations
    if (muts == rlen * fraction_mutations) {
        min_d = 0, s = 0;
        return { min_d, s, 1 };
    }

    uint threshold_muts = fraction_mutations * (rlen + ilen);
    uint remain_muts  = threshold_muts - muts;

    if (debug) {
        cout << "\n*** Insertion " << insert << " ***\n";
        cout << "*** Valid " << valid_insert << " ***\n";
        cout << "*** Remaining mutations: " << remain_muts << " ***\n";
    }

    string ext_motif = valid_insert.substr(valid_insert.length()-motif_size, motif_size);
    string prefix = valid_motif.substr(0, valid_insert.length()-motif_size);

    if (!terminate) {
        uint valid_lrlmt = ilen - remain_muts;
        if (remain_muts > ilen) { valid_lrlmt = 0; }
        uint valid_uplmt = ilen + remain_muts;
        uint ext_units = 0;
        if (valid_uplmt > prefix.length()) {
            ext_units = (valid_uplmt - prefix.length()) / motif_size;
        }
        string insert_comp = prefix + utils::expand_repeat(ext_motif, ext_units*motif_size);
        int x = insert_comp.length();
        int y = insert.length();
        int z = (x+1)*(y+1);
        int d_matrix[z];
        int i,j,substitution_cost;
        d_matrix[0+0*(y+1)] = 0;
        for ( i = 1; i <= x; i++ ) { d_matrix[i+0*(x+1)] = i; }
        for ( j = 1; j <= y; j++ ) { d_matrix[0+j*(x+1)] = j; }
        for ( j = 1; j <= y; j++ ) {
            for ( i = 1; i <= x; i++ ) {
                if ( insert_comp[i-1] == insert[j-1] ) { substitution_cost = 0; }
                else { substitution_cost = 1; }
                d_matrix[i+j*(x+1)] = i4_min ( d_matrix[i-1+j*(x+1)] + 1,
                                i4_min ( d_matrix[i+(j-1)*(x+1)] + 1, 
                                        d_matrix[i-1+(j-1)*(x+1)] + substitution_cost ) );
            }
        }
        uint d = d_matrix[((x+1)*(y+1))-1];
        if (debug) {
            cout << "Comparing:  " << insert << "\t" << insert_comp << "\tMutations: " << d << "\n";
            cout << "*** Levenshtein Matrix ***\n";
            for (int j=0; j<z; j++) {
                if (insert_comp.length()>0 && j > 0) {
                    if (j % (x+1) == 0) { cout << "\n"; }
                }
                cout << d_matrix[j] << "\t";
            }
            cout << "\n";
        }
        int idx = z-1;
        uint row_min = -1;
        int min_idx = z-1;
        for (; idx>=0; idx--) {
            int mut = d_matrix[idx];
            if (mut < row_min) { row_min = mut; min_idx = idx; }
            if (idx % (x+1) == 0) {
                uint dynamic_insert_len = idx / (x+1);
                if (dynamic_insert_len < row_min) { row_min = dynamic_insert_len; }
                uint dynamic_remain_muts = ((rlen + dynamic_insert_len )*fraction_mutations) - muts;
                if (row_min <= dynamic_remain_muts) {
                    break;
                } 
            }
        }
        uint valid_insert_len = min_idx / (x+1);
        uint valid_repeat_len = min_idx % (x+1);
        if (debug) {
            cout << "*** Index " << min_idx << " ***\n";
            cout << "*** Valid insert length: " << valid_insert_len << " ***\n";
            cout << "*** Valid repeat length: " << valid_repeat_len << " ***\n";
        }
        s = valid_insert_len;
        min_d = row_min;

        if (s < ilen) { terminate = 1; }
        else { terminate = 0; }
        return { min_d, s, terminate };
    }
    else {
        uint valid_lrlmt = ilen - remain_muts;
        if (remain_muts > ilen) { valid_lrlmt = 0; }
        uint valid_uplmt = ilen + remain_muts;
        string insert_comp = utils::expand_repeat(ext_motif, valid_uplmt);
        int x = insert_comp.length();
        int y = insert.length();
        int z = (x+1)*(y+1);
        int d_matrix[z];
        int i,j,substitution_cost;
        d_matrix[0+0*(y+1)] = 0;
        for ( i = 1; i <= x; i++ ) { d_matrix[i+0*(x+1)] = i; }
        for ( j = 1; j <= y; j++ ) { d_matrix[0+j*(x+1)] = j; }
        for ( j = 1; j <= y; j++ ) {
            for ( i = 1; i <= x; i++ ) {
                if ( insert_comp[i-1] == insert[j-1] ) { substitution_cost = 0; }
                else { substitution_cost = 1; }
                d_matrix[i+j*(x+1)] = i4_min ( d_matrix[i-1+j*(x+1)] + 1,
                                i4_min ( d_matrix[i+(j-1)*(x+1)] + 1, 
                                        d_matrix[i-1+(j-1)*(x+1)] + substitution_cost ) );
            }
        }
        uint d = d_matrix[((x+1)*(y+1))-1];
        if (debug) {
            cout << "Comparing:  " << insert << "\t" << insert_comp << "\tMutations: " << d << "\n";
            cout << "*** Levenshtein Matrix ***\n";
            for (int j=0; j<z; j++) {
                if (insert_comp.length()>0 && j > 0) {
                    if (j % (x+1) == 0) { cout << "\n"; }
                }
                cout << d_matrix[j] << "\t";
            }
            cout << "\n";
        }
        int idx = z-1;
        uint row_min = -1;
        int min_idx = z-1;
        for (; idx>=0; idx--) {
            int mut = d_matrix[idx];
            if (mut < row_min) { row_min = mut; min_idx = idx; }
            if (idx % (x+1) == 0) {
                uint dynamic_insert_len = idx / (x+1);
                if (dynamic_insert_len < row_min) { row_min = dynamic_insert_len; }
                uint dynamic_remain_muts = ((rlen + dynamic_insert_len )*fraction_mutations) - muts;
                if (row_min <= dynamic_remain_muts) {
                    break;
                } 
            }
        }
        uint valid_insert_len = min_idx / (x+1);
        uint valid_repeat_len = min_idx % (x+1);
        if (debug) {
            cout << "*** Index " << min_idx << " ***\n";
            cout << "*** Valid insert length: " << valid_insert_len << " ***\n";
            cout << "*** Valid repeat length: " << valid_repeat_len << " ***\n";
        }
        s = valid_insert_len;
        min_d = row_min;

        if (s == ilen && min_d < remain_muts) {
            if (debug) { cout << "*** Continue repeat ***\n"; }
            terminate = 0;
            uint c = valid_repeat_len % motif_size;
            return { min_d, s, terminate, c };
        }
    }

    if (debug) {
        cout << "\nLeast Mutations: " << min_d << "\nInsert Length: " << s << "\n\n";
        if (terminate) { cout << "*** Termination ***\n" ; }
    }
    return vector<uint> { min_d, s, terminate };
}

int main(int argc, char* argv[]) {
    string fin = argv[1];
    string line;
    ifstream ins(fin);
    string seq_name;
    utils::SequenceWindow window;

    uint position_offset = 2*motif_size;
    uint end_offset = motif_size;

    while(getline(ins, line)) {
        if (line[0] == '>') {
            seq_name = line.substr(1, line.find(' ')-1);
            window.reset();
        }
        else {
            for (const auto c: line) {
                char curr_nuc = toupper(c);
                window.update(curr_nuc);

                if (window.count >= position_offset) {
                    uint position = window.count - position_offset;
                    curr_nuc = window.sequence[motif_size-1];
                    
                    if (debug) {
                        cout << "\n\n******************************  Position: " << position << "  ******************************\n";
                        cout << "Nucleotide: " << curr_nuc << "\n" << "Motif: " << window.motif << "\n";
                    }

                    string curr_motif = window.motif;
                    string curr_rclass = utils::get_repeat_class(curr_motif, rClassMap);
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
                                if (position < globalRepeatTracker[rclass].end) {
                                    if (debug && (curr_motif != valid_motif)) {
                                        cout << "*** Cycle found mid motif ***\n";
                                        cout << "*** " << valid_motif << " encountered as " << window.motif << " ***\n\n";
                                    }
                                    globalRepeatTracker[rclass].interrupt = 0;
                                    globalRepeatTracker[rclass].insert = "";
                                    globalRepeatTracker[rclass].end = window.count - end_offset;
                                    globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                    globalRepeatTracker[rclass].valid_nuc = globalRepeatTracker[rclass].valid_motif[0];
                                }
                                else if (position == globalRepeatTracker[rclass].end) {
                                    if (debug && (curr_motif != valid_motif)) {
                                        cout << "*** Cycle found book ended ***\n";
                                        cout << "*** " << valid_motif << " encountered as " << curr_motif << " ***\n\n";
                                    }
                                    globalRepeatTracker[rclass].interrupt = 0;
                                    globalRepeatTracker[rclass].insert = "";
                                    globalRepeatTracker[rclass].end = window.count - end_offset;
                                    globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                    globalRepeatTracker[rclass].valid_nuc = globalRepeatTracker[rclass].valid_motif[0];
                                }
                                else {
                                    string insert = globalRepeatTracker[rclass].insert;
                                    uint insert_len = insert.length();
                                    string valid_insert = "";
                                    if (debug && (curr_motif != valid_motif)) {
                                        cout << "*** Cycle found after insertion: " << insert << "***\n";
                                        cout << "*** " << valid_motif << " encountered as " << curr_motif << " ***\n\n";
                                    }
                                    string cycle = "";
                                    uint c = 0;
                                    for (; c < motif_size; c++) {
                                        cycle = valid_motif.substr(c) + valid_motif.substr(0, c);
                                        if (cycle == window.motif) { break; }
                                    }
                                    valid_insert = valid_motif.substr(0,c) + cycle;
                                    
                                    if (debug) {
                                        cout << "*** Match " << insert.substr(0, insert_len-motif_size) << " with " << valid_insert << " ***\n";
                                    }
                                    uint terminate = 0;
                                    vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], valid_insert, terminate);
                                    uint d = insertion_result[0];
                                    uint s = insertion_result[1];
                                    terminate = insertion_result[2];
                                    if (terminate) {
                                        globalRepeatTracker[rclass].end += s;
                                        globalRepeatTracker[rclass].mutations += d;
                                        uint start = globalRepeatTracker[rclass].start;
                                        uint end = globalRepeatTracker[rclass].end;
                                        uint rlen = end - start;
                                        uint muts = globalRepeatTracker[rclass].mutations;
                                        if (rlen > 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass << "\t" << muts << "\t" << globalRepeatTracker.size() << "\n"; 
                                        }
                                        globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count - motif_size, window.next);
                                    }
                                    else {
                                        globalRepeatTracker[rclass].mutations += d;
                                        globalRepeatTracker[rclass].interrupt = 0;
                                        globalRepeatTracker[rclass].insert = "";
                                        globalRepeatTracker[rclass].end = window.count - end_offset;
                                        globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                        globalRepeatTracker[rclass].valid_nuc = globalRepeatTracker[rclass].valid_motif[0];
                                    }
                                }
                            }
                        }
                        
                        else {
                            uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                            string valid_motif = globalRepeatTracker[rclass].valid_motif;
                            char valid_nuc = globalRepeatTracker[rclass].valid_nuc;
                            globalRepeatTracker[rclass].curr_motif = window.motif;

                            if (rclass_continue) {
                                globalRepeatTracker[rclass].interrupt = 1;
                                globalRepeatTracker[rclass].insert = curr_nuc;
                                globalRepeatTracker[rclass].valid_motif = valid_motif[motif_size-1] + valid_motif.substr(0, motif_size-1);
                            }

                            else {
                                globalRepeatTracker[rclass].insert += curr_nuc;

                                string insert = globalRepeatTracker[rclass].insert;
                                uint start = globalRepeatTracker[rclass].start;
                                uint end = globalRepeatTracker[rclass].end;
                                uint ilen = insert.length();
                                uint iunits = ilen / motif_size;
                                uint rlen = end - start;
                                uint runits = rlen / motif_size;
                                uint muts = globalRepeatTracker[rclass].mutations;

                                if ((ilen >= 4*motif_size) || (ilen >= rlen)) {
                                    if (debug) { cout << "*** Terminate this repeat ***\n"; }
                                    uint terminate = 1;
                                    vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], valid_motif, terminate);
                                    uint d = insertion_result[0];
                                    uint s = insertion_result[1];
                                    terminate = insertion_result[2];
                                    globalRepeatTracker[rclass].end += s;
                                    globalRepeatTracker[rclass].mutations += d;
                                    start = globalRepeatTracker[rclass].start;
                                    end = globalRepeatTracker[rclass].end;
                                    rlen = end - start;
                                    muts = globalRepeatTracker[rclass].mutations;
                                    if (terminate) {
                                        if (rlen > 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass << "\t" << muts << "\t" << globalRepeatTracker.size() << "\n"; 
                                        } 
                                        drop_rclass = 1;
                                    }
                                    else {
                                        uint c = insertion_result[3];
                                        globalRepeatTracker[rclass].valid_motif = valid_motif.substr(0, c) + valid_motif.substr(c);
                                        globalRepeatTracker[rclass].insert = "";
                                    }
                                }
                            }
                        }
                        
                        if (drop_rclass) { drop_rclasses.push_back(rclass); }
                        else if (debug) { globalRepeatTracker[rclass].print(); }
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