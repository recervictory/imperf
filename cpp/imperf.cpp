#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "utils.h"

using namespace std;
using namespace utils;

// global variables
unordered_map<string, string> rClassMap;
unordered_map<string, utils::RepeatTracker> globalRepeatTracker;
const uint motif_size = 6;
const float fraction_mutations = 0.1;
uint debug = 0;


int i4_min ( int i1, int i2 ) {
  if ( i1 < i2 ) { return i1; }
  else { return i2; }
}


void optimal_insert(int x, int y, uint rlen, uint muts, int distance_matrix[], int &min_idx, uint &row_least_muts) {
    uint z = (x+1)*(y+1);       // matrix dimensions
    int idx = z-1;              // initialising index
    row_least_muts = -1;        // least mutations of the row i.e., a repeat length
    min_idx = z-1;              // index of matrix at least mutations
    for (; idx>=0; idx--) {
        
        int idx_muts = distance_matrix[idx];
        // storing the least mutations for a particular insert length
        if (idx_muts < row_least_muts) {
            row_least_muts = idx_muts;
            min_idx = idx;
        }

        // after completely scanning the mutations in a row i.e.,
        // identifying the least number of mutations for a particular insert length
        if (idx % (x+1) == 0) {
            uint dynamic_ilen = idx / (x+1);            // dynamic insert length
            
            // if the least mutations for the dynamic insert length is 
            // greater than dynamic insert length  
            if (row_least_muts > dynamic_ilen) {
                row_least_muts = dynamic_ilen;
            }

            // allowed mutations for dynamic insert length
            uint dynamic_remain_muts = ((rlen + dynamic_ilen )*fraction_mutations) - muts;
            if (debug) {
                cout << "\nDynamic insert length:     " << dynamic_ilen << "\n";
                cout << "Remaining mutations:         " << dynamic_remain_muts << "\n";
                cout << "Minimum mutations:           " << row_least_muts << "\n";
            }

            // The least mutations for the insert is in the allowed range of mutations
            if (row_least_muts <= dynamic_remain_muts) {
                if (debug) {
                    cout << "Least mutations:      " << row_least_muts << " at " << min_idx << "\n";
                    cout << "Greedy insert length: " << min_idx / (x+1) << "\n";
                    cout << "Greedy repeat length: " << min_idx % (x+1) << "\n\n";
                }

                // Trimming the insert to exclude mutations in the end of the insert
                int a = min_idx, b = row_least_muts;
                while (a >= (x+2)) {
                    b = distance_matrix[a];
                    if (debug) {
                        cout << "Comparing indexes " << a << " and " << (a - (x+2)) << "\n";
                        cout << "Comparing mutations " << b << " and " << distance_matrix[a - (x+2)] << "\n\n";
                    }

                    // If the edit distance for insert length is equal to 
                    // edit distance of the insert trimmed by one base
                    if (distance_matrix[a - (x+2)] == b) {
                        if (debug) {
                            cout << "Non-greedy insert length: " << a / (x+1) << "\n";
                            cout << "Non-greedy repeat length: " << a % (x+1) << "\n";
                        }
                        break;
                    }
                    else { a = a - (x + 2); }
                }
                
                if (a < (x+2)) { min_idx = 0; row_least_muts = 0; }
                else { min_idx = a; row_least_muts = b; }
                break;
            } 
        }
    }
}


vector<uint> insertion_mutations(utils::RepeatTracker rtracker, string valid_insert, uint terminate) {

    string valid_motif = rtracker.valid_motif;    
    uint start = rtracker.start;
    uint end   = rtracker.end;
    uint rlen  = end - start;
    uint muts  = rtracker.mutations;

    string insert = rtracker.insert;
    uint ilen = insert.length();
    uint m = motif_size;

    // Final minimum mutations and insert length
    uint min_d = -1, s = 0;

    // if repeat already reached threshold mutations
    if (muts == rlen * fraction_mutations) {
        min_d = 0, s = 0; terminate = 1;
        return { min_d, s, terminate };
    }

    uint threshold_muts = fraction_mutations * (rlen + ilen);
    int remain_muts  = threshold_muts - muts;

    if (debug) {
        cout << "*** Insertion: " << insert << " ***";
        cout << "\n*** Valid: " << valid_insert << " ***";
        cout << "\n*** Remaining mutations: " << remain_muts << " ***\n";
    }

    if (remain_muts <= 0) {
        min_d = 0, s = 0; terminate = 1;
        return { min_d, s, terminate };
    }

    string prefix = valid_motif.substr(0, valid_insert.length()-motif_size);                // perfix part of valid repeat
    string ext_motif = valid_insert.substr(valid_insert.length()-motif_size, motif_size);   // motif part of valid repeat

    uint valid_uplmt = ilen + remain_muts;
    string ins_repeat = utils::expand_repeat(ext_motif, valid_uplmt);
    
    int x = ins_repeat.length();        // the x-axis of distance matrix is valid insert repeat
    int y = insert.length();            // the y-axis of distance matrix is insert sequence
    int z = (x+1)*(y+1);                // size of the distance matrix
    int d_matrix[z];                    // distance matrix
    int i,j,substitution_cost;
    d_matrix[0+0*(y+1)] = 0;
    for ( i = 1; i <= x; i++ ) { d_matrix[i+0*(x+1)] = i; }
    for ( j = 1; j <= y; j++ ) { d_matrix[0+j*(x+1)] = j; }
    for ( j = 1; j <= y; j++ ) {
        for ( i = 1; i <= x; i++ ) {
            if ( ins_repeat[i-1] == insert[j-1] ) { substitution_cost = 0; }
            else { substitution_cost = 1; }
            d_matrix[i+j*(x+1)] = i4_min ( d_matrix[i-1+j*(x+1)] + 1,
                            i4_min ( d_matrix[i+(j-1)*(x+1)] + 1, 
                                    d_matrix[i-1+(j-1)*(x+1)] + substitution_cost ) );
        }
    }
    
    if (debug) {
        cout << "\n*** Levenshtein Matrix ***\n";
        for (int j=0; j<z; j++) {
            if (ins_repeat.length()>0 && j > 0) {
                if (j % (x+1) == 0) { cout << "\n"; }
            }
            cout << d_matrix[j] << "\t";
        }
        cout << "\n";
    }

    // if the sequence is an insertion
    if (!terminate) {

        // Calculate the edit distance of the insert from valid insert sequence
        int idx = z-1, min_idx = z-1;
        uint row_min = -1;
        for(; idx>=(x+1)*y; idx=idx-motif_size) {
            int mut = d_matrix[idx];
            if (mut < row_min) { row_min = mut; min_idx = idx; }
        }

        // If the minimum number of mutations is greater than the insert length
        // then minimum number of mutations equal the insert length
        if (row_min > ilen) { row_min = ilen; }

        // If the minimum mutations is greater than the remaining mutations
        // Terminate the repeat within allowed number of mutations
        if (row_min > remain_muts) {
            min_idx = z-1, row_min = -1;
            optimal_insert(x, y, rlen, muts, d_matrix, min_idx, row_min);
        }
        uint valid_insert_len = min_idx / (x+1);
        uint valid_repeat_len = min_idx % (x+1);
        if (debug) {
            cout << "*** Valid insert length: " << valid_insert_len << " ***\n";
            cout << "*** Valid repeat length: " << valid_repeat_len << " ***\n\n";
        }
        s = valid_insert_len;
        min_d = row_min;
        if (debug) {
            cout << "*** Insert:    " << s << " ***\n";
            cout << "*** Mutations: " << min_d << " ***\n";
        }
        if (s < ilen) {
            if (debug) { cout << "-xxx- Mid Termination -xxx-\n"; }
            terminate = 1;
        }
        else {
            if (debug) { cout << "-+++- Continuation -+++-\n"; }
            terminate = 0;
        }
        return { min_d, s, terminate };
    }
    
    // if the sequence is an extension
    else {
        int min_idx = z-1;
        uint row_min = -1;
        optimal_insert(x, y, rlen, muts, d_matrix, min_idx, row_min);
        uint valid_insert_len = min_idx / (x+1);
        uint valid_repeat_len = min_idx % (x+1);
        if (debug) {
            cout << "*** Valid insert length: " << valid_insert_len << " ***\n";
            cout << "*** Valid repeat length: " << valid_repeat_len << " ***\n\n";
        }
        s = valid_insert_len;
        min_d = row_min;
        
        if (debug) {
            cout << "*** Insert:    " << s << " ***\n";
            cout << "*** Mutations: " << min_d << " ***\n";
        }
        if (s == ilen && min_d < remain_muts) {
            if (debug) { cout << "-+++- Continuation -+++-\n"; }
            terminate = 0;
            uint c = valid_repeat_len % motif_size;
            return { min_d, s, terminate, c };
        }
        else {
            if (debug) { cout << "-xxx- Termination -xxx-\n"; }
            return { min_d, s, terminate };
        }
    }
}

int main(int argc, char* argv[]) {
    string fin = argv[1];
    string line;
    ifstream ins(fin);
    string seq_name;
    utils::SequenceWindow window(motif_size);

    while(getline(ins, line)) {
        if (line[0] == '>') {
            seq_name = line.substr(1, line.find(' ')-1);
            window.reset(motif_size);
        }
        else {
            for (const auto c: line) {
                window.update(toupper(c));

                if (window.count >= motif_size) {
                    uint position = window.count - motif_size;      // start of the motif
                    
                    // if (position >= 26250 && position <= 26280) {
                    //     debug = 1;
                    // }
                    // else {
                    //     debug = 0;
                    // }
                    
                    string curr_motif = window.motif;               // current motif
                    char curr_nuc = window.nuc;                     // current nucleotide i.e., last nucleotide of the motif
                    
                    if (debug) {
                        cout << "\n\n************************  Position: " << position << " Motif: " << curr_motif << "  ************************\n";
                    }

                    string curr_rclass = utils::get_repeat_class(curr_motif, rClassMap);    // repeat class of current motif
                    uint curr_rclass_first_check = 0;       // occurrence of repeat class for first time

                    // if repeat class is not encoutered so far
                    if (globalRepeatTracker.find(curr_rclass) == globalRepeatTracker.end()) {
                        curr_rclass_first_check = 1;
                        globalRepeatTracker[curr_rclass] = utils::RepeatTracker();
                        globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count);
                    }

                    vector<string> drop_rclasses;           // list of repeat classes that should be dropped after this iteration
                    std::unordered_map<string, utils::RepeatTracker>::iterator iter = globalRepeatTracker.begin();
                    std::unordered_map<string, utils::RepeatTracker>::iterator end_iter = globalRepeatTracker.end();
                    for(; iter != end_iter; ++iter) {
                        uint drop_rclass = 0;
                        string rclass = iter->first;

                        if (debug) { cout << "\n==========  " << rclass << "  ==========\n\n"; }
                        
                        if (rclass == curr_rclass) {
                            
                            // if repeat class is already encountered
                            if (!curr_rclass_first_check) {
                                uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                                string valid_motif = globalRepeatTracker[rclass].valid_motif;
                                
                                // either valid continuation or cyclical variation is found
                                // overlapping with the current continuation
                                if (position < globalRepeatTracker[rclass].end) {

                                    if (debug) {
                                        if (rclass_continue) {
                                            cout << "*** Valid continuation ***\n";
                                        }
                                    }
                                    uint d = 0;
                                    if (!rclass_continue) {
                                        uint overlap = globalRepeatTracker[rclass].end - position;
                                        uint c = 0; string cycle = "";
                                        for (; c < motif_size; c++) {
                                            cycle = valid_motif.substr(c) + valid_motif.substr(0, c);
                                            if (cycle == window.motif) { break; }
                                        }
                                        d = (overlap + c) % motif_size;
                                        if (debug) {
                                            cout << "*** Motif found overlapping ***\n";
                                            cout << "Motif: " << curr_motif << "\n";
                                            cout << "Valid: " << valid_motif << "\n";
                                            cout << "Overlap: " << overlap << "\n";
                                            cout << "Cycle: " << c << "\n";
                                            if (d > (motif_size-d)) { 
                                                cout << "Mutations: " << (motif_size - d) << "\n";
                                            }
                                            else {
                                                cout << "Mutations: " << d << "\n";
                                            }
                                            cout << "\n";
                                        }
                                        if (d > (motif_size-d)) { d = motif_size - d; }
                                    }

                                    uint start = globalRepeatTracker[rclass].start;
                                    uint end = globalRepeatTracker[rclass].end;
                                    uint rlen = end - start;
                                    uint muts = globalRepeatTracker[rclass].mutations;
                                    int remain_muts = ((window.count - start) * fraction_mutations) - muts;

                                    if (remain_muts < d) {
                                        uint terminate = 1;
                                        vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], valid_motif, terminate);
                                        uint d = insertion_result[0];
                                        uint s = insertion_result[1];
                                        globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, s);
                                        uint start = globalRepeatTracker[rclass].start;
                                        uint end = globalRepeatTracker[rclass].end;
                                        uint rlen = end - start;
                                        uint muts = globalRepeatTracker[rclass].mutations;
                                        string repeat = globalRepeatTracker[rclass].repeat;
                                        if (rlen >= 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass << "\t" << muts << "\t" << repeat << "\n"; 
                                        }

                                        // reinitialse the repeat from here
                                        globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count);
                                    }
                                    else {
                                        globalRepeatTracker[rclass].interrupt = 0;
                                        globalRepeatTracker[rclass].repeat += (globalRepeatTracker[rclass].insert + curr_nuc);
                                        globalRepeatTracker[rclass].insert = "";
                                        globalRepeatTracker[rclass].mutations += d;
                                        globalRepeatTracker[rclass].end = window.count;
                                        globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                    }

                                }

                                else if (position == globalRepeatTracker[rclass].end) {
                                    if (debug) {
                                        cout << "*** Motif found book ended ***\n";
                                    }
                                    globalRepeatTracker[rclass].interrupt = 0;
                                    globalRepeatTracker[rclass].repeat += (globalRepeatTracker[rclass].insert + curr_nuc);
                                    globalRepeatTracker[rclass].insert = "";
                                    globalRepeatTracker[rclass].end = window.count;
                                    globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                }

                                // valid continuation or cyclical variation is found
                                // after an interruption insertion
                                else {
                                    if (debug) {
                                        if (curr_motif == valid_motif) { cout << "*** Motif found after insertion ***\n"; }
                                        else { cout << "*** Cycle " << curr_motif << " found inplace of " << valid_motif << " ***\n"; }
                                    }

                                    // The last m-1 bp are removed from the insert
                                    // The m-1 bp are part of current valid continuation
                                    string insert = globalRepeatTracker[rclass].insert;
                                    globalRepeatTracker[rclass].insert = insert.substr(0, insert.length() - (motif_size-1));
                                    string cycle = ""; uint c = 0;
                                    // Check which cycle of the motif is found
                                    for (; c < motif_size; c++) {
                                        cycle = valid_motif.substr(c) + valid_motif.substr(0, c);
                                        if (cycle == window.motif) { break; }
                                    }
                                    string valid_insert = valid_motif.substr(0,c) + cycle;

                                    uint terminate = 0;
                                    vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], valid_insert, terminate);
                                    uint d = insertion_result[0];
                                    uint s = insertion_result[1];
                                    terminate = insertion_result[2];
                                    globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, s);
                                    if (terminate) {
                                        globalRepeatTracker[rclass].end += s;
                                        globalRepeatTracker[rclass].mutations += d;
                                        uint start = globalRepeatTracker[rclass].start;
                                        uint end = globalRepeatTracker[rclass].end;
                                        uint rlen = end - start;
                                        uint muts = globalRepeatTracker[rclass].mutations;
                                        string repeat = globalRepeatTracker[rclass].repeat;
                                        if (rlen >= 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass << "\t" << muts << "\t" << repeat << "\n"; 
                                        }

                                        // reinitialse the repeat from here
                                        globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count);
                                    }
                                    else {
                                        globalRepeatTracker[rclass].mutations += d;
                                        globalRepeatTracker[rclass].interrupt = 0;
                                        globalRepeatTracker[rclass].insert = "";
                                        globalRepeatTracker[rclass].end = window.count;
                                        globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                        globalRepeatTracker[rclass].repeat += curr_nuc;
                                    }
                                }
                            }
                        }
                        
                        else {
                            uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                            string valid_motif = globalRepeatTracker[rclass].valid_motif;

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
                                    if (debug) { cout << "*** Threshold insert length ***\n"; }
                                    uint terminate = 1;
                                    vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], valid_motif, terminate);
                                    uint d = insertion_result[0];
                                    uint s = insertion_result[1];
                                    terminate = insertion_result[2];
                                    globalRepeatTracker[rclass].end += s;
                                    globalRepeatTracker[rclass].mutations += d;
                                    globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, s);
                                    start = globalRepeatTracker[rclass].start;
                                    end = globalRepeatTracker[rclass].end;
                                    rlen = end - start;
                                    muts = globalRepeatTracker[rclass].mutations;
                                    string repeat = globalRepeatTracker[rclass].repeat;

                                    if (terminate) {
                                        if (rlen >= 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass << "\t" << muts << "\t" << repeat << "\n"; 
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
