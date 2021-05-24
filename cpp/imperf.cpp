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
const bool debug = 1;

/*
 *  Indentifies the optimal length of insert within allowed mutations by looking
 *  at the distance matrix of insert and valid repeat sequence
 *  @param x length of the valid repeat sequence (x-axis of distance matrix)
 *  @param y length of the insert sequence (y-axis of distance matrix)
 *  @param rlen length of the repeat so far
 *  @param muts number of mutations in the repeat so far
 *  @param fraction_mutations float (< 1) allowed mutations as fraction of repeat length
 *  @param distance_matrix the distance matrix of insert versus repeat sequence
 *  @param min_idx index in the distance matrix with optimal insert length and mutations (by ref)
 *  @param row_least_mutations least mutations corresponding to the repeat length (by ref)
 *  @return minimum of the two integers
*/
void optimal_insert(int x, int y, uint rlen, uint muts, float fraction_mutations, int distance_matrix[], int &min_idx, uint &row_least_muts) {
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


/*
 *  Inspects the insert sequence and identifies optimal length of insert sequence
 *  within allowed mutations
 *  @param rtracker RepeatTracker object of a particular repeat
 *  @param m motif size
 *  @param valid_insert valid repeat sequence against which insert is compared
 *  @param terminate indicates if the insert is at the terminal end of the repeat
 *  @return vector of insert length, least mutations, signal for repeat termination
*/
vector<uint> insertion_mutations(utils::RepeatTracker rtracker, string vir_motif, uint terminal, uint m, float fraction_mutations) {

    string valid_motif = rtracker.valid_motif;    
    uint start = rtracker.start;
    uint end   = rtracker.end;
    uint rlen  = end - start;
    uint muts  = rtracker.mutations;

    string insert = rtracker.insert;
    uint ilen = insert.length();

    // Final minimum mutations and insert length
    uint least_muts = -1, final_ilen = 0, terminate = 0;

    // if repeat already reached threshold mutations
    if (muts == rlen * fraction_mutations) {
        least_muts = 0, final_ilen = 0; terminate = 1;
        return { least_muts, final_ilen, terminate };
    }

    uint threshold_muts = fraction_mutations * (rlen + ilen);
    int remain_muts  = threshold_muts - muts;

    if (debug) {
        cout << "*** Insertion: " << insert << " ***";
        cout << "\n*** Valid: " << vir_motif << " ***";
        cout << "\n*** Remaining mutations: " << remain_muts << " ***\n";
    }

    if (remain_muts <= 0) {
        least_muts = 0; final_ilen = 0; terminate = 1;
        return { least_muts, final_ilen, terminate };
    }

    // virm - valid insert repeat motif
    uint virm_len = vir_motif.length();
    string vir_prefix   = vir_motif.substr(0, virm_len - m);     // perfix part of valid repeat
    string vir_extmotif = vir_motif.substr(virm_len - m, m);     // motif part of valid repeat

    uint vir_maxlen = ilen + remain_muts;
    string vi_repeat = utils::expand_repeat(vir_extmotif, vir_maxlen);
    
    int x = vi_repeat.length();         // the x-axis of distance matrix is valid insert repeat
    int y = insert.length();            // the y-axis of distance matrix is insert
    int z = (x+1)*(y+1);                // size of the distance matrix
    int distance_matrix[z];             // distance matrix
    int sub_cost;                       // substitution cost
    int i, j;
    distance_matrix[0] = 0;
    for ( i = 1; i <= x; i++ ) { distance_matrix[i+0*(x+1)] = i; }
    for ( j = 1; j <= y; j++ ) { distance_matrix[0+j*(x+1)] = j; }
    for ( j = 1; j <= y; j++ ) {
        for ( i = 1; i <= x; i++ ) {
            if ( vi_repeat[i-1] == insert[j-1] ) { sub_cost = 0; }
            else { sub_cost = 1; }
            distance_matrix[i+j*(x+1)] = utils::min_int (
                                            distance_matrix[i-1+j*(x+1)] + 1,
                                            utils::min_int (
                                                distance_matrix[i+(j-1)*(x+1)] + 1,
                                                distance_matrix[i-1+(j-1)*(x+1)] + sub_cost
                                            )
                                         );
        }
    }
    
    if (debug) {
        cout << "\n*** Levenshtein Matrix ***\n";
        for (int j=0; j<z; j++) {
            if (vi_repeat.length()>0 && j > 0) {
                if (j % (x+1) == 0) { cout << "\n"; }
            }
            cout << distance_matrix[j] << "\t";
        }
        cout << "\n";
    }

    // if the sequence is an insertion
    if (!terminal) {

        // Calculate the edit distance of the insert from valid repeat sequence
        int idx = z-1, min_idx = z-1;
        for(; idx>=(x+1)*y; idx=idx-motif_size) {
            int mut = distance_matrix[idx];
            if (mut < least_muts) { least_muts = mut; min_idx = idx; }
        }

        // If the minimum number of mutations is greater than the insert length
        // then minimum number of mutations equal the insert length
        if (least_muts > ilen) { least_muts = ilen; }

        // If the minimum mutations is greater than the remaining mutations
        // Terminate the repeat within allowed number of mutations
        if (least_muts > remain_muts) {
            min_idx = z-1, least_muts = -1;
            optimal_insert(x, y, rlen, muts, fraction_mutations, distance_matrix, min_idx, least_muts);
        }
        uint vir_len = min_idx % (x+1);
        final_ilen = min_idx / (x+1);
        if (debug) {
            cout << "*** Insert:    " << final_ilen << " ***\n";
            cout << "*** Valid repeat length: " << vir_len << " ***\n";
            cout << "*** Mutations: " << least_muts << " ***\n";
        }
        if (final_ilen < ilen) {
            if (debug) { cout << "-xxx- Mid Termination -xxx-\n"; }
            terminate = 1;
        }
        else {
            if (debug) { cout << "-+++- Continuation -+++-\n"; }
            terminate = 0;
        }
        return { least_muts, final_ilen, terminate };
    }
    
    // if the sequence is an extension
    else {
        int min_idx = z-1;
        optimal_insert(x, y, rlen, muts, fraction_mutations, distance_matrix, min_idx, least_muts);
        uint vir_len = min_idx % (x+1);
        final_ilen = min_idx / (x+1);
        if (debug) {
            cout << "*** Insert:    " << final_ilen << " ***\n";
            cout << "*** Valid repeat length: " << vir_len << " ***\n";
            cout << "*** Mutations: " << least_muts << " ***\n";
        }
        // if (final_ilen == ilen && least_muts < remain_muts) {
        //     if (debug) { cout << "-+++- Continuation -+++-\n"; }
        //     terminate = 0;
        //     uint c = vir_len % motif_size;
        //     return { least_muts, final_ilen, terminate, c };
        // }
        // else {
            if (debug) { cout << "-xxx- Termination -xxx-\n"; }
            terminate = 1;
            return { least_muts, final_ilen, terminate };
        // }
    }
}


/*
 *  Handle the termination of sequence, report all valid repeats
 *  @param seq_name name/id of the sequence
 *  @return none prints out all occurrences of valid repeats
*/
void sequence_termination(string seq_name) {
    vector<string> drop_rclasses;
    std::unordered_map<string, utils::RepeatTracker>::iterator iter = globalRepeatTracker.begin();
    std::unordered_map<string, utils::RepeatTracker>::iterator end_iter = globalRepeatTracker.end();
    for(; iter != end_iter; ++iter) {
        string rclass = iter->first;
        string valid_motif = globalRepeatTracker[rclass].valid_motif;
        uint terminal = 1;
        vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], valid_motif, terminal, motif_size, fraction_mutations);
        uint least_muts = insertion_result[0];
        uint final_ilen = insertion_result[1];
        uint terminate = insertion_result[2];
        globalRepeatTracker[rclass].end += final_ilen;
        globalRepeatTracker[rclass].mutations += least_muts;
        globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, final_ilen);
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
        drop_rclasses.push_back(rclass);
    }
    for (uint j=0; j<drop_rclasses.size(); j++) {
        globalRepeatTracker.erase(drop_rclasses[j]);
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
            // terminate all current repeats
            sequence_termination(seq_name);

            // new sequence initiation
            seq_name = line.substr(1, line.find(' ')-1);
            window.reset(motif_size);
        }
        else {
            for (const auto c: line) {
                window.update(toupper(c));

                if (window.count >= motif_size) {
                    uint position = window.count - motif_size;      // start of the motif
                    string curr_motif = window.motif;               // current motif
                    char curr_nuc = window.nuc;                     // current nucleotide i.e., last nucleotide of the motif
                    string curr_rclass = "";
                    uint curr_rclass_first_check = 1;

                    if (debug) {
                        cout << "\n\n************************  Position: " << position << " Motif: " << curr_motif << "  ************************\n";
                    }

                    if (curr_nuc == 'N') {
                        curr_rclass = "N";
                        curr_rclass_first_check = 0;
                    }
                    else {
                        curr_rclass = utils::get_repeat_class(curr_motif, rClassMap);    // repeat class of current motif
                        uint curr_rclass_first_check = 0;       // occurrence of repeat class for first time

                        // if repeat class is not encoutered so far
                        if (globalRepeatTracker.find(curr_rclass) == globalRepeatTracker.end()) {
                            curr_rclass_first_check = 1;
                            globalRepeatTracker[curr_rclass] = utils::RepeatTracker();
                            globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count);
                        }
                    }

                    vector<string> drop_rclasses;           // list of repeat classes that should be dropped after this iteration
                    
                    // Update all repeat sequences with current window
                    std::unordered_map<string, utils::RepeatTracker>::iterator iter = globalRepeatTracker.begin();
                    std::unordered_map<string, utils::RepeatTracker>::iterator end_iter = globalRepeatTracker.end();
                    for(; iter != end_iter; ++iter) {
                        string rclass = iter->first;
                        uint drop_rclass = 0;

                        if (debug) { cout << "\n==========  " << rclass << "  ==========\n\n"; }
                        
                        // updating currently encountered repeat class
                        if (rclass == curr_rclass) {
                            
                            // if repeat class is already encountered
                            if (!curr_rclass_first_check) {
                                uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                                string valid_motif   = globalRepeatTracker[rclass].valid_motif;
                                
                                // either valid continuation or cyclical variation is found
                                // overlapping with the current continuation
                                if (position < globalRepeatTracker[rclass].end) {

                                    if (debug) {
                                        if (rclass_continue) { cout << "*** Valid continuation ***\n"; }
                                    }

                                    uint cycle_muts = 0;        // mutations when overlapping cyclical variant is found
                                    if (!rclass_continue) {
                                        uint overlap = globalRepeatTracker[rclass].end - position;
                                        uint c = 0; string cycle = "";
                                        for (; c < motif_size; c++) {
                                            cycle = valid_motif.substr(c) + valid_motif.substr(0, c);
                                            if (cycle == window.motif) { break; }
                                        }
                                        cycle_muts = (overlap + c) % motif_size;
                                        if (cycle_muts > (motif_size - cycle_muts) ) {
                                            cycle_muts = motif_size - cycle_muts;
                                        }
                                        if (debug) {
                                            cout << "*** Motif found overlapping ***\n";
                                            cout << "Motif: " << curr_motif << "\n";
                                            cout << "Valid: " << valid_motif << "\n";
                                            cout << "Overlap: " << overlap << "\n";
                                            cout << "Cycle: " << c << "\n";
                                            cout << "Mutations: " << cycle_muts << "\n\n";
                                        }
                                    }

                                    uint start = globalRepeatTracker[rclass].start;
                                    uint end = globalRepeatTracker[rclass].end;
                                    uint rlen = end - start;
                                    uint muts = globalRepeatTracker[rclass].mutations;
                                    int remain_muts = ((window.count - start) * fraction_mutations) - muts;

                                    // if mutations introduced by cyclical variation is greater than allowed
                                    // remaining mutations
                                    if (cycle_muts > remain_muts) {
                                        uint terminal = 1;
                                        vector<uint> result_vector = insertion_mutations(globalRepeatTracker[rclass], valid_motif, terminal, motif_size, fraction_mutations);
                                        uint least_muts = result_vector[0];
                                        uint final_ilen = result_vector[1];
                                        globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, final_ilen);
                                        globalRepeatTracker[rclass].mutations += least_muts;
                                        globalRepeatTracker[rclass].end += final_ilen;
                                        start = globalRepeatTracker[rclass].start;
                                        end = globalRepeatTracker[rclass].end;
                                        rlen = end - start;
                                        muts = globalRepeatTracker[rclass].mutations;
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
                                        globalRepeatTracker[rclass].mutations += cycle_muts;
                                        globalRepeatTracker[rclass].end = window.count;
                                        globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                    }

                                }

                                // if a cyclical variation is found at the end of the previous motif
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
                                // after an interruption by insertion
                                else {
                                    if (debug) {
                                        if (curr_motif == valid_motif) { cout << "*** Motif found after insertion ***\n"; }
                                        else { cout << "*** Cycle " << curr_motif << " found inplace of " << valid_motif << " ***\n"; }
                                    }

                                    // last m-1 bp are removed from the insert as they are part of current valid continuation
                                    string insert = globalRepeatTracker[rclass].insert;
                                    globalRepeatTracker[rclass].insert = insert.substr(0, insert.length() - (motif_size-1));
                                    string cycle = ""; uint c = 0;
                                    // Check which cycle of the motif is found
                                    for (; c < motif_size; c++) {
                                        cycle = valid_motif.substr(c) + valid_motif.substr(0, c);
                                        if (cycle == window.motif) { break; }
                                    }
                                    string vir_motif = valid_motif.substr(0,c) + cycle;     // valid insert repeat motif
                                    uint terminal = 0;      // non terminal insertion
                                    vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], vir_motif, terminal, motif_size, fraction_mutations);
                                    uint least_muts = insertion_result[0];
                                    uint final_ilen = insertion_result[1];
                                    uint terminate = insertion_result[2];
                                    globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, final_ilen);
                                    globalRepeatTracker[rclass].mutations += least_muts;

                                    if (terminate) {
                                        globalRepeatTracker[rclass].end += final_ilen;
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
                                        globalRepeatTracker[rclass].insert = "";
                                        globalRepeatTracker[rclass].end = window.count;
                                        globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                        globalRepeatTracker[rclass].repeat += curr_nuc;
                                    }
                                }
                            }
                        }
                        
                        // updating all other repeat classes
                        else {
                            uint rclass_continue = !(globalRepeatTracker[rclass].interrupt);
                            string valid_motif = globalRepeatTracker[rclass].valid_motif;
                            uint start = globalRepeatTracker[rclass].start;
                            uint end = globalRepeatTracker[rclass].end;
                            uint rlen = end - start;
                            uint muts = globalRepeatTracker[rclass].mutations;

                            if (rclass_continue) {
                                // if mutations already reached allowed number
                                if (muts == (rlen * fraction_mutations)) {
                                    string repeat = globalRepeatTracker[rclass].repeat;
                                    if (rlen >= 12) {
                                        if (debug) { cout << "*** Valid repeat ***\n"; }
                                        cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                        rclass << "\t" << muts << "\t" << repeat << "\n"; 
                                    } 
                                    drop_rclass = 1;
                                }

                                else {
                                    globalRepeatTracker[rclass].interrupt = 1;
                                    globalRepeatTracker[rclass].insert = curr_nuc;
                                    globalRepeatTracker[rclass].valid_motif = valid_motif[motif_size-1] + valid_motif.substr(0, motif_size-1);
                                }
                            }

                            else {
                                globalRepeatTracker[rclass].insert += curr_nuc;
                                string insert = globalRepeatTracker[rclass].insert;
                                uint ilen = insert.length();

                                // if no cyclical variation is found for greater than 4 motif lengths
                                // or for greater than current repeat length
                                if ((ilen >= 5*motif_size) || (ilen >= rlen)) {
                                    if (debug) { cout << "*** Threshold insert length ***\n"; }
                                    uint terminal = 1;
                                    vector<uint> insertion_result = insertion_mutations(globalRepeatTracker[rclass], valid_motif, terminal, motif_size, fraction_mutations);
                                    uint least_muts = insertion_result[0];
                                    uint final_ilen = insertion_result[1];
                                    uint terminate = insertion_result[2];
                                    globalRepeatTracker[rclass].end += final_ilen;
                                    globalRepeatTracker[rclass].mutations += least_muts;
                                    globalRepeatTracker[rclass].repeat += globalRepeatTracker[rclass].insert.substr(0, final_ilen);
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

                    // discontinuing all dropped repeat class
                    for (uint j=0; j<drop_rclasses.size(); j++) {
                        globalRepeatTracker.erase(drop_rclasses[j]);
                    }
                }
            }
        }
    }
    
    
    // End of the file
    // terminate all current repeats
    sequence_termination(seq_name);

    ins.close();
}
