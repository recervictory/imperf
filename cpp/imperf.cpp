#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "utils.h"
#include "mut_utils.h"

using namespace std;
using namespace utils;
using namespace mut_utils;

// global variables
unordered_map<string, string> rClassMap;
unordered_map<string, utils::RepeatTracker> globalRepeatTracker;
// running individually for each motif is faster than
// running than in parllel
const uint motif_size = 4;
const float fraction_mutations = 0.1;
const bool debug = 0;

/*
 * In this version the repeat is extend on the left side by backtracking. 
*/



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
        uint atomicity = utils::check_atomicity(rclass);
        if (rlen >= 12) {
            if (debug) { cout << "*** Valid repeat ***\n"; }
            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
            rclass.substr(0, atomicity) << "\t" << muts << "\t" << repeat << "\n"; 
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
                    // uint curr_rclass_first_check = 1;

                    if (debug) {
                        cout << "\n\n************************  Position: " << position << " Motif: " << curr_motif << "  ************************\n";
                    }

                    
                    curr_rclass = utils::get_repeat_class(curr_motif, rClassMap);    // repeat class of current motif
                    uint curr_rclass_first_check = 0;       // occurrence of repeat class for first time

                    // if repeat class is not encoutered so far
                    if (globalRepeatTracker.find(curr_rclass) == globalRepeatTracker.end()) {
                        curr_rclass_first_check = 1;
                        globalRepeatTracker[curr_rclass] = utils::RepeatTracker();
                        globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count, window.upstream);
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
                                if (position <= globalRepeatTracker[rclass].end) {

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
                                        uint atomicity = utils::check_atomicity(rclass);
                                        if (rlen >= 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass.substr(0, atomicity) << "\t" << muts << "\t" << repeat << "\n"; 
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
                                        uint atomicity = utils::check_atomicity(rclass);
                                        if (rlen >= 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass.substr(0, atomicity) << "\t" << muts << "\t" << repeat << "\n"; 
                                        }

                                        // reinitialse the repeat from here
                                        globalRepeatTracker[curr_rclass].initialise(window.motif, position, window.count);
                                    }
                                    else {
                                        globalRepeatTracker[rclass].interrupt = 0;
                                        globalRepeatTracker[rclass].insert = "";
                                        globalRepeatTracker[rclass].end = window.count;
                                        globalRepeatTracker[rclass].valid_motif = curr_motif.substr(1) + curr_motif[0];
                                        globalRepeatTracker[rclass].repeat += curr_motif;
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
                                    uint atomicity = utils::check_atomicity(rclass);
                                    if (rlen >= 12) {
                                        if (debug) { cout << "*** Valid repeat ***\n"; }
                                        cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                        rclass.substr(0, atomicity) << "\t" << muts << "\t" << repeat << "\n"; 
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
                                        uint atomicity = utils::check_atomicity(rclass);
                                        if (rlen >= 12) {
                                            if (debug) { cout << "*** Valid repeat ***\n"; }
                                            cout << seq_name << "\t" << start << "\t" << end << "\t" << rlen << "\t" << 
                                            rclass.substr(0, atomicity) << "\t" << muts << "\t" << repeat << "\n"; 
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
