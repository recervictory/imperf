/*
    Utils: Auxiliary functions for imperf.cpp
    @file utils.h
    @author Akshay Avvaru
    @version 0.1 06/08/2020
*/
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <bitset>
#include <fstream>
#include <chrono>
#include <string>
#inlcude "utils.h"

using namespace std;
using namespace chrono;
using namespace utils;

namespace mut_utils {
    
    /*
    *  Inspects the insert sequence and identifies optimal length of insert sequence
    *  within allowed mutations
    *  @param rtracker RepeatTracker object of a particular repeat
    *  @param m motif size
    *  @param vir_motif valid repeat sequence against which insert is compared
    *  @param terminal indicates if the insert is terminal
    *  @param fraction_mutations threshold number of mutations as a fraction of repeat length
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

        // if repeat already reached threshold mutations and the insert 
        // doesn't increase the threshold number of mutations
        if ((muts == rlen * fraction_mutations) && (ilen * fraction_mutations < 1)) {
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

        // As we have already checked if the number of mutations are at the threshold level
        // and if the insert length allows addition of new mutataions the remaining mutations
        // should always be greater than 0.
        // if (remain_muts <= 0) {
        //     least_muts = 0; final_ilen = 0; terminate = 1;
        //     return { least_muts, final_ilen, terminate };
        // }

        // virm - valid insert repeat motif
        uint virm_len = vir_motif.length();
        string vir_prefix   = vir_motif.substr(0, virm_len - m);     // perfix part of valid repeat
        string vir_extmotif = vir_motif.substr(virm_len - m, m);     // motif part of valid repeat

        // insert is compared with the maximum allowed length of the repeat
        // the valid repeat sequence should also contain the prefix sequence
        uint vir_maxlen = ((ilen + remain_muts)/m) + 1;
        string vi_repeat = vir_prefix;
        for (uint v=0; v<vir_maxlen; v++) { vi_repeat += vir_extmotif; }
        
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
            
            if (debug) { cout << "-xxx- Termination -xxx-\n"; }
            terminate = 1;
            return { least_muts, final_ilen, terminate };
                
        }
    }

    vector<uint> backtrack_repeat(utils::RepeatTracker rtracker, float fraction_mutations, uint m) {
        string valid_motif = ""
        for (int i=m-1; i>=0; i--) {
            valid_motif += rtracker.repeat[i];
        }

        uint start = rtracker.start;
        uint end   = rtracker.end;
        uint rlen  = end - start;
        uint muts  = rtracker.mutations;
        string upstream = rtracker.upstream;
        uint ulen = upstream.length();

        // Final minimum mutations and upstream length
        uint least_muts = -1, final_ulen = 0;

        uint threshold_muts = fraction_mutations * (rlen + ilen);
        int remain_muts  = threshold_muts - muts;

        if (debug) {
            cout << "*** Upstream: " << upstream << " ***";
            cout << "\n*** Valid: " << vir_motif << " ***";
            cout << "\n*** Remaining mutations: " << remain_muts << " ***\n";
        }

        // insert is compared with the maximum allowed length of the repeat
        // the valid repeat sequence should also contain the prefix sequence
        uint vir_maxlen = ulen + remain_muts;
        string vi_repeat = utils::expand_repeat(valid_motif, vir_maxlen);
        
        int x = vi_repeat.length();         // the x-axis of distance matrix is valid insert repeat
        int y = upstream.length();            // the y-axis of distance matrix is insert
        int z = (x+1)*(y+1);                // size of the distance matrix
        int distance_matrix[z];             // distance matrix
        int sub_cost;                       // substitution cost
        int i, j;
        distance_matrix[0] = 0;
        for ( i = 1; i <= x; i++ ) { distance_matrix[i+0*(x+1)] = i; }
        for ( j = 1; j <= y; j++ ) { distance_matrix[0+j*(x+1)] = j; }
        for ( j = 1; j <= y; j++ ) {
            for ( i = 1; i <= x; i++ ) {
                if ( vi_repeat[i-1] == upstream[j-1] ) { sub_cost = 0; }
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

        int min_idx = z-1;
        optimal_insert(x, y, rlen, muts, fraction_mutations, distance_matrix, min_idx, least_muts);
        uint vir_len = min_idx % (x+1);
        final_ulen = min_idx / (x+1);
        if (debug) {
            cout << "*** Upstream:    " << final_ulen << " ***\n";
            cout << "*** Valid repeat length: " << vir_len << " ***\n";
            cout << "*** Mutations: " << least_muts << " ***\n";
        }
        
        if (debug) { cout << "-xxx- Termination -xxx-\n"; }
        return { least_muts, final_ulen };
    }

    
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

}