//
//  smithWaterman.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 15/7/2023.
//

#include "smithWaterman.hpp"
#include "ssw_cpp.h"
#include <sstream>
/*
     COMPUTE SCORE AND LENGTH (from CIGAR)
 */

extern unsigned int KMERSIZE;
extern bool MATCH_ALL; // find all mathces flag
extern std::vector<std::string> sampleSequences;

std::pair<int, int> computeScoreAndLength(const std::string& cigarString) {
    /*
     The SAM spec offers us this table of CIGAR operations which indicates which ones "consume" the query or the reference, complete with explicit instructions on how to calculate sequence length from a CIGAR string:
                                                                  Consumes  Consumes
     Op  BAM Description                                             query  reference
     M   0   alignment match (can be a sequence match or mismatch)   yes   yes
     I   1   insertion to the reference                              yes   no
     D   2   deletion from the reference                             no    yes
     N   3   skipped region from the reference                       no    yes
     S   4   soft clipping (clipped sequences present in SEQ)        yes   no
     H   5   hard clipping (clipped sequences NOT present in SEQ)    no    no
     P   6   padding (silent deletion from padded reference)         no    no
     =   7   sequence match                                          yes   yes
     X   8   sequence mismatch                                       yes   yes
     “Consumes query” and “consumes reference” indicate whether the CIGAR operation causes the alignment to step along the query sequence and the reference sequence respectively.
     ...
     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
     */
    int score = 0;
    int length = 0;
    std::stringstream ss(cigarString);
    char c;
    int currentNumber = 0;

    while (ss >> currentNumber >> c) {
        // this could be a simple if statement, but using a switch to allow
        //     different penalty/reward for alignment operations in score calculations
        switch (c) {
            case '=':
                score += currentNumber;
                length += currentNumber; // Add length for the matches to the total length
                break;
            case 'M':
            case 'I':
            case 'S':
            case 'X':
            case 'D':
            case 'N':
            case 'H':
            case 'P':
                score -= currentNumber; // penalise all variations due to misalignment
            default:
                break;
        }
    }
    return {score, length};
}
 
/*
     SSW Aligner
 */
void alignSequences(const std::string ref, const std::string query,std::string& cigar, int& sw_score, uint32_t& refBegin, uint32_t& refEnd) {
    int32_t maskLen = (int32_t) strlen(query.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;

    // Declares a default Aligner
    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the query to the ref
    aligner.Align(query.c_str(), ref.c_str(), (int) ref.size(), filter, &alignment, maskLen);

/*
// the following is an example from the source code of the aligner.
// it is here only to document the different values that may be extracted.
//     std::cout << "===== SSW result =====" << std::endl;
//     std::cout << "Best Smith-Waterman score:\t" << alignment.sw_score << std::endl
//         << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << std::endl
//         << "Reference start:\t" << alignment.ref_begin << std::endl
//         << "Reference end:\t" << alignment.ref_end << std::endl
//         << "Query start:\t" << alignment.query_begin << std::endl
//         << "Query end:\t" << alignment.query_end << std::endl
//         << "Next-best reference end:\t" << alignment.ref_end_next_best << std::endl
//         << "Number of mismatches:\t" << alignment.mismatches << std::endl
//         << "Cigar: " << alignment.cigar_string << std::endl;
//     std::cout << "======================" << std::endl;
*/
    sw_score = alignment.sw_score;
    cigar = alignment.cigar_string;
    refBegin = alignment.ref_begin;
    refEnd = alignment.ref_end;
}

int max(int a, int b, int c) {
    return std::max(std::max(a, b), c);
}

/*
    REVERSE COMPLEMENT (string)
 */
std::string reverseComplement(const std::string& sequence) {
    std::string reverseSeq;
    for (char nucleotide : sequence) {
        switch (nucleotide) {
            case 'A':
                reverseSeq += 'T';
                break;
            case 'T':
                reverseSeq += 'A';
                break;
            case 'C':
                reverseSeq += 'G';
                break;
            case 'G':
                reverseSeq += 'C';
                break;
            default:
                reverseSeq += nucleotide;
        }
    }
    std::reverse(reverseSeq.begin(), reverseSeq.end());
    return reverseSeq;
}

/*
     GETSW
 */

// Function to compute the Smith-Waterman scores between pairs in the results
// In this code, the resultsSW vector is defined as std::vector<std::tuple<uint32, uint32, uint32>>
// to store the index pair and the Smith-Waterman score. The assignment of sampleIndex
// and referenceIndex is done directly from the result tuple, and the Smith-Waterman score
// is computed using the sample and reference sequences. Finally, the index pair
// and the SW score are stored in the resultsSW vector using std::make_tuple.

std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>> getSW(std::string genome, std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>>& candidate_results, int matchScore, int mismatchScore, int gapScore, int cutoff)
{
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>> localResultsSW;

    // Iterate through each result in the results vector

//    std::cout << "getSW over " << candidate_results.size() << " candidate results, " << std::endl;

    std::string cigar;
    uint32_t currentSample = UINT_MAX;// sentinel
    bool skipMulti = false; // used to skip multiple matches per read when MATCH_ALL is false
    
    for (const auto& result : candidate_results) {
        uint32_t sampleIndex, referenceIndex, minPos, maxPos;
        std::tie(sampleIndex, referenceIndex, minPos, maxPos) = result;

// for debugging, print out the result
// std::cout << "Candidate Result: (" << sampleIndex << ", " << referenceIndex << ", " << minPos << ", " << maxPos << ")" << std::endl;

        
        if (currentSample == sampleIndex) {
            if (skipMulti) {
                continue;
            }
        } else {
            currentSample = sampleIndex;
        }
        
        // Retrieve the sample and reference sequences
        std::string sampleSequence = sampleSequences[sampleIndex];

        // this code is doing alignment on a limited region in the refrence
        // comment it out for a full alignment with the reference (beware long ref sequences)
        // TODO: should be able to work out from the positions of seeds in which direction to align
        uint32_t expand = 50;
        minPos = (minPos <= expand) ? 0 : (minPos - expand);
        maxPos = std::min(maxPos+expand,(uint32_t) genome.size());
        std::string referenceSequence = genome.substr(minPos,maxPos-minPos+1);
        int swScore=0;
        int swLength=0;
        uint32_t refBegin;
        uint32_t refEnd;

        // use SSW lib code to compute SW
        // match in forward direction
        std::pair<int, int> res;
        alignSequences(referenceSequence, sampleSequence, cigar, swScore, refBegin, refEnd);
        res = computeScoreAndLength(cigar);
        // swScore = res.first; // note - we replace the swScore from alignSequences() with that derived from cigar
        if (swScore < cutoff) {
            // match in reverse direction
            alignSequences(referenceSequence, reverseComplement(sampleSequence), cigar, swScore, refBegin, refEnd);
            res = computeScoreAndLength(cigar); // note - we replace the swScore from alignSequences()
            swScore = res.first;
        }
        swLength  = res.second;
        refBegin += minPos; // adjust to genome position
        refEnd += minPos;   //    - " -
 
        // Store the result index pair and the SW score in the resultsSW vector

        if (swScore >= cutoff) {
            double score = ((double) swScore)/swLength;
            localResultsSW.push_back(std::make_tuple(sampleIndex, refBegin, refEnd, score, cigar));
            if (!MATCH_ALL) {
                skipMulti = true;
                currentSample = sampleIndex;
            }
        }
    }
    
//    std::cout << "Passed SW " << localResultsSW.size() << std::endl;

    return localResultsSW;
}
