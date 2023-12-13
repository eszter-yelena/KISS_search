//
//  smithWaterman.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 15/7/2023.
//

#ifndef smithWaterman_hpp
#define smithWaterman_hpp

#include <stdio.h>
#include <string>
#include <tuple>
#include <map>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>
#include <climits>

//int smithWaterman(const std::string& seq1, const std::string& seq2, int MATCH_SCORE, int MISMATCH_SCORE, int GAPSCORE, int BAND_WIDTH);

std::string reverseComplement(const std::string& sequence);

//int smithWatermanScore(const std::string& sequence1, const std::string& sequence2,
//          int matchScore, int mismatchScore, int gapScore);

std::vector<std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>>
    getSW(std::string genome, std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>>& candidate_results,
          int matchScore, int mismatchScore, int gapScore, int cutoff);

uint64_t reverse_complement(uint64_t kmer, size_t bases);

struct SortByScore {
    bool operator()(const std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>& a, const std::tuple<uint32_t, uint32_t, uint32_t, double, std::string>& b) const {
        if (std::get<3>(a) != std::get<3>(b)) {
            return std::get<3>(a) > std::get<3>(b); // Sort by swScore (third value)
        } else if (std::get<0>(a) != std::get<0>(b)) {
            return std::get<0>(a) < std::get<0>(b); // Sort by first value (uint32_t)
        } else {
            return std::get<1>(a) < std::get<1>(b); // Sort by second value (uint32_t)
        }
    }
};
//struct SortByScore {
//    bool operator()(const std::tuple<int, int, int>& a, const std::tuple<int, int, int>& b) const {
//        if (std::get<2>(a) != std::get<2>(b)) {
//            return std::get<2>(a) > std::get<2>(b); // Sort by swScore (third value)
//        } else {
//            return std::get<0>(a) < std::get<0>(b); // Sort by sampleIndex (first value)
//        }
//    }
//};

#endif /* smithWaterman_hpp */

