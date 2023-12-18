//
//  kmerUtilities.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#ifndef kmerUtilities_hpp
#define kmerUtilities_hpp

#include <stdio.h>
#include <sys/stat.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
#include <cstdint>
#include <string>
#include <set>

struct Span {
    uint32_t start;
    uint32_t end;
};

struct Position {
    uint32_t value;
    size_t inputSetIndex;
};

uint64_t packKmer(const char *sequence);
std::string unpackKmer(uint64_t packed_sequence);
std::vector<uint32_t> findLongestSequenceWithinRange(std::vector<uint32_t>& nums, int range);
char *read_entire_file(const char *filename, uint64_t& fileSize);
std::vector<std::vector<uint32_t>> findAllSequencesWithinRange(std::vector<uint32_t>& nums,int range);
uint32_t calculateCoverage(const std::vector<uint32_t>& pos, uint32_t span);
std::vector<uint32_t> findLongestSubsequence(const std::vector<uint32_t>& pos, int N);
uint32_t getSpan(const std::vector<uint32_t>& pos, uint32_t span);
std::vector<uint32_t> cleanupVector(const std::vector<uint32_t>& pos, uint32_t length);
std::vector<std::pair<uint32_t, uint32_t>> validSpans(const std::vector<std::set<uint32_t>>& inputSets, uint32_t N, uint32_t S);
uint32_t getSpan(const std::vector<Position>& pos, uint32_t span);
std::vector<std::pair<uint32_t, uint32_t>> validate_sets(std::vector<std::set<uint32_t>>& input_sets, uint32_t min_matches, uint32_t query_length);
#endif /* kmerUtilities_hpp */

