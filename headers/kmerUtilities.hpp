//
//  kmerUtilities.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//
#pragma once

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

struct store_position {
    const uint32_t* ptr;
    std::size_t position;
    std::size_t size;
    uint32_t hash;
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
std::vector<std::pair<uint32_t, uint32_t>> validate_sets(std::vector<store_position>& input_sets, uint32_t min_matches, uint32_t query_length);
void process_duplicates(std::vector<store_position> &input_positions);