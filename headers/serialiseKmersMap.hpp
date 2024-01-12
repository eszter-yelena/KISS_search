//
//  serialiseKmersMap.hpp
//  indexReference
//
//  Created by Shlomo Geva on 19/7/2023.
//

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <sstream>

#include "kmerUtilities.hpp"

struct QueryTable {
    std::map<uint32_t, std::string> entries;

    uint32_t query(uint32_t queryKey) const {
        auto it = entries.lower_bound(queryKey);
        if (it != entries.begin()) {
            --it; // Move the iterator one step back to the highest entry not exceeding the queryKey
            return (uint32_t) std::distance(entries.begin(), it); // Return the index
        } else {
            if (!entries.empty()) {
                // If the queryKey is smaller than the first entry, return 0 (index of the first entry)
                return 0;
            }
        }
        return -1; // Return -1 to indicate key not found for an empty table
    }
 
    // Function to print out the keys of the 'entries' map
    void printKeys() const {
        std::cout << "Table entries" << std::endl;
        int n=0;
        for (const auto& entry : entries) {
            std::cout << n++ << ":  " << entry.first << " "
                      << entry.second << std::endl;
        }
    }
    // Function to get the key from the index
    uint32_t getKeyFromIndex(int index) const {
        if (index >= 0 && index < entries.size()) {
            auto it = entries.begin();
            std::advance(it, index);
            return it->first;
        } else {
            // Invalid index, return some default value or throw an exception
            return 0;
        }
    }
    // Function to access the elements of the underlying map
    const std::string& getValue(uint32_t key) const {
        // Use the at() method of the map to access the value
        // It throws an exception if the key is not found
        return entries.at(key);
    }

};

//struct QueryTable {
//    std::map<uint32_t, std::string> entries;
//
//    // Update the function to return the index instead of the key
//    int query(uint32_t queryKey) const {
//        auto it = entries.lower_bound(queryKey);
//        if (it != entries.end()) {
//            return std::distance(entries.begin(), it); // Return the index
//        } else {
//            // Check if queryKey exceeds the last key in the map
//            auto lastEntryIt = entries.end();
//            if (lastEntryIt != entries.begin()) {
//                --lastEntryIt; // Get the iterator to the last entry
//                if (queryKey > lastEntryIt->first) {
//                    return std::distance(entries.begin(), lastEntryIt); // Return the index of the last entry
//                }
//            }
//            return -1; // Return -1 to indicate key not found
//        }
//    }
//};

QueryTable readTableFromFile(const std::string& filename);

void serializeMap(const std::vector<std::vector<uint32_t>>& kmersMap, const std::string& innerMapFilename, const std::string& outerMapFilename);

void deserializeMap(const std::string& innerMapFilename, const std::string& outerMapFilename, std::vector<uint32_t>& innerMapBlob, std::vector<uint32_t>& outerMapBlob);

void getInnerVector(std::vector<store_position> &result, const std::vector<uint32_t> &innerMapBlob, const std::vector<uint32_t> &outerMapBlob, uint32_t index);

bool writeTextBlobToFile(const char* text, std::size_t length, const std::string& filename);

std::string readTextBlobFromFile(const std::string& filename);

