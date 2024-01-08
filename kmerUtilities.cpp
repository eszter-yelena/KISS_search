//
//  kmerUtilities.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#include "kmerUtilities.hpp"
#include <climits>

extern unsigned int KMERSIZE;
extern int CUTOFF;
extern uint64_t kmer_encoding_table[];

std::vector<uint32_t> findLongestSubsequence(const std::vector<uint32_t>& pos, int N) {
    uint32_t start = 0;
    uint32_t maxLength = 0;
    uint32_t maxSpan = 0;
    uint32_t endIndex = 0;

    for (int end = 0; end < pos.size(); ++end) {
        while (pos[end] - pos[start] > N) {
            ++start;
        }
        
        if (end - start + 1 > maxLength) {
            maxLength = end - start + 1;
            maxSpan = pos[end] - pos[start];
            endIndex = end;
        } else if (end - start + 1 == maxLength) {
            int currentSpan = pos[end] - pos[start];
            if (currentSpan > maxSpan) {
                maxSpan = currentSpan;
                endIndex = end;
            }
        }
    }

    std::vector<uint32_t> longestSubsequence;
    for (uint32_t i = endIndex - maxLength + 1; i <= endIndex; ++i) {
        longestSubsequence.push_back(pos[i]);
    }

    return longestSubsequence;
}


std::vector<uint32_t> cleanupVector(const std::vector<uint32_t>& pos, uint32_t length) {
    // walk along the hits and delete all positions that cannot be part of a match
    //    because they are too far from any other hit position
    std::vector<uint32_t> cleanedPos;
    int start=-1;
    if (pos.empty()) {
        return cleanedPos; // Return an empty vector
    }

    for (int i=1; i<pos.size(); i++ )
    {
        // skip while next is too far ahead
        if (pos[i]-pos[i-1]<=length) {
            start=i-1;
            break; // found two positions within
        }
    }
    if (start == -1 || start == pos.size()) {
        return cleanedPos; // Return an empty vector
    }
        
    cleanedPos.push_back(pos[start]);

    for (size_t i = start+1; i < pos.size(); ++i) {
        if (pos[i] - cleanedPos.back() <= length) {
            cleanedPos.push_back(pos[i]);
        } else if (i + 1 < pos.size() && pos[i + 1] - pos[i] <= length) {
            cleanedPos.push_back(pos[i]);
        }
    }

    return cleanedPos;
}

//struct Span {
//    uint32_t start;
//    uint32_t end;
//};

uint32_t getSpan(const std::vector<uint32_t>& pos, uint32_t span) {
    // calculate the span of matches covered by seeds
    if (pos.empty()) {
        return 0;
    }

    std::vector<Span> spans;

    for (const uint32_t p : pos) {
        spans.push_back({p, p + span - 1});
    }

    std::vector<Span> mergedSpans;
    std::sort(spans.begin(), spans.end(), [](const Span& a, const Span& b) {
        return a.start < b.start;
    });

    for (const Span& span : spans) {
        if (mergedSpans.empty() || span.start > mergedSpans.back().end) {
            mergedSpans.push_back(span);
        } else {
            mergedSpans.back().end = std::max(mergedSpans.back().end, span.end);
        }
    }

    uint32_t totalSpan = 0;
    for (const Span& span : mergedSpans) {
        totalSpan += (span.end - span.start + 1);
    }

    return totalSpan;
}


/*
    PACK KMER
    -----------
*/
uint64_t packKmer(const char *sequence)
    {
    uint64_t packed = 0;

    for (int pos = 0; pos < KMERSIZE; pos++)
        packed = (packed << 2) | kmer_encoding_table[(size_t)sequence[pos]];

    return packed;
    }

/*
    UNPACK KMER
    -------------
*/
std::string unpackKmer(uint64_t packed_sequence)
    {
    std::string sequence;
    for (int32_t pos = KMERSIZE-1; pos >= 0; pos--)
        {
        switch ((packed_sequence >> (pos * 2)) & 3)
            {
            case 0:
                sequence += 'A';
                break;
            case 1:
                sequence += 'C';
                break;
            case 2:
                sequence += 'G';
                break;
            case 3:
                sequence += 'T';
                break;
            }
        }
    return sequence;
    }

std::vector<uint32_t> findLongestSequenceWithinRange(std::vector<uint32_t>& nums, int range) {
//    std::vector<uint32_t> sortedNums = nums;
//    std::sort(sortedNums.begin(), sortedNums.end());

    int start = 0;
    int maxLength = 0;// number of elements in max length sequence
    int maxStart = 0;// position where max length sequence begins

    for (int end = 1; end < nums.size(); ++end) {
        while (nums[end] - nums[start] > range)
        {
            // skip while the next is too far from start (beyond range).
            ++start;
        }
        // need to maintain the max span sequence
        if (end - start + 1 > maxLength) {
            maxLength = end - start + 1;
            maxStart = start;
        }
    }
    // Check if the longest sequence is valid (more than one element)
    if (maxLength <= 1) {
        return {}; // Return an empty vector to indicate no valid subsequence was found
    }

    // Extract the longest sequence from the sorted vector
    std::vector<uint32_t> longestSequence(nums.begin() + maxStart, nums.begin() + maxStart + maxLength);
    return longestSequence;
}

std::vector<std::vector<uint32_t>> findAllSequencesWithinRange(std::vector<uint32_t>& nums, int range) {
    std::vector<std::vector<uint32_t>> result;
//    std::sort(nums.begin(), nums.end());

    while (true) {
        std::vector<uint32_t> longestSequence = findLongestSequenceWithinRange(nums, range);
        if (longestSequence.empty()) {
            break; // No valid subsequence found, exit the loop
        }
        if (getSpan(longestSequence,range)>=CUTOFF)
            result.push_back(longestSequence);

        // Remove the longest sequence from the original sequence to find the next subsequences
        auto it = std::search(nums.begin(), nums.end(), longestSequence.begin(), longestSequence.end());
        if (it == nums.end()) {
            break; // Longest sequence not found in nums, break the loop
        }

        nums.erase(it, it + longestSequence.size());
    }

    return result;
}
/*
    READ_ENTIRE_FILE
 */
char *read_entire_file(const char *filename, uint64_t& fileSize)
  {
  FILE *fp;
  struct stat details;
  char *contents = NULL;
  fileSize = 0;
  if ((fp = fopen(filename, "rb")) != NULL)
  {
      if (fstat(fileno(fp), &details) == 0)
      {
          if (details.st_size != 0 || details.st_size > UINT32_MAX)
          {
              contents = (char *)malloc(details.st_size);
              if (fread(contents, details.st_size, 1, fp) != 1)
              {
                  free(contents);
                  contents = NULL;
              }
              else
              {
                  fileSize = details.st_size;
              }
          }
      }
      fclose(fp);
   }
   return contents;
}

uint32_t calculateCoverage(const std::vector<uint32_t>& pos, uint32_t span) {
    // this function calculates the number of bases covered by the matching seeds.
    // seeds can cover overlapping parts of the reference
    if (pos.empty()) {
        return 0;
    }

    uint32_t coverage = 0;
    uint32_t start = pos[0];
    uint32_t end = pos[0] + span;

    for (size_t i = 1; i < pos.size(); ++i) {
        if (pos[i] <= end) {
            // Current range overlaps with the previous one, extend the coverage.
            end = std::max(end, pos[i] + span);
        } else {
            // Non-overlapping range found, update the coverage and start a new range.
            coverage += end - start;
            start = pos[i];
            end = pos[i] + span;
        }
    }

    // Add the coverage of the last range.
    coverage += end - start;

    return coverage;
}


void displayProgress(uint64_t current, uint64_t total, int desiredUpdateInterval) {
    uint64_t percent = (current * 100) / total;
    static uint64_t lastDisplayedPercent = -desiredUpdateInterval; // Initialize to a value that will trigger the first update
    if (percent - lastDisplayedPercent >= desiredUpdateInterval) {
        std::cout << "Progress: " << std::setw(3) << percent << "%" << std::endl;
        std::cout.flush();
        lastDisplayedPercent = percent;
    }
    /*
     // example of use:
     int totalIterations = 1000;
     int percentUpdateInterval = 10; // Update progress every 10%
     
     for (int i = 0; i < totalIterations; ++i) {
     // Your loop processing here
     
     displayProgress(i + 1, totalIterations, percentUpdateInterval);
     }
     */
}

//struct Position {
//    uint32_t value;
//    size_t inputSetIndex;
//};

/*
    isVALID
 */
bool isValid(const std::vector<int>& counts, uint32_t minMatches) {
    size_t uniqueSets = 0;
    for (int count : counts) {
        if (count > 0) {
            uniqueSets++;
        }
    }
    return uniqueSets >= minMatches;
}

/*
    VALID SPANS
 */
std::vector<std::pair<uint32_t, uint32_t>> validSpans(const std::vector<std::set<uint32_t>>& inputSets, uint32_t minMatches, uint32_t lastPos) {
    
    std::vector<std::pair<uint32_t, uint32_t>> result;
    size_t numSets = inputSets.size();


    // Create a vector to store positions along with their set index
    std::vector<Position> positions;
    for (size_t i = 0; i < numSets; ++i) {
        for (uint32_t pos : inputSets[i]) {
            positions.push_back({pos, i});
           // std:: cout << "positions: " << pos << " : " << i << "\n";
        }
    }

    // Sort positions based on their values
    std::sort(positions.begin(), positions.end(), [](const Position& a, const Position& b) {
        return a.value < b.value;
    });

    // for (const auto &p : positions ){
    //     std:: cout << "positions: " << p.value << ": " << p.inputSetIndex << "\n";
    // }

    // Initialize counts for each set to 0
    std::vector<int> counts(numSets, 0);
    size_t startPos = 0;
    size_t endPos = 0;

    // Iterate over the sorted positions
    while (startPos < positions.size()) {
        if (endPos < positions.size() && positions[endPos].value - positions[startPos].value <= lastPos) {
            counts[positions[endPos].inputSetIndex]++;
            // Check if the current configuration of counts satisfies the condition N
            if (isValid(counts, minMatches)) {        
                // If valid, collect the start and end positions of the range
                if (positions[endPos].value-positions[startPos].value+KMERSIZE >= CUTOFF) {
                    // TODO: check coverage against cutoff
                    std::vector<Position> pos(endPos-startPos+1);
                    for (size_t k=startPos; k<=endPos; k++)
                        pos[k-startPos] = positions[k];
                   
                    if (getSpan(pos,KMERSIZE)>=CUTOFF)
                        // collecting first and last position in range (for SW alignment)
                        result.push_back({positions[startPos].value, positions[endPos].value});
                }

                // Skip the entire sequence
                while (startPos < positions.size() && positions[startPos].value <= positions[endPos].value) {
                    counts[positions[startPos].inputSetIndex]--;
                    startPos++;
                }
            }
            endPos++;
        } else {
            counts[positions[startPos].inputSetIndex]--;
            startPos++;
        }
    }

// for debugging checking results are the same
// std::cout << "Results: ";
// for (const auto& span : result) {
//     std::cout << "[" << span.first << ", " << span.second << "] ";
// }
// std:: cout << "\n";

    return result;
}

/*
    GET SPAN
 */
uint32_t getSpan(const std::vector<Position>& pos, uint32_t span) {
    if (pos.empty()) {
        return 0;
    }

    std::vector<Span> spans;

    for (const Position p : pos) {
        spans.push_back({p.value, p.value + span - 1});
    }

    std::vector<Span> mergedSpans;
    std::sort(spans.begin(), spans.end(), [](const Span& a, const Span& b) {
        return a.start < b.start;
    });

    for (const Span& span : spans) {
        if (mergedSpans.empty() || span.start > mergedSpans.back().end) {
            mergedSpans.push_back(span);
        } else {
            mergedSpans.back().end = std::max(mergedSpans.back().end, span.end);
        }
    }

    uint32_t totalSpan = 0;
    for (const Span& span : mergedSpans) {
        totalSpan += (span.end - span.start + 1);
    }

    return totalSpan;
}

/* 
    SORT VECTORS FUNCTION
*/

// Function to sort sets based on their first elements
std::vector<std::set<uint32_t>> sort_sets(const std::vector<std::set<uint32_t>>& input_sets) {
    std::vector<std::set<uint32_t>> sorted_sets;
    
    for (const auto& s : input_sets) {
        if (!s.empty()) {
            sorted_sets.push_back(s);
        }
    }

    std::sort(sorted_sets.begin(), sorted_sets.end(), [](const std::set<uint32_t>& a, const std::set<uint32_t>& b) {
        return *a.begin() < *b.begin();
    });

    
    return sorted_sets;
}

std::vector<std::pair<uint32_t, uint32_t>> validate_sets(std::vector<std::set<uint32_t>>& input_sets, uint32_t min_matches, uint32_t query_length) {
    std::vector<std::pair<uint32_t, uint32_t>> results_vector;

    if (input_sets.size() >= min_matches) {
     
        // sort the sets based on the first element
        std::vector<std::set<uint32_t>> sorted_sets = sort_sets(input_sets);
        
        uint32_t end_span_value = 0;

        // while the value of the sorted list at min_matches is not empty
        while (sorted_sets.size() >= min_matches && !sorted_sets[min_matches - 1].empty()) {  
            sorted_sets.erase(std::remove_if(sorted_sets.begin(), sorted_sets.end(), [](const std::set<uint32_t>& s) { return s.empty(); }), sorted_sets.end());

            // sorted_sets.erase(std::remove_if(sorted_sets.begin(), sorted_sets.end(), [](const std::set<uint32_t>& s) { return s.empty(); }), sorted_sets.end());
            uint32_t lastValue = *sorted_sets.back().begin();
            uint32_t minmatchesValue = *std::next(sorted_sets.begin(), min_matches-1)->begin();
            uint32_t firstValue = *sorted_sets.front().begin();
            
            // if the first value is less than the previous largest value remove that set
            if (firstValue < end_span_value){
                for (auto& set : sorted_sets) {
                    set.erase(set.begin());
                }
            sorted_sets.erase(std::remove_if(sorted_sets.begin(), sorted_sets.end(), [](const std::set<uint32_t>& s) { return s.empty(); }), sorted_sets.end());
            }

            // check if last minus first is in range and add it to the results, and remove the range once it has been added
             if (minmatchesValue - firstValue <= query_length ) {
                results_vector.emplace_back(firstValue, lastValue);
                
                for (auto& set : sorted_sets) {
                    set.erase(set.begin());
                }
            }

            // else remove the first value of the first ordered set and review
            else {
                sorted_sets.front().erase(sorted_sets.front().begin());
            }
            //set the end value and resort the sets
            end_span_value = minmatchesValue;
            sorted_sets = sort_sets(sorted_sets);
        }
    }

    return results_vector;
}