//
//  kmerUtilities.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#include "kmerUtilities.hpp"
#include <climits>
#include <limits>

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

// struct Span {
//    uint32_t start;
//    uint32_t end;
// };

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
}

/* 
    SORT VECTORS FUNCTION
*/
void sort_positions(std::vector<store_position> &input_positions, uint32_t size)
    {
    std::sort(&input_positions[0], &input_positions[size],
              [](const store_position &a, const store_position &b) {
                  return (*a.position < *b.position);
              });
    }

/* 
    PROCESS DUPLICATE VECTORS FUNCTION
*/
void process_duplicates(std::vector<store_position> &input_positions){
    
    sort_positions(input_positions, input_positions.size());
    size_t duplicates = 1;
    size_t start = 0;
    for (std::size_t index = 1; index < input_positions.size(); index++)
        {
        if (*input_positions[index].position == *input_positions[index - 1].position) 
            duplicates++;
        else
            {
            if (duplicates != 1)
                {
                size_t size = input_positions[start].size;
                const uint32_t *ptr = input_positions[start].position;

                // Allocate space for the split lists
                for (size_t which = start; which < start + duplicates; which++)
                    {
                    // TO DO : Add to store_position the knowledge that this was dynamically allocated so that we can deallocate it later.
                    input_positions[which].position = new uint32_t[((size + duplicates - 1) / duplicates) + 1];
                    input_positions[which].size = 0;
                    }

                // Split the list n ways (round robin) where n is the number of duplicates
                for (size_t pos = 0; pos < size; pos++)  // -1 because we don't want the sentinal (which is added later)
                    {
                    store_position *which = &input_positions[start + (pos % duplicates)];
                    *((uint32_t *)(which->position) + (which->size++)) = ptr[pos];
                    }

                // Put the sentinal on the end
                for (size_t pos = start; pos < start + duplicates; pos++)
                    {
                    store_position *which = &input_positions[pos];
                    *((uint32_t *)(which->position) + (which->size++)) = UINT32_MAX;
                    }
                 }
            start = index;
            duplicates = 1;
            }
        }
 }

/* 
    VALIDATE VECTORS FUNCTION
*/
std::vector<std::pair<uint32_t, uint32_t>> validate_sets(std::vector<store_position> &positions, uint32_t min_matches, uint32_t query_length){
    std::vector<std::pair<uint32_t, uint32_t>> results_vector;
    uint32_t number_of_lists = positions.size();

    
    while (number_of_lists >= min_matches)
    {
        sort_positions(positions, number_of_lists);

        uint32_t firstValue = *positions[0].position;
        uint32_t minMatchesValue = *positions[min_matches - 1].position;

        if (minMatchesValue - firstValue <= query_length) 
        {
            for (uint32_t pos = 0; *positions[pos].position <= minMatchesValue; pos++)
            {
                while (*positions[pos].position <= minMatchesValue)
                    positions[pos].position++;
            
                if (*positions[pos].position == UINT32_MAX)
                    std::swap(positions[--number_of_lists], positions[pos]);
            }

            results_vector.emplace_back(firstValue, minMatchesValue);
        }
        else
        {

            for (uint32_t pos = 0; pos < min_matches; pos++)
            {
                while (*positions[pos].position < minMatchesValue)
                        positions[pos].position++;
                if (*positions[pos].position == UINT32_MAX)
                    std::swap(positions[--number_of_lists], positions[pos]); 
            }
        }
        
    } 
    return results_vector;
}

// example of printing out valid vectors
// for (const auto &sp : positions)
//     std::cout << "A.value: " << sp.ptr[sp.position] << std::endl;