//
//  hash.hpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#ifndef hash_hpp
#define hash_hpp

#include <stdio.h>
#include <cstdint>

uint32_t murmurHash3(uint64_t key);
uint32_t xorHash(uint64_t packedKmer);
uint64_t canonical_kmer(uint64_t kmer,int KMERSIZE);
uint64_t reverse_complement(uint64_t kmer, size_t bases);

#endif /* hash_hpp */


////
////  hash.hpp
////  indexKISS
////
////  Created by Shlomo Geva on 13/7/2023.
////
//
//#ifndef hash_hpp
//#define hash_hpp
//
//#include <stdio.h>
//#include <cstdint>
//
//uint32_t murmurHash3(uint64_t key);
//uint32_t xorHash(uint64_t packedKmer);
//
//#endif /* hash_hpp */

