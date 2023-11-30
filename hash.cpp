//
//  hash.cpp
//  indexKISS
//
//  Created by Shlomo Geva on 13/7/2023.
//

#include "hash.hpp"

/*
    REVERSE COMPLEMENT (uint64)
 */
uint64_t reverse_complement(uint64_t kmer, size_t bases)
    {
        // in-situ reverese complement of a kmer packed into uint64
        uint64_t result = 0;
        /*
            Compute the complement (A<->T, C<->G), which is a bit-fit because A=00, T=11, C=01, G=10
        */
        uint64_t complement = ~kmer;

        /*
            Reverse the order of the bases
        */
        for (uint64_t base = 0; base < bases; base++, complement >>= 2)
            result = (result << 2) | (complement & 0x03);

        return result;
    }

/*
    MURMURHASH3
 */
uint32_t murmurHash3(uint64_t key) {
    // hash a 64bit value to 32 bits.
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return static_cast<uint32_t>(key);
}

/*
    XOR_HASH
*/
uint32_t xorHash(uint64_t packedKmer) {
    uint32_t hash = static_cast<uint32_t>(packedKmer >> 32) ^ static_cast<uint32_t>(packedKmer);
    return hash;
}

/*
     CANONICAL_KMER
 */
uint64_t canonical_kmer(uint64_t kmer,int KMERSIZE) {
    return kmer^reverse_complement(kmer,KMERSIZE);
}

