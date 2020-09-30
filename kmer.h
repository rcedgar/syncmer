#ifndef kmer_h
#define kmer_h

#include "murmur.h"

uint64 WordToKmer(const byte *Seq, uint k);
const byte *KmerToWord(uint64 Kmer, uint k, byte *Word);

uint GetMinSubkmerPos(uint64 Kmer, uint k, uint m);
uint GetMinSubkmerPos_Hash(uint64 Kmer, uint k, uint m);
uint GetMinSubkmerPos_Rotate(uint64 Kmer, uint k, uint m);

uint64 GetSubkmer(uint64 Kmer, uint k, uint m, uint i);
uint64 GetSubkmer_Rotate(uint64 Kmer, uint k, uint m, uint i);

#endif // kmer_h
