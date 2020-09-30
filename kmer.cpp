#include "myutils.h"
#include "kmer.h"
#include "alpha.h"

uint64 WordToKmer(const byte *Seq, uint k)
	{
	asserta(k <= 32);

	uint64 Kmer = 0;
	for (uint i = 0; i < k; ++i)
		{
		byte c = Seq[i];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter > 3)
			return UINT64_MAX;
		Kmer = (Kmer << uint64(2)) | Letter;
		}
	return Kmer;
	}

const byte *KmerToWord(uint64 Kmer, uint k, byte *Word)
	{
	asserta(k <= 32);

	Word[k] = 0;
	for (uint i = 0; i < k; ++i)
		{
		byte Letter = Kmer & 3;
		byte c = g_LetterToCharNucleo[Letter];
		Word[k-i-1] = c;
		Kmer >>= uint64(2);
		}
	return Word;
	}

uint GetMinSubkmerPos_Rotate(uint64 Kmer, uint k, uint s)
	{
	asserta(k <= 32);
	assert(k >= s);
	uint BestPos = 0;
	uint64 BestMer = GetSubkmer_Rotate(Kmer, k, s, 0);
	for (uint i = 1; i <= k; ++ i)
		{
		uint64 Mer = GetSubkmer_Rotate(Kmer, k, s, i);
		if (Mer < BestMer)
			{
			BestMer = Mer;
			BestPos = i;
			}
		}
	return BestPos;
	}

uint GetMinSubkmerPos_Hash(uint64 Kmer, uint k, uint s)
	{
	asserta(k <= 32);
	assert(k >= s);
	uint BestPos = 0;
	uint64 BestMer = murmur64(GetSubkmer(Kmer, k, s, 0));
	for (uint Pos = 1; Pos + s <= k; ++Pos)
		{
		uint64 Mer = murmur64(GetSubkmer(Kmer, k, s, Pos));
		if (Mer < BestMer)
			{
			BestMer = Mer;
			BestPos = Pos;
			}
		}
	return BestPos;
	}

uint GetMinSubkmerPos(uint64 Kmer, uint k, uint s)
	{
	asserta(k <= 32);
	assert(k >= s);
	uint BestPos = 0;
	uint64 BestMer = GetSubkmer(Kmer, k, s, 0);
	for (uint Pos = 1; Pos + s <= k; ++Pos)
		{
		uint64 Mer = GetSubkmer(Kmer, k, s, Pos);
		if (Mer < BestMer)
			{
			BestMer = Mer;
			BestPos = Pos;
			}
		}
	return BestPos;
	}

uint64 GetSubkmer_Rotate(uint64 Kmer, uint k, uint s, uint Pos)
	{
	asserta(k <= 32);
	asserta(s <= k);
	uint64 Mer = 0;
	for (uint j = 0; j < s; ++j)
		{
		uint64 Shift = k - (Pos + s)%k + j;
		Shift = Shift%k;
		byte Letter = byte(((Kmer >> (2*Shift)) & 3));
		Mer = Mer | (Letter << 2*j);
		}
	return Mer;
	}

uint64 GetSubkmer(uint64 Kmer, uint k, uint s, uint Pos)
	{
	asserta(k <= 32);
	asserta(Pos + s <= k);
	uint64 Mer = 0;
	Kmer >>= 2*(k - Pos - s);
	for (uint j = 0; j < s; ++j)
		{
		byte Letter = byte(((Kmer >> 2*j) & 3));
		Mer = Mer | (Letter << 2*j);
		}
	return Mer;
	}
