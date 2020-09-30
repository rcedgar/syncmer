#include "myutils.h"
#include "kmer.h"
#include "randseq.h"

void TestSubmers_Old()
	{
	const uint L = 100;
	byte *Seq = myalloc(byte, L);
	MakeRandSeq(Seq, L);

	const uint k = 5;
	const uint m = 2;

	for (uint Pos = 0; Pos < 64; ++Pos)
		{
		uint64 Kmer = WordToKmer(Seq + Pos, k);
		Log("\n");
		Log("[%3u] %*.*s", Pos, k, k, Seq + Pos);
		Log("\n");
		for (uint j = 0; j < k; ++j)
			{
			byte SubWord_Rotate[64];
			uint64 Submer_Rotate = GetSubkmer_Rotate(Kmer, k, m, j);
			KmerToWord(Submer_Rotate, m, SubWord_Rotate);


			Log("   ");
			for (uint iii = 0; iii < j; ++iii)
				Log(" ");
			Log("%u  %*.*s", j, m, m, SubWord_Rotate);

			byte SubWord[64];
			if (j + m <= k)
				{
				uint64 Submer = GetSubkmer(Kmer, k, m, j);
				KmerToWord(Submer, m, SubWord);
				}
			else
				{
				memset(SubWord, '.', m);
				SubWord[m] = 0;
				}
			Log("  %s", SubWord);
			Log("  <%" PRIu64 ">", Submer_Rotate);
			Log("\n");
			}
		uint BestPos_Rotate = GetMinSubkmerPos_Rotate(Kmer, k, m);
		Log("  best pos %u\n", BestPos_Rotate);
		}
	}

void TestSubmers()
	{
	const uint L = 100;
	byte *Seq = myalloc(byte, L);
	MakeRandSeq(Seq, L);

	const uint k = 32;

	Log(" Seq=%*.*s\n", L, L, Seq);

#if 0
	for (uint Pos = 0; Pos + k <= L; ++Pos)
		{
		uint64 Kmer = WordToKmer(Seq + Pos, k);
		byte Word[65];
		KmerToWord(Kmer, k, Word);
		Word[k] = 0;
		Log("Word=");
		for (uint j = 0; j < Pos; ++j)
			Log(" ");
		Log("%s\n", Word);
		}
#endif

	uint64 Kmer = WordToKmer(Seq, k);
	uint m = 4;
	byte Word[65];
	KmerToWord(Kmer, k, Word);
	Word[k] = 0;
	Log("Word=");
	Log("%s\n", Word);
	for (uint Pos = 0; Pos < 4; ++Pos)
		{
		uint64 Sub = GetSubkmer(Kmer, k, m, Pos);
		byte SubWord[64];
		KmerToWord(Sub, m, SubWord);
		SubWord[m] = 0;
		Log("Sub %s\n", SubWord);
		}
	}
