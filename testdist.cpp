#include "myutils.h"
#include "kmer.h"
#include "randseq.h"
#include "syncmerindex.h"

void TestDist()
	{
	const uint L = 10000;
	byte *Seq = myalloc(byte, L);
	MakeRandSeq(Seq, L);

	asserta(optset_k);
	asserta(optset_t);
	const uint k = opt(k);
	const uint t = opt(t);
//	const uint m = ::GetSubmerLength(k, t);
	const uint m = 2;
	Log("k %u, t %u, m %u, 1/t=%.4f\n", k, t, m, 1.0/t);
	vector<uint> Counts(t);
	const uint K = L - k + 1;
	for (uint Pos = 20; Pos < K; ++Pos)
		{
		uint64 Kmer = WordToKmer(Seq + Pos, k);
		uint MaxPos = GetMinSubkmerPos_Rotate(Kmer, k, m);
//		Log("%u\n", MaxPos);
		asserta(MaxPos < t);
		++(Counts[MaxPos]);
		}

	Log("Dist:\n");
	for (uint i = 0; i < t; ++i)
		{
		double f = double(Counts[i])/K;
		Log("%u	%.4f\n", i, f);
		}
	}
