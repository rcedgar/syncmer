#include "myutils.h"
#include "kmer.h"
#include "randseq.h"
#include "syncmerindex.h"

void TestMinimizers()
	{
	const uint L = 100;
	byte *Seq = myalloc(byte, L);
	MakeRandSeq(Seq, L);

	const uint k = 5;
	const uint w = 3;
	const uint t = 0;

	//SyncmerIndex &SI = SyncmerIndex::Create(ST_Minimizer1, k, t, w, Seq, L);

	for (uint Pos = 20; Pos < 21; ++Pos)
		{
		uint64 Kmer = WordToKmer(Seq + Pos, k);
		Log("\n");
		Log("[%3u] %*.*s", Pos, k, k, Seq + Pos);
		Log("\n");
		for (uint j = 0; j < k; ++j)
			{
			}
		}
	}
