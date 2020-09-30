#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

uint64 GetTotalCount(const vector<uint> &Vec)
	{
	uint64 Sum = 0;
	for (vector<uint>::const_iterator p = Vec.begin(); p != Vec.end(); ++p)
		{
		uint Count = *p;
		Sum += uint64(Count);
		}
	return Sum;
	}

static byte *Seq;
static FILE *fTab;

static void DepthDist(SyncmerType ST, uint k, uint t)
	{
	asserta(t < k);

	uint w = SyncmerIndex::EstimateWindow(k, t);
	SyncmerIndex SI;
	SI.Create(ST, k, t, w, Seq, BENCHL);

	vector<uint> DepthToCount;
	uint Max = SI.GetDepthToCount(DepthToCount);
	asserta(Max <= 255);
	uint64 Total = GetTotalCount(DepthToCount);
	fprintf(fTab, "k=%u,t=%u,w=%u\n", k, t, w);
	for (unsigned Depth = 0; Depth <= w + 2; ++Depth)
		{
		uint n = DepthToCount[Depth];
		double Freq = double(n)/double(Total);
		fprintf(fTab, "%u", Depth);
		fprintf(fTab, "\t%.4f", Freq);
		fprintf(fTab, "\n");
		}
	}

void cmd_depthdist()
	{
	const string Name = opt(depthdist);
	asserta(optset_tabbedout);
	fTab = CreateStdioFile(opt(tabbedout));

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);

	SyncmerType ST = StrToST(Name);

	if (opt(tile))
		{
		for (uint k = 8; k <= 16; ++k)
			{
			Progress("k=%u\n", k);
			DepthDist(ST, k, k-1);
			}
		}
	else
		{
		for (uint k = 8; k <= 16; ++k)
			{
			Progress("k=%u\n", k);
			for (uint t = 2; t <= 10; ++t)
				DepthDist(ST, k, t);
			}
		}

	CloseStdioFile(fTab);
	}
