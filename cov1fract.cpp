#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

static byte *Seq;
static FILE *fTab;

static void Cov1Fract(const vector<SyncmerType> STs, uint k, uint t)
	{
	asserta(t < k);

	uint w = SyncmerIndex::EstimateWindow(k, t);

	fprintf(fTab, "k=%u	t=%u	w=%u", k, t, w);
	for (uint i = 0; i < SIZE(STs); ++i)
		{
		SyncmerIndex SI;
		SyncmerType ST = STs[i];
		if (!SyncmerIndex::ValidParams(ST, k, t, w))
			{
			fprintf(fTab, "	%s=invalid_params\n", STToStr(ST));
			continue;
			}
		SI.Create(ST, k, t, w, Seq, BENCHL);

		double Fract = SI.GetCov1Fract();
		fprintf(fTab, "	%s=%.4f", STToStr(ST), Fract);
		}
	fprintf(fTab, "\n");
	}

void cmd_cov1fract()
	{
	const string Names = opt(cov1fract);
	asserta(optset_tabbedout);
	fTab = CreateStdioFile(opt(tabbedout));
	asserta(!optset_tile);

	vector<string> NameVec;
	Split(Names, NameVec, '+');

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);

	vector<SyncmerType> STs;
	for (uint i = 0; i < SIZE(NameVec); ++i)
		{
		const string &Name = NameVec[i];
		SyncmerType ST = StrToST(Name);
		STs.push_back(ST);
		}

	for (uint k = 8; k <= 20; ++k)
		{
		Progress("k=%u\n", k);
		Cov1Fract(STs, k, k-1);
		}

	CloseStdioFile(fTab);
	}
