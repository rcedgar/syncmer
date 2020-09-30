#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

uint GetCount(const map<uint, uint> &Map, uint n);
uint64 GetTotalCount(const map<uint, uint> &Map);

static byte *Seq;
static FILE *fTab;

static void Spacing(SyncmerType ST, uint k, uint t)
	{
	if (ST == ST_Syncmer2 && t + 1 >= k)
		return;

	uint w;
	if (optset_w)
		w = opt(w);
	else
		w = SyncmerIndex::EstimateWindow(k, t);

	SyncmerIndex SI;
	SI.Create(ST, k, t, w, Seq, BENCHL);

	map<uint, uint> SpaceToCount;
	uint Max = SI.GetSpaceToCount(SpaceToCount);
	uint64 Total = GetTotalCount(SpaceToCount);
	uint SyncmerCount = SI.GetSyncmerCount();
	uint KmerCount = SI.GetKmerCount();
	double MaxFreq = 1e-6;
	for (unsigned Space = 1; Space <= Max; ++Space)
		{
		uint n = GetCount(SpaceToCount, Space);
		double Freq = double(n)/double(Total);
		MaxFreq = max(MaxFreq, Freq);
		}

	double Compress = double(KmerCount)/double(SyncmerCount);

	fprintf(fTab, "k=%u", k);
	if (t != UINT_MAX)
		fprintf(fTab, ",t=%u", t);
	if (w != UINT_MAX)
		fprintf(fTab, ",w=%u", w);
	fprintf(fTab, ",Submers=%u", SyncmerCount);
	fprintf(fTab, ",Kmers=%u", KmerCount);
	fprintf(fTab, ",Compress=%.1f", Compress);
	fprintf(fTab, ",max_space=%u", Max);
	fprintf(fTab, "\n");
	
	const uint H = 32;
	for (unsigned Space = 1; Space <= Max; ++Space)
		{
		uint n = GetCount(SpaceToCount, Space);
		double Freq = double(n)/double(Total);
		uint h = uint((Freq*H)/MaxFreq);
		fprintf(fTab, "%u", Space);
		fprintf(fTab, "\t%.4f", Freq);
		fprintf(fTab, "\t");
		for (uint i = 0; i < h; ++i)
			fprintf(fTab, "=");
		fprintf(fTab, "\n");
		}

	if (opt(log_syncmers))
		SI.LogSyncmers();
	}

void cmd_spacings()
	{
	const string Name = opt(spacings);
	asserta(optset_tabbedout);
	fTab = CreateStdioFile(opt(tabbedout));

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);

	SyncmerType ST = StrToST(Name);

	if (opt(tile))
		{
		for (uint k = 8; k <= 20; k += 2)
			{
			Progress("k=%u\n", k);
			Spacing(ST, k, k-1);
			}
		}
	else if (optset_k && optset_w)
		{
		Spacing(ST, opt(k), UINT_MAX);
		}
	else if (optset_klo)
		{
		asserta(optset_khi);
		asserta(optset_tlo);
		asserta(optset_thi);
		for (uint k = opt(klo); k <= opt(khi); ++k)
			{
			for (uint t = opt(tlo); t <= opt(thi); ++t)
				{
				Progress("k=%u,t=%u\n", k, t);
				Spacing(ST, k, t);
				}
			}
		}
	else
		{
		for (uint k = 8; k <= 20; k += 2)
			{
			Progress("k=%u\n", k);
			for (uint t = 2; t < k; ++t)
				Spacing(ST, k, t);
			}
		}

	CloseStdioFile(fTab);
	}
