#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

extern uint BENCHL;

static FILE *fConsTab;
static FILE *fWinTab;
static byte *Seq;
static uint PctIds[] = { 99 , 98, 97, 95, 90, 85, 80 };
static const uint IDS = sizeof(PctIds)/sizeof(PctIds[0]);
static byte *MutatedSeqs[IDS];
static uint TestCount;
static uint WindowSize = 64;

static uint ConsBetterCount1;
static uint ConsBetterCount2;
static uint ConsTieCount;

static uint CovBetterCount1;
static uint CovBetterCount2;
static uint CovTieCount;

static uint SizeBetterCount1;
static uint SizeBetterCount2;
static uint SizeTieCount;

static uint WinBetterCount1;
static uint WinBetterCount2;
static uint WinTieCount;

static vector<double> CovDeltas;

static vector<double> ConsDeltas;
static vector<double> SizeDeltas;

static vector<double> ConsDeltasPctIds[IDS];

void Bench1(uint k, uint t, uint w, const string &Name)
	{
	//if (t + 1 >= k)
	//	return;

	SyncmerType ST = StrToST(Name);

	SyncmerIndex SI;
	SI.Create(ST, k, t, w, Seq, BENCHL);

	vector<SyncmerIndex> SIs(IDS);
	vector<SyncmerIndex> SI2s(IDS);
	for (uint i = 0; i < IDS; ++i)
		{
		const byte *MutatedSeq = MutatedSeqs[i];
		SIs[i].Create(ST, k, t, w, MutatedSeq, BENCHL);
		}

	uint Size = SI.GetIndexSize();

	double t1 = 1.0/SI.GetFractKmersIndexed();
	double Cov = SI.GetFractBasesCovered();
	Log("%s", Name.c_str());
	Log(", t=%.1f", double(t1));
	Log(", Cov=%.4f", double(t1));

	static bool HdrDone = false;
	if (!HdrDone)
		{
		Pf(fConsTab, "k");
		Pf(fConsTab, "\tt");
		for (uint i = 0; i < IDS; ++i)
			{
			uint PctId = PctIds[i];
			Pf(fConsTab, "\t%u%%", PctId);
			}
		Pf(fConsTab, "\tall\n");

		Pf(fWinTab, "k");
		Pf(fWinTab, "\tt");
		for (uint i = 0; i < IDS; ++i)
			{
			uint PctId = PctIds[i];
			Pf(fWinTab, "\t%u%%", PctId);
			}
		Pf(fWinTab, "\tall\n");

		HdrDone = true;
		}

	Pf(fConsTab, "%u", k);
	Pf(fConsTab, "\t%u", t);

	Pf(fWinTab, "%u", k);
	Pf(fWinTab, "\t%u", t);

	Log("k=%u", k);
	if (t != UINT_MAX)
		Log(", t=%u", t);
	if (w != UINT_MAX)
		Log(", w=%u", w);
	Log(", t=%.1f", double(t));
	Log(", Size=%.3g", double(Size));
	Log(", Cov=%.4f", Cov);

	for (uint i = 0; i < IDS; ++i)
		{
		uint PctId = PctIds[i];
		double C = GetConservedSyncmerFract(SI, SIs[i]);
		++TestCount;
		double W = GetFractWindowsWithSyncmer(SI, SIs[i], WindowSize);
		Log(", Cons_%u=%.4f", PctId, C);
		Log(", W_%u=%.4f", PctId, W);
		}

	Log("\n");

	SI.Clear();
	for (uint i = 0; i < IDS; ++i)
		SIs[i].Clear();
	}

void cmd_bench1()
	{
	const string Name = opt(bench1);
	fConsTab = CreateStdioFile(opt(constabbedout));
	fWinTab = CreateStdioFile(opt(wintabbedout));
	if (optset_seqlength)
		{
		extern uint BENCHL;
		BENCHL = opt(seqlength);
		}
	if (optset_window)
		WindowSize = opt(window);

	ResetRand(1);

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);
	for (uint i = 0; i < IDS; ++i)
		{
		uint PctId = PctIds[i];
		MutatedSeqs[i] = myalloc(byte, BENCHL);
		MutateSeq(Seq, BENCHL, PctId, MutatedSeqs[i]);
		}

	asserta(optset_k);
	uint k = opt(k);
	if (optset_t)
		Bench1(k, opt(t), UINT_MAX, Name);
	else if (optset_w)
		Bench1(k, UINT_MAX, opt(w), Name);
	else if (optset_tlo)
		{
		asserta(optset_thi);
		for (uint t = opt(tlo); t <= opt(thi); ++t)
			Bench1(k, t, UINT_MAX, Name);
		}
	else if (optset_wlo)
		{
		asserta(optset_whi);
		for (uint w = opt(wlo); w <= opt(whi); ++w)
			Bench1(k, UINT_MAX, w, Name);
		}
	else
		asserta(false);

	CloseStdioFile(fConsTab);
	CloseStdioFile(fWinTab);
	}
