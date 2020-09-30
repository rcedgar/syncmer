#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

uint BENCHL = 1024*1024;

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

void InitBench2(uint WinSiz)
	{
	WindowSize = WinSiz;
	fConsTab = CreateStdioFile(opt(tabbedout));
	if (optset_seqlength)
		{
		extern uint BENCHL;
		BENCHL = opt(seqlength);
		}

	ResetRand(1);

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);
	for (uint i = 0; i < IDS; ++i)
		{
		uint PctId = PctIds[i];
		MutatedSeqs[i] = myalloc(byte, BENCHL);
		MutateSeq(Seq, BENCHL, PctId, MutatedSeqs[i]);
		}
	}

void Bench2(uint k, uint t, uint w, const string &Name1, const string &Name2)
	{
	//if (t + 1 >= k)
	//	return;

	SyncmerType ST1 = StrToST(Name1);
	SyncmerType ST2 = StrToST(Name2);

	if (!SyncmerIndex::ValidParams(ST1, k, t, w))
		return;
	if (!SyncmerIndex::ValidParams(ST2, k, t, w))
		return;

	SyncmerIndex SI1;
	SyncmerIndex SI2;
	SI1.Create(ST1, k, t, w, Seq, BENCHL);
	SI2.Create(ST2, k, t, w, Seq, BENCHL);

	vector<SyncmerIndex> SI1s(IDS);
	vector<SyncmerIndex> SI2s(IDS);
	for (uint i = 0; i < IDS; ++i)
		{
		const byte *MutatedSeq = MutatedSeqs[i];
		SI1s[i].Create(ST1, k, t, w, MutatedSeq, BENCHL);
		SI2s[i].Create(ST2, k, t, w, MutatedSeq, BENCHL);
		}

	uint Size1 = SI1.GetIndexSize();
	uint Size2 = SI2.GetIndexSize();
	int dSize = int(Size1) - int(Size2);

	double t1 = 1.0/SI1.GetFractKmersIndexed();
	double t2 = 1.0/SI2.GetFractKmersIndexed();

	double Cov1 = SI1.GetFractBasesCovered();
	double Cov2 = SI2.GetFractBasesCovered();
	double dCov = Cov1 - Cov2;

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
	Log(", t=%u", t);
	Log(", w=%u", w);

	Log(", t1=%.1f", double(t1));
	Log(", t2=%.1f", double(t2));

	Log(", Size1=%.3g", double(Size1));
	Log(", Size2=%.3g", double(Size2));
	Log(", dSize=%+.1f%%", double(dSize*200)/double(Size1+Size2));

	Log(", Cov1=%.4f", Cov1);
	Log(", Cov2=%.4f", Cov2);
	Log(", dCov=%.4f", dCov);

	CovDeltas.push_back(dCov);
	if (Cov1 > Cov2)
		++CovBetterCount1;
	else if (Cov2 > Cov1)
		++CovBetterCount2;
	else
		{
		asserta(Cov1 == Cov2);
		++CovTieCount;
		}

	if (Size1 < Size2)
		++SizeBetterCount1;
	else if (Size2 < Size1)
		++SizeBetterCount2;
	else
		{
		asserta(Size1 == Size2);
		++SizeTieCount;
		}

	double SumDeltaCons = 0.0;
	double SumDeltaWin = 0.0;
	for (uint i = 0; i < IDS; ++i)
		{
		uint PctId = PctIds[i];
		double C1 = GetConservedSyncmerFract(SI1, SI1s[i]);
		double C2 = GetConservedSyncmerFract(SI2, SI2s[i]);
		++TestCount;
		if (C1 > C2)
			++ConsBetterCount1;
		else if (C2 > C1)
			++ConsBetterCount2;
		else
			{
			asserta(C1 == C2);
			++ConsTieCount;
			}
		double DeltaCons = C1 - C2;
		SumDeltaCons += DeltaCons;
		Pf(fConsTab, "\t%.4f", DeltaCons);

		ConsDeltas.push_back(DeltaCons);
		ConsDeltasPctIds[i].push_back(DeltaCons);
		Log(", Cons1_%u=%.4f", PctId, C1);
		Log(", Cons2_%u=%.4f", PctId, C2);
		Log(", dCons_%u=%.4f", PctId, DeltaCons);


		double W1 = GetFractWindowsWithSyncmer(SI1, SI1s[i], WindowSize);
		double W2 = GetFractWindowsWithSyncmer(SI2, SI2s[i], WindowSize);
		double DeltaW = W1 - W2;
		SumDeltaWin += DeltaW;
		Pf(fWinTab, "\t%.4f", DeltaW);

		Log(", W1_%u=%.4f", PctId, W1);
		Log(", W2_%u=%.4f", PctId, W2);
		Log(", dW_%u=%.4f", PctId, DeltaW);

		if (W1 > W2)
			++WinBetterCount1;
		else if (W2 > W1)
			++WinBetterCount2;
		else
			{
			asserta(W1 == W2);
			++WinTieCount;
			}
		}

	double MeanDeltaCons = SumDeltaCons/IDS;
	Pf(fConsTab, "\t%.4f", MeanDeltaCons);
	Pf(fConsTab, "\n");

	double MeanDeltaWin = SumDeltaWin/IDS;
	Pf(fWinTab, "\t%.4f", MeanDeltaWin);
	Pf(fWinTab, "\n");

	Log("\n");

	SI1.Clear();
	SI2.Clear();
	for (uint i = 0; i < IDS; ++i)
		{
		SI1s[i].Clear();
		SI2s[i].Clear();
		}
	}

void cmd_bench()
	{
	const string Names = opt(bench);
	vector<string> Fields;
	Split(Names, Fields, '+');
	asserta(SIZE(Fields) == 2);
	const string &Name1 = Fields[0];
	const string &Name2 = Fields[1];
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

	if (optset_k)
		{
		uint k = opt(k);
		uint t = opt(t);
		uint w = opt(w);
		Progress("name1=%s, name2=%s, k=%u, t=%u, w=%u\n",
		  Name1.c_str(), Name2.c_str(), k, t, w);
		Bench2(k, t, w, Name1, Name2);
		}
	else if (opt(tile))
		{
		for (uint k = 8; k <= 16; ++k)
			{
			uint t = k - 1;
			Progress("tile k=%u\n", k);
			uint w = (optset_w ? opt(w) : SyncmerIndex::EstimateWindow(k, t));
			Bench2(k, t, w, Name1, Name2);
			}
		}
	else
		{
		for (uint k = 8; k <= 16; ++k)
			{
			Progress("k=%u\n", k);
			for (uint t = 2; t <= 10; ++t)
				{
				uint w = (optset_w ? opt(w) : SyncmerIndex::EstimateWindow(k, t));
				Bench2(k, t, w, Name1, Name2);
				}
			}
		}

	ProgressLog("Cons: %s better %u (%.1f%%)",
	  Name1.c_str(), ConsBetterCount1, GetPct(ConsBetterCount1, TestCount));

	ProgressLog(", %s better %u (%.1f%%)",
	  Name2.c_str(), ConsBetterCount2, GetPct(ConsBetterCount2, TestCount));

	if (ConsTieCount > 0)
		ProgressLog(", ties %u (%.1f%%)",
		  ConsTieCount, GetPct(ConsTieCount, TestCount));
	else
		ProgressLog(", no ties");
	ProgressLog("\n");

	ProgressLog("Cov: %s better %u (%.1f%%)",
	  Name1.c_str(), CovBetterCount1, GetPct(CovBetterCount1, TestCount));

	ProgressLog(", %s better %u (%.1f%%)",
	  Name2.c_str(), CovBetterCount2, GetPct(CovBetterCount2, TestCount));

	if (CovTieCount > 0)
		ProgressLog(", ties %u (%.1f%%)",
		  CovTieCount, GetPct(CovTieCount, TestCount));
	else
		ProgressLog(", no ties");
	ProgressLog("\n");

	ProgressLog("Size: %s better %u (%.1f%%)",
	  Name1.c_str(), SizeBetterCount1, GetPct(SizeBetterCount1, TestCount));

	ProgressLog(", %s better %u (%.1f%%)",
	  Name2.c_str(), SizeBetterCount2, GetPct(SizeBetterCount2, TestCount));

	if (SizeTieCount > 0)
		ProgressLog(", ties %u (%.1f%%)",
		  SizeTieCount, GetPct(SizeTieCount, TestCount));
	else
		ProgressLog(", no ties");
	ProgressLog("\n");

	ProgressLog("Win: %s better %u (%.1f%%)",
	  Name1.c_str(), WinBetterCount1, GetPct(WinBetterCount1, TestCount));

	ProgressLog(", %s better %u (%.1f%%)",
	  Name2.c_str(), WinBetterCount2, GetPct(WinBetterCount2, TestCount));

	if (WinTieCount > 0)
		ProgressLog(", ties %u (%.1f%%)",
		  WinTieCount, GetPct(WinTieCount, TestCount));
	else
		ProgressLog(", no ties");
	ProgressLog("\n");

	sort(CovDeltas.begin(), CovDeltas.end());
	double MedianCovDelta = CovDeltas[SIZE(CovDeltas)/2];
	ProgressLog("Cov difference : median %+.2f%%\n", MedianCovDelta*100.0);

	sort(ConsDeltas.begin(), ConsDeltas.end());
	double MedianConsDelta = ConsDeltas[TestCount/2];
	ProgressLog("Cons difference: median %+.2f%%\n", MedianConsDelta*100.0);
	ProgressLog("Cons by pctid: ");
	for (uint i = 0; i < IDS; ++i)
		{
		uint PctId = PctIds[i];
		vector<double> &ConsDeltas = ConsDeltasPctIds[i];
		uint n = SIZE(ConsDeltas);
		sort(ConsDeltas.begin(), ConsDeltas.end());
		MedianConsDelta = ConsDeltas[n/2];
		if (i > 0)
			ProgressLog(", ");
		ProgressLog("%u%%=%+.2f%%", PctId, MedianConsDelta*100.0);
		}
	ProgressLog("\n");

	CloseStdioFile(fConsTab);
	}
