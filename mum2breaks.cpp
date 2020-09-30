#include "myutils.h"
#include "linereader.h"
#include "seqdb.h"
#include "quarts.h"

/***
   1 [SR] start of the alignment region in the reference sequence
   2 [ER] end of the alignment region in the reference sequence
   3 [SQ] start of the alignment region in the query sequence
   4 [EQ] end of the alignment region in the query sequence
   5 [LEN R] length of the alignment region in the reference sequence
   6 [LEN Q] length of the alignment region in the query sequence
   7 [% IDY] percent identity of the alignment 
  -- [% SIM] percent similarity of the alignment (as determined by the BLOSUM scoring matrix) 
  -- [% STP] percent of stop codons in the alignment 
   8 [LEN R] length of the reference sequence 
   9 [LEN Q] length of the query sequence 
  10 [COV R] percent alignment coverage in the reference sequence 
  11 [COV Q] percent alignment coverage in the query sequence
  -- [FRM] reading frame for the reference and query sequence alignments respectively
  12,13 [TAGS] the reference and query FastA IDs respectively. 
  
  All output coordinates and lengths are relative to the forward 
  strand of the reference DNA sequence.

     0       1          2       3    4       5        6         7       8    9      10            11              12
    SR      ER         SQ      EQ LENR    LENQ      IDY      LENR    LENQ COVR    COVQ        LABELR          LABELQ
  1200    7993    1854386 1847610 6794    6777    87.40   2521574 2222370 0.27    0.30    CP013334.1      CP024698.1
  9209    13177   1846393 1842405 3969    3989    86.73   2521574 2222370 0.16    0.18    CP013334.1      CP024698.1
***/

static FILE *g_fTab;

uint GetOverlap(uint Lo1, uint Hi1, uint Lo2, uint Hi2)
	{
	uint MaxLo = max(Lo1, Lo2);
	uint MinHi = min(Hi1, Hi2);
	if (MaxLo > MinHi)
		return 0;
	return MinHi - MaxLo + 1;
	}

uint InRange(uint Pos, uint Lo, uint Hi, uint Min)
	{
	return Pos >= Lo + Min && Pos + Min <= Hi;
	}

static uint const NBREAKS = 4;
static uint const MIN_ALN = 100;
static uint g_ReadCount;
static uint g_NoAlnCount;
static vector<uint> g_BreakCounts;
static vector<uint> g_CovQuarts;
static vector<uint> g_AlnLengths;
static vector<uint> g_ReadCovPcts;

static void GetMeanStdDev(const vector<uint> &v, uint &Mean, uint &StdDev)
	{
	Quarts Q;
	GetQuarts(v, Q);
	Mean = uint(Q.Avg);
	StdDev = uint(Q.StdDev);
	}

static void GetMedian(const vector<uint> &v, uint &Med)
	{
	Quarts Q;
	GetQuarts(v, Q);
	Med = uint(Q.Med);
	}

static void QuartVecToFreqs(const string &Name,
  const vector<uint> &Quarts, vector<double> &Freqs)
	{
	Freqs.clear();
	Freqs.resize(4, 0.0);
	const uint N = SIZE(Quarts);
	if (N == 0)
		return;

	vector<uint> Counts(4, 0);
	for (uint i = 0; i < N; ++i)
		{
		uint q = Quarts[i];
		asserta(q < 4);
		++Counts[q];
		}
	for (uint q = 0; q < 4; ++q)
		Freqs[q] = double(Counts[q])/double(N);

	Log("%10.10s", Name.c_str());
	for (uint q = 0; q < 4; ++q)
		Log("  %.4f", Freqs[q]);
	Log("\n");

	double SumFreqs = 0.0;
	for (uint q = 0; q < 4; ++q)
		Pf(g_fTab, "\t%.4f", Freqs[q]);
	}

static void Sim(uint ReadLength, uint Iters,
  const vector<uint> &Los, const vector<uint> &His, uint L)
	{
	const unsigned N = SIZE(Los);
	asserta(SIZE(His) == N);
	uint MinFlank = 100;
	uint Margin = 1000;
	asserta(L > 5*ReadLength + 2*Margin);
	uint W = L - ReadLength - 2*Margin;
	for (uint Iter = 0; Iter < Iters; ++Iter)
		{
		uint BreakCount = 0;
		uint ReadLo = Margin + randu32()%W;
		uint ReadHi = ReadLo + ReadLength - 1;
		asserta(ReadHi < L - Margin);
		uint TotalCov = 0;
		for (uint i = 0; i < N; ++i)
			{
			uint Lo = Los[i];
			uint Hi = His[i];
			uint Cov = GetOverlap(ReadLo, ReadHi, Lo, Hi);
			if (Cov >= MIN_ALN)
				{
				g_AlnLengths.push_back(Cov);
				TotalCov += Cov;
				}

			if (InRange(Lo, ReadLo, ReadHi, MinFlank) ||
			    InRange(Hi, ReadLo, ReadHi, MinFlank))
				++BreakCount;
			}

		++g_ReadCount;
		uint CovPct = (TotalCov*100)/ReadLength;
		if (CovPct > 100)
			CovPct = 100;
		g_ReadCovPcts.push_back(CovPct);
		Log("Iter %u breaks %u cov %u%%\n", Iter, BreakCount, CovPct);
		uint CovQuart = (TotalCov*4)/ReadLength;
		if (CovQuart >= 4)
			CovQuart = 3;
		if (TotalCov >= MIN_ALN)
			{
			if (BreakCount > 3)
				BreakCount = 3;
			g_BreakCounts.push_back(BreakCount);
			}
		else
			++g_NoAlnCount;
		g_CovQuarts.push_back(CovQuart);
		}
	}

void cmd_mum2breaks()
	{
	const string &MummerFileName = opt(mum2breaks);
	const string &QueryFileName = opt(input);
	const string &DBFileName = opt(db);
	const string &TabbedFileName = opt(tabbedout);
	uint ReadLength = 20000;
	uint Iters = 100;
	if (optset_readlength)
		ReadLength = opt(readlength);
	if (optset_readlength)
		Iters = opt(iters);

	g_fTab = CreateStdioFile(TabbedFileName);

	SeqDB QDB;
	SeqDB RDB;

	QDB.FromFastx(QueryFileName);
	RDB.FromFasta(DBFileName);

	const uint QSeqCount = QDB.GetSeqCount();
	const uint RSeqCount = RDB.GetSeqCount();

	uint QL = 0;
	uint QSeqIndex = 0;
	for (uint SeqIndex = 0; SeqIndex < QSeqCount; ++SeqIndex)
		{
		uint L = QDB.GetSeqLength(SeqIndex);
		if (L > QL)
			{
			QL = L;
			QSeqIndex = SeqIndex;
			}
		}

	uint RL = 0;
	uint RSeqIndex = 0;
	for (uint SeqIndex = 0; SeqIndex < RSeqCount; ++SeqIndex)
		{
		uint L = RDB.GetSeqLength(SeqIndex);
		if (L > RL)
			{
			RL = L;
			RSeqIndex = SeqIndex;
			}
		}

	const string &QLabel = string(QDB.GetLabel(QSeqIndex));
	const string &RLabel = string(RDB.GetLabel(RSeqIndex));

	ProgressLog("\n");
	ProgressLog("Read length %u\n", ReadLength);
	ProgressLog("Q>%s length %s\n", QLabel.c_str(), MemBytesToStr(QL));
	ProgressLog("R>%s length %s\n", RLabel.c_str(), MemBytesToStr(RL));
	ProgressLog("\n");

	vector<uint> QLos;
	vector<uint> QHis;
	vector<uint> RLos;
	vector<uint> RHis;

	vector<string> Fields;
	t_LineBuff LB;
	LineReader LR;
	LR.Open(MummerFileName);
	while (LR.ReadLine(LB))
		{
		const string Line = string(LB.Data);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) >= 13);

		const string &QLabel2 = Fields[12];
		const string &RLabel2 = Fields[11];

		if (QLabel2 != QLabel || RLabel2 != RLabel)
			continue;

		uint QLo = StrToUint(Fields[2]);
		uint QHi = StrToUint(Fields[3]);

		uint RLo = StrToUint(Fields[0]);
		uint RHi = StrToUint(Fields[1]);

		bool Plus = QLo < QHi;
		if (!Plus)
			swap(QLo, QHi);

		asserta(QLo < QHi);
		asserta(RLo < RHi);

		QLos.push_back(QLo);
		QHis.push_back(QHi);

		RLos.push_back(RLo);
		RHis.push_back(RHi);
		}

	Sim(ReadLength, Iters, QLos, QHis, QL);
	Sim(ReadLength, Iters, RLos, RHis, RL);

	vector<double> BreakFreqs;
	vector<double> CovFreqs;

	uint AlnLength_Med;
	GetMedian(g_AlnLengths, AlnLength_Med);

	uint ReadCovPct_Med;
	GetMedian(g_ReadCovPcts, ReadCovPct_Med);

	uint BreakCount_Med;
	GetMedian(g_BreakCounts, BreakCount_Med);

	Pf(g_fTab, "%s", BaseName(MummerFileName.c_str()));
	Pf(g_fTab, "\t%.4f", double(g_NoAlnCount)/g_ReadCount);

	QuartVecToFreqs("Breaks", g_BreakCounts, BreakFreqs);
	QuartVecToFreqs("Cov", g_CovQuarts, CovFreqs);

	Pf(g_fTab, "\t%u", AlnLength_Med);
	Pf(g_fTab, "\t%u", ReadCovPct_Med);
	Pf(g_fTab, "\t%u", BreakCount_Med);

	Pf(g_fTab, "\n");

	CloseStdioFile(g_fTab);

	Log("f_noaln %.4f\n", double(g_NoAlnCount)/g_ReadCount);
	}
