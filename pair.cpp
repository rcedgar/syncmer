#include "myutils.h"
#include "seqdb.h"
#include "syncmerindex.h"

static void MakeStr(const SyncmerIndex &SI1, const SyncmerIndex &SI2,
  char c, string &s)
	{
	s.clear();
	uint K = SI1.GetKmerCount();
	for (uint i = 0; i < K; ++i)
		{
		if (SI1.IsSyncmer(i))
			{
			if (SI2.IsSyncmer(i))
				s += toupper(c);
			else
				s += tolower(c);
			}
		else
			s += ' ';
		}
	}

double Correl(const vector<double> &x, const vector<double> &y)
	{
	const uint N = SIZE(x);
	asserta(N > 0 && SIZE(y) == N);
	double Sumx = 0.0;
	double Sumy = 0.0;
	double Sumx2  = 0.0;
	double Sumy2  = 0.0;
	double Sumxy  = 0.0;
	for (uint i = 0; i < N; ++i)
		{
		double X = x[i];
		double Y = y[i];

		Sumx += X;
		Sumy += Y;
		Sumx2 += X*X;
		Sumy2 += Y*Y;
		Sumxy += X*Y;
		}

	double Top = Sumxy - (Sumx*Sumy)/N;
	double BottomL = Sumx2 - (Sumx*Sumx)/N;
	double BottomR = Sumy2 - (Sumy*Sumy)/N;
	double r = Top/sqrt(BottomL*BottomR);
	return r;
	}

static void DoPair(SyncmerType ST1, SyncmerType ST2,
  const char *LabelA, const char *LabelB,
  const byte *SeqA, const byte *SeqB, uint L,
  uint k, uint t, uint w, double &F1, double &F2, double &PctId)
	{
	F1 = -1.0;
	F2 = -1.0;
	PctId = -1.0;

	uint Ids = 0;
	for (uint i = 0; i < L; ++i)
		if (toupper(SeqA[i]) == toupper(SeqB[i]))
			++Ids;
	PctId = (100.0*Ids)/L;

	SyncmerIndex SIA1;
	SyncmerIndex SIA2;
	SyncmerIndex SIB1;
	SyncmerIndex SIB2;

	SIA1.Create(ST1, k, t, w, SeqA, L);
	SIA2.Create(ST2, k, t, w, SeqA, L);

	SIB1.Create(ST1, k, t, w, SeqB, L);
	SIB2.Create(ST2, k, t, w, SeqB, L);

#if 0
	{
	string StrA1;
	string StrA2;
	string StrB1;
	string StrB2;

	MakeStr(SIA1, SIB1, 'x', StrA1);
	MakeStr(SIB1, SIA1, 'x', StrB1);

	MakeStr(SIA2, SIB2, 'y', StrA2);
	MakeStr(SIB2, SIA2, 'y', StrB2);

	Log("  A1: %s\n", StrA1.c_str());
	Log("  B1: %s\n", StrB1.c_str());
	Log("SeqA: %*.*s\n", L, L, SeqA);
	Log("SeqB: %*.*s\n", L, L, SeqB);
	Log("  A2: %s\n", StrA2.c_str());
	Log("  B2: %s\n", StrB2.c_str());
	}
#endif

	F1 = GetConservedSyncmerFract(SIA1, SIB1);
	F2 = GetConservedSyncmerFract(SIA2, SIB2);

	//Log("%.4f", F1);
	//Log("\t%.4f", F2);
	//Log("\t%.1f", PctId);
	//Log("\t%s", LabelA);
	//Log("\t%s", LabelB);
	//Log("\n");
	}

void cmd_pair()
	{
	const string Names = opt(pair);
	vector<string> Fields;
	Split(Names, Fields, '+');
	asserta(SIZE(Fields) == 2);
	const string &Name1 = Fields[0];
	const string &Name2 = Fields[1];
	asserta(optset_k && optset_t);
	uint k = opt(k);
	uint t = opt(t);
	uint w = (optset_w ? opt(w) : SyncmerIndex::EstimateWindow(k, t));
	Log("k=%u, t=%u, w=%u\n", k, t, w);
	ResetRand(1);

	SyncmerType ST1 = StrToST(Name1);
	SyncmerType ST2 = StrToST(Name2);

	const string &FastaFileName = opt(input);

	SeqDB Pairs;
	Pairs.FromFasta(FastaFileName);
	const uint SeqCount = Pairs.GetSeqCount();

	const uint L = Pairs.GetSeqLength(0);
	for (uint i = 0; i < SeqCount; ++i)
		{
		const char *LabelA = Pairs.GetLabel(i);
		const byte *SeqA = Pairs.GetSeq(i);
		asserta(Pairs.GetSeqLength(i) == L);
		for (uint j = i; j < SeqCount; ++j)
			{
			const byte *SeqB = Pairs.GetSeq(j);
			const char *LabelB = Pairs.GetLabel(j);
			asserta(Pairs.GetSeqLength(j) == L);

			double F1, F2, PctId;
			DoPair(ST1, ST2, LabelA, LabelB, SeqA, SeqB, L, k, t, w,
			  F1, F2, PctId);
			}
		}
	}

void cmd_pairs()
	{
	const string Names = opt(pairs);
	vector<string> Fields;
	Split(Names, Fields, '+');
	asserta(SIZE(Fields) == 2);
	const string &Name1 = Fields[0];
	const string &Name2 = Fields[1];
	ResetRand(1);
	FILE *fTab = CreateStdioFile(opt(tabbedout));

	SyncmerType ST1 = StrToST(Name1);
	SyncmerType ST2 = StrToST(Name2);

	const string &FastaFileName = opt(input);

	SeqDB Pairs;
	Pairs.FromFasta(FastaFileName);
	const uint SeqCount = Pairs.GetSeqCount();

	for (uint k = 8; k <= 20; ++k)
		{
		uint t = k;
		uint w = SyncmerIndex::EstimateWindow(k, t);
		ProgressLog("k=%u, w=%u", k, w);
		vector<double> PctIds;
		vector<double> F1s;
		vector<double> F2s;
		const uint L = Pairs.GetSeqLength(0);
		for (uint i = 0; i < SeqCount; ++i)
			{
			const char *LabelA = Pairs.GetLabel(i);
			const byte *SeqA = Pairs.GetSeq(i);
			asserta(Pairs.GetSeqLength(i) == L);
			for (uint j = i; j < SeqCount; ++j)
				{
				const byte *SeqB = Pairs.GetSeq(j);
				const char *LabelB = Pairs.GetLabel(j);
				asserta(Pairs.GetSeqLength(j) == L);

				double F1, F2, PctId;
				DoPair(ST1, ST2, 
				  LabelA, LabelB, SeqA, SeqB, L, 
				  k, t, w,
				  F1, F2, PctId);
				PctIds.push_back(PctId);
				F1s.push_back(F1);
				F2s.push_back(F2);
				Pf(fTab, "@P");
				Pf(fTab, "\tk=%u", k);
				Pf(fTab, "\tt=%u", t);
				Pf(fTab, "\tw=%u", w);
				Pf(fTab, "\tA=%s", LabelA);
				Pf(fTab, "\tB=%s", LabelB);
				Pf(fTab, "\tPctId=%.1f", PctId);
				Pf(fTab, "\tF1=%.4f", F1);
				Pf(fTab, "\tF2=%.4f", F2);
				Pf(fTab, "\n");
				}
			}

		double r1 = Correl(F1s, PctIds);
		double r2 = Correl(F2s, PctIds);

		Pf(fTab, "@R");
		Pf(fTab, "\tk=%u", k);
		Pf(fTab, "\tt=%u", t);
		Pf(fTab, "\tw=%u", w);
		Pf(fTab, "\tr1=%.4f", r1);
		Pf(fTab, "\tr2=%.4f", r2);
		Pf(fTab, "\n");

		ProgressLog(", r1=%.2f, r2=%.2f\n", r1, r2);
		}
	CloseStdioFile(fTab);
	}
