#include "myutils.h"
#include "syncmerindex2.h"
#include "randseq.h"

static byte *Seq;
static byte *MutatedSeq;
static FILE *f;
static FILE *f2;
static FILE *fPctId;
static FILE *fSp;

static void GetConsFracts(const SyncmerIndex2 &SI1,
  const SyncmerIndex2 &SI2, double *ptrSubmerFract, double *ptrBaseFract)
	{
	asserta(SI1.m_L == SI2.m_L);
	asserta(SI1.m_k == SI2.m_k);
	asserta(SI1.m_s == SI2.m_s || SI1.m_w == SI2.m_w);

	const uint L = SI1.m_L;
	const uint k = SI1.m_k;
	vector<bool> CvgVect(L, false);

	uint N = 0;
	uint n = 0;
	uint K = SI1.GetKmerCount();
	for (uint Pos = 0; Pos < K; ++Pos)
		{
		if (SI1.CalcIsSubmer(Pos))
			{
			++N;
			if (SI1.m_Kmers[Pos] != SI2.m_Kmers[Pos])
				continue;
			if (SI2.CalcIsSubmer(Pos))
				{
				for (uint j = 0; j < k; ++j)
					CvgVect[Pos+j] = true;
				++n;
				}
			}
		}
	if (N == 0)
		{
		asserta(n == 0);
		*ptrSubmerFract = 0;
		*ptrBaseFract = 0;
		}

	*ptrSubmerFract = double(n)/N;

	uint m = 0;
	for (uint j = 0; j < L; ++j)
		if (CvgVect[j])
			++m;
	*ptrBaseFract = double(m)/L;
	}

static void GetConsFracts_PctId(const SyncmerIndex2 &SI, uint PctId,
  double *ptrSubmerFract, double *ptrBaseFract)
	{
	SyncmerIndex2 SI2;
	SI2.CopyParams(SI);

	MutateSeq(Seq, BENCHL, PctId, MutatedSeq);
	SI2.Create(MutatedSeq, BENCHL);

	GetConsFracts(SI, SI2, ptrSubmerFract, ptrBaseFract);
	}

static uint PctIdLo;
static uint PctIdHi;
static uint PctIdInc;

static void Bench(SyncmerIndex2 &SI)
	{
	string ParamStr;
	SI.GetParamStr(ParamStr);

	SI.Create(Seq, BENCHL);

	const uint k = SI.m_k;
	const uint s = SI.m_s;
	double c = SI.GetCompressionFactor();
	double f_zero;
	double Depth = SI.GetMeanDepth(&f_zero);
	double f1 = 1.0 - f_zero;

	Pf(fSp, "%s\n", ParamStr.c_str());
	Pf(fPctId, "%s\n", ParamStr.c_str());
	Progress("%s c=%.2f depth=%.2g f1=%.4f\n",
	  SI.GetParamStr(ParamStr, '\t'), c, Depth, f1);

	static bool HdrDone = false;
	if (!HdrDone)
		{
		Pf(f2, "k");
		Pf(f2, "\tw");
		Pf(f2, "\ts");
		Pf(f2, "\td");
		Pf(f2, "\tc");
		Pf(f2, "\tdepth");
		Pf(f2, "\tf1");
		for (uint PctId = PctIdLo; PctId <= PctIdHi; PctId += PctIdInc)
			{
			Pf(f2, "\tFS%u", PctId);
			Pf(f2, "\tFL%u", PctId);
			}
		Pf(f2, "\n");
		HdrDone = true;
		}

	Pf(f2, "%u", SI.m_k);
	Pf(f2, "\t%u", SI.m_w == UINT_MAX ? 0 : SI.m_w);
	Pf(f2, "\t%u", SI.m_s == UINT_MAX ? 0 : SI.m_s);
	Pf(f2, "\t%u", SI.m_d);
	Pf(f2, "\t%.3g", c);
	Pf(f2, "\t%.3g", Depth);
	Pf(f2, "\t%.4f", f1);

	Pf(f, "%s", ParamStr.c_str());
	Pf(f, "\tc=%.2f", c);
	Pf(f, "\tdepth=%.2f", Depth);
	Pf(f, "\tf1=%.4f", f1);

	uint wp = SI.m_w;
	bool IsSyn = (SI.m_w == UINT_MAX);
	if (IsSyn)
		{
		wp = k - s;
		Pf(f, "\twp=%u", wp);
		}

	Pf(fPctId, "PctId\tSmers\tBases\n");
	for (uint PctId = PctIdLo; PctId <= PctIdHi; PctId += PctIdInc)
		{
		double SF;
		double BF;
		GetConsFracts_PctId(SI, PctId, &SF, &BF);

		Pf(f2, "\t%.3g", SF);
		Pf(f2, "\t%.3g", BF);

		Pf(f, "\tFS%u=%.3g", PctId, SF);
		Pf(f, "\tFL%u=%.3g", PctId, BF);

		Pf(fPctId, "%u", PctId);
		Pf(fPctId, "\t%.3g", SF);
		Pf(fPctId, "\t%.3g", BF);
		Pf(fPctId, "\n");
		}
	Pf(f, "\n");
	Pf(f2, "\n");

	vector<uint> SpaceToCount;
	uint Maxd = SI.GetSpaceToCount(SpaceToCount);
	asserta(SIZE(SpaceToCount) == Maxd + 1);
	asserta(Maxd > 0);

	Pf(fSp, "maxd=%u\n", Maxd);
	uint MaxCount = 0;
	uint TotalCount = 0;
	for (uint d = 1; d <= Maxd; ++d)
		{
		uint n = SpaceToCount[d];
		MaxCount = max(n, MaxCount);
		TotalCount += n;
		}
	for (uint d = 1; d <= Maxd; ++d)
		{
		uint n = SpaceToCount[d];
		double freq = double(n)/TotalCount;
		double frel = double(n)/MaxCount;
		Pf(fSp, "%u\t%u\t%.3g\t%.3g\n", d, n, freq, frel);
		}

//	SI.LogMaxDist();
	}

static void BigTable()
	{
	PctIdLo = 80;
	PctIdHi = 90;
	PctIdInc = 10;

	SyncmerIndex2 SI;
	for (uint k = 5; k <= 32; ++k)
		{
		for (uint s = 2; s < k - 2; ++s)
			{
			uint w = k - s;

			SI.m_d = 0;

		// Minimizer
			SI.m_k = k;
			SI.m_w = w;
			SI.m_s = UINT_MAX;
			Bench(SI);

		// Syncmer
			SI.m_w = UINT_MAX;
			SI.m_s = s;
			Bench(SI);

		// Downsampled
			for (SI.m_d = 2; SI.m_d <= 4; ++SI.m_d)
				Bench(SI);
			}
		}
	}

void cmd_newbench()
	{
	const string TabbedFileName = opt(newbench);
	f = CreateStdioFile(TabbedFileName);
	f2 = CreateStdioFile(opt(tabbedout2));
	if (f != 0)
		setbuf(f, 0);
	if (f2 != 0)
		setbuf(f2, 0);
	fPctId = CreateStdioFile(opt(pctidout));
	fSp = CreateStdioFile(opt(spacingout));
	if (optset_seqlength)
		BENCHL = opt(seqlength);
	ResetRand(1);

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);
	MutatedSeq = myalloc(byte, BENCHL);

	if (opt(bigtable))
		{
		BigTable();
		return;
		}

	SyncmerIndex2 SI;
	SI.m_s = UINT_MAX;
	SI.m_w = UINT_MAX;
	SI.m_ds = 0;
	SI.m_d = 0;

	SI.m_ds = 0;
	SI.m_Open = opt(open);
	if (optset_ds)
		SI.m_ds = opt(ds);
	if (optset_d)
		SI.m_d = opt(d);

	PctIdLo = 80;
	PctIdHi = 80;
	PctIdInc = 2;
	if (optset_pctidinc)
		PctIdInc = opt(pctidinc);
	if (optset_pctid)
		{
		PctIdLo = opt(pctid);
		PctIdHi = opt(pctid);
		}
	else if (optset_pctidlo && optset_pctidhi)
		{
		PctIdLo = opt(pctidlo);
		PctIdHi = opt(pctidhi);
		}

	uint klo = opt(klo);
	uint khi = opt(khi);
	if (optset_k)
		{
		klo = opt(k);
		khi = opt(k);
		}

	for (SI.m_k = klo; SI.m_k <= khi; SI.m_k++)
		{
		if (optset_s && !optset_slo && !optset_shi)
			{
			SI.m_s = opt(s);
			Bench(SI);
			}
		else if (!optset_s && optset_slo && optset_shi)
			{
			for (uint s = opt(slo); s <= opt(shi); ++s)
				{
				SI.m_s = s;
				Bench(SI);
				}
			}
		else if (optset_w && !optset_wlo && !optset_whi)
			{
			SI.m_w = opt(w);
			Bench(SI);
			}
		else if (!optset_w && optset_wlo && optset_whi)
			{
			for (uint w = opt(wlo); w <= opt(whi); ++w)
				{
				SI.m_w = w;
				Bench(SI);
				}
			}
		else if (optset_d)
			{
			SI.m_d = opt(d);
			Bench(SI);
			}
		}

	CloseStdioFile(f);
	CloseStdioFile(f2);
	CloseStdioFile(fPctId);
	CloseStdioFile(fSp);
	}
