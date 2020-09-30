#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

static byte *Seq;
static bool Trace = false;

static void FindEquiv(SyncmerType ST, uint k, uint t)
	{
	SyncmerIndex SIS;
	SIS.Create(ST, k, t, 0, Seq, BENCHL);
	uint K = SIS.GetKmerCount();
	uint SSize = SIS.GetIndexSize();
	double SStride = K/double(SSize);

	uint EstimatedWindow = SyncmerIndex::EstimateWindow(k, t);
	uint wHi = UINT_MAX;
	uint wLo = UINT_MAX;
	double MStrideLo = 0.0;
	double MStrideHi = 999.9;
	for (int iw = int(EstimatedWindow); iw > 1; --iw)
		{
		uint w = uint(iw);
		SyncmerIndex SIM;
		SIM.Create(ST_Minimizer1, k, 0, w, Seq, BENCHL);
		uint KM = SIM.GetKmerCount();
		asserta(KM == K);
		uint MSize = SIM.GetIndexSize();
		double MStride = K/double(MSize);
		if (Trace)
			Log("(-) w %u, stride %.1f\n", w, MStride);

		if (MStride <= t && MStride > MStrideLo)
			{
			MStrideLo = MStride;
			wLo = w;
			}
		if (MStride >= t && MStride < MStrideHi)
			{
			MStrideHi = MStride;
			wHi = w;
			}

		if (MStride <= t)
			break;
		}

	for (uint w = EstimatedWindow; w < EstimatedWindow + 32; ++w)
		{
		SyncmerIndex SIM;
		SIM.Create(ST_Minimizer1, k, 0, w, Seq, BENCHL);
		uint KM = SIM.GetKmerCount();
		asserta(KM == K);
		uint MSize = SIM.GetIndexSize();
		double MStride = K/double(MSize);
		if (Trace)
			Log("(+) w %u, stride %.1f\n", w, MStride);

		if (MStride <= t && MStride > MStrideLo)
			{
			MStrideLo = MStride;
			wLo = w;
			}
		if (MStride >= t && MStride < MStrideHi)
			{
			MStrideHi = MStride;
			wHi = w;
			}

		if (MStride >= t)
			break;
		}

	double DiffHi = MStrideHi - double(t);
	double DiffLo = MStrideLo - double(t);
	uint wBest;
	double DiffBest;
	if (fabs(DiffHi) <= fabs(DiffLo))
		{
		wBest = wHi;
		DiffBest = DiffHi;
		}
	else
		{
		wBest = wLo;
		DiffBest = DiffLo;
		}

	ProgressLog("k=%u", k);
	ProgressLog(", t=%u", t);
	ProgressLog(", SStride=%.2f", SStride);
	ProgressLog(", wLo=%u", wLo);
	ProgressLog(", MStrideLo=%.2f", MStrideLo);
	ProgressLog(", wHi=%u", wHi);
	ProgressLog(", MStrideHi=%.2f", MStrideHi);
	ProgressLog(", MStrideHi=%.2f", MStrideHi);
	ProgressLog(", wBest=%u", wBest);
	ProgressLog(", dBest=%.2g", DiffBest);
	ProgressLog(", wEst=%u", EstimatedWindow);
	ProgressLog(", dw=%+d", int(wBest - EstimatedWindow));
	ProgressLog("\n");
	}

void cmd_equiv()
	{
	const string Name = opt(equiv);

	ResetRand(1);
	SyncmerType ST = StrToST(Name);

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);
	for (uint k = 8; k <= 16; ++k)
		for (uint t = 2; t <= 10; ++t)
			FindEquiv(ST, k, t);
	}

void cmd_equiv1()
	{
	const string Name = opt(equiv1);

	ResetRand(1);
	SyncmerType ST = StrToST(Name);

	asserta(optset_k && optset_t);
	const uint k = opt(k);
	const uint t = opt(t);

	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);
	Trace = true;
	FindEquiv(ST, k, t);
	}
