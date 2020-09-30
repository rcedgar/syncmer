#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

static void Compress1(const byte *Seq, const string &Name,
  uint k, uint t, uint w)
	{
	SyncmerType ST = StrToST(Name);

	SyncmerIndex SI;
	SI.Create(ST, k, t, w, Seq, BENCHL);

	double f = SI.GetFractKmersIndexed();
	double c = 1.0/f;
	ProgressLog("%s", Name.c_str());
	ProgressLog(" k=%u", k);
	if (t != UINT_MAX)
		ProgressLog(" t=%u", t);
	else
		ProgressLog(" w=%u", w);
	ProgressLog(" c=%.1f", c);
	ProgressLog("\n");
	}


void cmd_compress()
	{
	const string Name = opt(compress);

	ResetRand(1);

	byte *Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);

	asserta(optset_k);
	uint k = opt(k);

	if (optset_tlo)
		{
		asserta(optset_thi);
		for (uint t = opt(tlo); t <= opt(thi); ++t)
			Compress1(Seq, Name, k, t, UINT_MAX);
		}
	else if (optset_wlo)
		{
		asserta(optset_whi);
		for (uint w = opt(wlo); w <= opt(whi); ++w)
			Compress1(Seq, Name, k, UINT_MAX, w);
		}
	else
		asserta(false);
	}
