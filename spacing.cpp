#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

static byte *Seq;

uint GetCount(const map<uint, uint> &Map, uint Space)
	{
	map<uint, uint>::const_iterator p = Map.find(Space);
	if (p == Map.end())
		return 0;
	return p->second;
	}

static double GetMeanSpace(const map<uint, uint> &Map)
	{
	uint64 N = 0;
	uint64 Sum = 0;
	for (map<uint, uint>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		uint Space = p->first;
		uint Count = p->second;
		N += Count;
		Sum += uint64(Space)*uint64(Count);
		}
	double Mean = double(Sum)/N;
	return Mean;
	}

uint64 GetTotalCount(const map<uint, uint> &Map)
	{
	uint64 Sum = 0;
	for (map<uint, uint>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		uint Count = p->second;
		Sum += uint64(Count);
		}
	return Sum;
	}

uint64 GetTotalSpace(const map<uint, uint> &Map)
	{
	uint64 Sum = 0;
	for (map<uint, uint>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		uint Space = p->first;
		uint Count = p->second;
		Sum += uint64(Space)*uint64(Count);
		}
	return Sum;
	}

static uint GetMedianSpace(const map<uint, uint> &Map)
	{
	vector<uint> v;
	for (map<uint, uint>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		uint Space = p->first;
		uint Count = p->second;
		for (uint i = 0; i < Count; ++i)
			v.push_back(Space);
		}
	uint Median = v[SIZE(v)/2];
	return Median;
	}

static uint GetModeSpace(const map<uint, uint> &Map)
	{
	uint MaxCount = 0;
	uint MaxSpace = 0;
	for (map<uint, uint>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		uint Count = p->second;
		if (Count > MaxCount)
			{
			uint Space = p->first;
			MaxCount = Count;
			MaxSpace = Space;
			}
		}
	return MaxSpace;
	}

static void Spacing(const vector<string> &Names, uint k, uint t, uint w)
	{
	asserta(t < k);

	const uint NameCount = SIZE(Names);
	asserta(NameCount > 0);

//	uint w = SyncmerIndex::EstimateWindow(k, t);

	Log("k=%u, t=%u, w=%u\n", k, t, w);
	vector<map<uint, uint> > SpaceToCountVec;
	vector<uint64> TotalSpaceVec;
	uint MaxMax = 0;
	for (uint i = 0; i < NameCount; ++i)
		{
		const string &Name = Names[i];
		SyncmerType ST = StrToST(Name);
		SyncmerIndex SI;
		SI.Create(ST, k, t, w, Seq, BENCHL);

		map<uint, uint> SpaceToCount;

		uint Max = SI.GetSpaceToCount(SpaceToCount);
		uint64 TotalSpace = GetTotalSpace(SpaceToCount);
		asserta(TotalSpace > 0);
		double Mean = GetMeanSpace(SpaceToCount);
		uint Median = GetMedianSpace(SpaceToCount);
		uint Mode = GetModeSpace(SpaceToCount);

		uint K = SI.GetKmerCount();
		uint M = SI.GetSyncmerCount();
		double c = double(K)/double(M);

		MaxMax = max(Max, MaxMax);
		SpaceToCountVec.push_back(SpaceToCount);
		TotalSpaceVec.push_back(TotalSpace);

		Log("%s", STToStr(SI.m_ST));
		Log(", max %u", Max);
		Log(", mean %.1f", Mean);
		Log(", median %u", Median);
		Log(", mode %u", Mode);
		Log(", c %.1f", c);
		Log("\n");
		}

	Log("Space");
	for (uint i = 0; i < NameCount; ++i)
		Log("\t%s.f", Names[i].c_str());
	for (uint i = 0; i < NameCount; ++i)
		Log("\t%s.n", Names[i].c_str());
	Log("\n");
	asserta(SIZE(SpaceToCountVec) == NameCount);
	for (unsigned Space = 1; Space <= MaxMax; ++Space)
		{
		Log("%u", Space);
		for (uint i = 0; i < NameCount; ++i)
			{
			uint n = GetCount(SpaceToCountVec[i], Space);
			double Freq = double(n)/double(TotalSpaceVec[i]);
			Log("\t%.4f", Freq);
			}
		for (uint i = 0; i < NameCount; ++i)
			{
			uint n = GetCount(SpaceToCountVec[i], Space);
			Log("\t%u", n);
			}
		Log("\n");
		}
	}

void cmd_spacing()
	{
	vector<string> Names;
	Split(opt(spacing), Names, '+');

	ResetRand(1);
	Seq = myalloc(byte, BENCHL);
	MakeRandSeq(Seq, BENCHL);

	asserta(optset_k);
	asserta(optset_t);
	asserta(optset_w);
	uint t = opt(t);
	uint k = opt(k);
	uint w = opt(w);

	Spacing(Names, k, t, w);
	}
