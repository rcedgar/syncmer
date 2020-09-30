#include "myutils.h"
#include "syncmerindex.h"
#include "alpha.h"
#include "murmur.h"

SyncmerType StrToST(const string &s)
	{
#define t(x)	if (s == #x) return ST_##x;
#include "stypes.h"
	Die("Unknown type %s", s.c_str());
	return ST_None;
	}

const char *STToStr(SyncmerType ST)
	{
	switch (ST)
		{
#define t(x)	case ST_##x: return #x;
#include "stypes.h"
		}
	asserta(false);
	return "ST_??";
	}

void SyncmerIndex::Create(SyncmerType ST, uint k, uint t, uint w,
  const byte *Seq, uint L)
	{
	Create(ST, k, t, w, "_nolabel_", Seq, L);
	}

bool SyncmerIndex::ValidParams(SyncmerType ST, uint k, uint t, uint w)
	{
	if (k > 32)
		return false;

	if (ST == ST_Syncmer2)
		{
	//	asserta(k > t + 1);
		if (k <= t + 1)
			return false;
		}

	if (ST != ST_Syncmer1)
		{
		if (k <= t)
			return false;
		}

	if (ST != ST_Syncmer5)
		{
		//int m = (int) k - (int) t - 1;
		//if (m < 1)
		//	return false;
		}

	return true;
	}

void SyncmerIndex::Create(SyncmerType ST, uint k, uint t, uint w,
  const char *Label, const byte *Seq, uint L)
	{
	Clear();

	if (ST != ST_Syncmer1 && ST != ST_Minimizer1)
		asserta(t < k);

	m_ST = ST;
	m_k = k;
	m_t = t;
	m_w = w;
	m_Label = Label;
	m_Seq = Seq;
	m_L = L;

	SetKmers();
	SetSyncmers();
	if (opt(validate))
		Validate();
	}

bool SyncmerIndex::CalcIsSyncmer(uint32 Coord) const
	{
	switch (m_ST)
		{
#define t(x)	case ST_##x: return Is##x(Coord);
#include "stypes.h"
		}
	asserta(false);
	return false;
	}

uint64 SyncmerIndex::WordToKmer(const byte *Seq) const
	{
	return ::WordToKmer(Seq, m_k);
	}

const byte *SyncmerIndex::KmerToWord(uint64 Kmer, byte *Word) const
	{
	return ::KmerToWord(Kmer, m_k, Word);
	}

uint64 SyncmerIndex::KmerToHash(uint64 Kmer) const
	{
	uint64 h = murmur64(Kmer);
	return h;
	}

void SyncmerIndex::SetKmers()
	{
	m_Kmers.clear();
	m_Hashes.clear();

	asserta(m_L > 0 && m_k > 0);
	asserta(m_L > m_k);
	uint K = m_L - m_k + 1;

	m_Kmers.reserve(K);
	m_Hashes.reserve(K);

	for (uint i = 0; i < K; ++i)
		{
		uint64 Kmer = WordToKmer(m_Seq + i);
		uint64 Hash = KmerToHash(Kmer);

		m_Kmers.push_back(Kmer);
		m_Hashes.push_back(Hash);
		}
	}

bool SyncmerIndex::IsNone(uint Coord) const
	{
	asserta(false);
	return false;
	}

bool SyncmerIndex::IsStride(uint Coord) const
	{
	return Coord%m_t == 0;
	}

bool SyncmerIndex::IsSyncmer(uint Coord) const
	{
	asserta(Coord < SIZE(m_CoordToSyncmerIx));
	uint32 Ix = m_CoordToSyncmerIx[Coord];
	return Ix != UINT32_MAX;
	}

bool SyncmerIndex::IsSyncmer1(uint32 Coord) const
	{
	assert(Coord < SIZE(m_Hashes));
	uint64 h = m_Hashes[Coord];
	bool Yes = (h%m_t == 0);
	return Yes;
	}

uint SyncmerIndex::GetSubmerLength() const
	{
	asserta(m_k > m_t + 1);
	uint Length = m_k - m_t;
	return Length;
	}

bool SyncmerIndex::IsSyncmer2(uint32 Coord) const
	{
	assert(Coord < SIZE(m_Hashes));
	const uint m = GetSubmerLength();
	uint64 Kmer = m_Kmers[Coord];
	uint MinPos = GetMinSubkmerPos(Kmer, m_k, m);
	bool Yes = (MinPos%m_t == 0);
	return Yes;
	}

bool SyncmerIndex::IsSyncmer3(uint32 Coord) const
	{
	assert(Coord < SIZE(m_Hashes));
	uint64 Kmer = m_Kmers[Coord];
	uint MinPos = GetMinSubkmerPos_Rotate(Kmer, m_k, 2);
	bool Yes = (MinPos%m_t == 0);
	return Yes;
	}

bool SyncmerIndex::IsSyncmer4(uint32 Coord) const
	{
	assert(Coord < SIZE(m_Hashes));
	uint64 Kmer = m_Kmers[Coord];
	uint MinPos = GetMinSubkmerPos_Rotate(Kmer, m_k, 2);
	bool Yes = (MinPos%m_t == 1);
	return Yes;
	}

bool SyncmerIndex::IsSyncmer5(uint32 Coord) const
	{
	assert(Coord < SIZE(m_Hashes));
	uint64 Kmer = m_Kmers[Coord];
	uint MinPos = GetMinSubkmerPos_Hash(Kmer, m_k, m_t);
	bool Yes = (MinPos == 0 || MinPos == m_k - m_t - 1);
	return Yes;
	}

/***
k=6, w=3
Three consecutive k-mers, indow has 8 bases:
  123456--
  -123456-
  --123456

A kmer appears in 3 8-base windows:
   WW123456
    W123456W
	 123456WW
***/
bool SyncmerIndex::IsMinimizer1(uint32 Coord) const
	{
	uint K = GetKmerCount();
	asserta(Coord < K);

	bool Yes = false;
	for (int w = 0; w < int(m_w); ++w)
		{
		int WindowStart = int(Coord) - w;
		if (WindowStart < 0)
			continue;
		uint64 MinHash = m_Hashes[WindowStart];
		for (uint i = 1; i < m_w; ++i)
			{
			uint Pos = uint(WindowStart) + i;
			if (Pos < K)
				{
				uint64 h = m_Hashes[Pos];
				if (h <= MinHash)
					MinHash = h;
				}
			}
		if (m_Hashes[Coord] == MinHash)
			{
			Yes = true;
			break;
			}
		}
	return Yes;
	}

bool SyncmerIndex::IsMinimizer2(uint32 Coord) const
	{
	uint K = GetKmerCount();
	asserta(Coord < K);

	bool Yes = false;
	for (int w = 0; w < int(m_w); ++w)
		{
		int WindowStart = int(Coord) - w;
		if (WindowStart < 0)
			continue;
		uint64 MinKmer = m_Kmers[WindowStart];
		for (uint i = 1; i < m_w; ++i)
			{
			uint64 h = m_Kmers[WindowStart + int(i)];
			if (h <= MinKmer)
				MinKmer = h;
			}
		if (m_Kmers[Coord] == MinKmer)
			{
			Yes = true;
			break;
			}
		}
	return Yes;
	}

void SyncmerIndex::SetSyncmers()
	{
	m_CoordToSyncmerIx.clear();
	uint K = GetKmerCount();

	m_CoordToSyncmerIx.resize(K, UINT32_MAX);
	uint Ix = 0;
	for (uint Coord = 0; Coord < K; ++Coord)
		{
		if (CalcIsSyncmer(Coord))
			{
			m_CoordToSyncmerIx[Coord] = Ix++;
			m_SyncmerCoords.push_back(Coord);
			}
		}
	}

void SyncmerIndex::LogSyncmers() const
	{
	const uint M = SIZE(m_SyncmerCoords);
	for (uint i = 0; i < M; ++i)
		{
		uint Pos = m_SyncmerCoords[i];
		Log("%*.*s\n", m_k, m_k, m_Seq + Pos);
		}
	}

void SyncmerIndex::Validate() const
	{
	uint K = GetKmerCount();
	for (uint i = 0; i < K; ++i)
		{
		uint64 Kmer = WordToKmer(m_Seq + i);
		asserta(Kmer == m_Kmers[i]);

		uint64 Hash = KmerToHash(Kmer);
		asserta(Hash == m_Hashes[i]);
		asserta(Hash == KmerToHash(Kmer));

		if (CalcIsSyncmer(i))
			{
			uint SyncmerIx = m_CoordToSyncmerIx[i];
			asserta(SyncmerIx < SIZE(m_SyncmerCoords));
			asserta(m_SyncmerCoords[SyncmerIx] == i);
			}
		else
			asserta(m_CoordToSyncmerIx[i] == UINT32_MAX);
		}

	uint N = GetSyncmerCount();
	for (uint i = 1; i < N; ++i)
		{
		uint Pos1 = m_SyncmerCoords[i-1];
		uint Pos2 = m_SyncmerCoords[i];
		asserta(Pos1 < m_L);
		asserta(Pos2 < m_L);
		asserta(Pos1 < Pos2);
		}
	}

void SyncmerIndex::LogRange(uint32 Coord, uint32 n) const
	{
	byte Word[64];

	uint Hi = Coord + n;
	uint K = GetKmerCount();
	if (Hi > K)
		Hi = K;
	for (uint i = Coord; i < Hi; ++i)
		{
		uint64 Kmer = m_Kmers[i];
		char YN = ' ';
		if (IsSyncmer(i))
			YN = '<';
		Log("%5u", i);
		Log("  %*.*s", m_k, m_k, m_Seq + i);
		Log("  %s", KmerToWord(Kmer, Word));
		Log("  %c", YN);
		Log("\n");
		}
	}

double SyncmerIndex::GetStride() const
	{
	uint K = GetKmerCount();
	asserta(SIZE(m_CoordToSyncmerIx) == K);
	uint n = 0;
	for (uint i = 0; i < K; ++i)
		{
		if (m_CoordToSyncmerIx[i] != UINT32_MAX)
			++n;
		}
	asserta(n > 0);
	double Stride = double(K)/double(n);
	return Stride;
	}

double SyncmerIndex::GetFractKmersIndexed() const
	{
	uint n = GetIndexSize();
	uint K = GetKmerCount();
	asserta(n <= K);
	double Fract = double(n)/double(K);
	return Fract;
	}

double SyncmerIndex::GetFractBasesCovered() const
	{
	uint n = GetBasesCovered();
	asserta(n <= m_L);
	double Fract = double(n)/double(m_L);
	return Fract;
	}

uint SyncmerIndex::GetBasesCovered() const
	{
	uint K = GetKmerCount();
	asserta(SIZE(m_CoordToSyncmerIx) == K);

	vector<bool> Cov(m_L, false);

	for (uint i = 0; i < K; ++i)
		{
		if (m_CoordToSyncmerIx[i] != UINT32_MAX)
			{
			for (uint j = 0; j < m_k; ++j)
				Cov[i+j] = true;
			}
		}

	uint n = 0;
	for (uint i = 0; i < m_L; ++i)
		if (Cov[i])
			++n;
	asserta(n <= m_L);

	return n;
	}

uint SyncmerIndex::EstimateWindow(uint k, uint t)
	{
	uint w = 2*t - 1;
	return w;
	}

uint SyncmerIndex::EstimateStep(uint k, uint w)
	{
	uint t = (w + 1)/2;
	return t;
	}

double GetConservedSyncmerFract(const SyncmerIndex &SI1,
  const SyncmerIndex &SI2)
	{
	asserta(SI1.m_L == SI2.m_L);
	asserta(SI1.m_ST == SI2.m_ST);
	asserta(SI1.m_k == SI2.m_k);
	asserta(SI1.m_t == SI2.m_t || SI1.m_w == SI2.m_w);

	uint N = 0;
	uint n = 0;
	uint K = SI1.GetKmerCount();
	for (uint i = 0; i < K; ++i)
		{
		if (SI1.IsSyncmer(i))
			{
			++N;
			if (SI1.m_Kmers[i] == SI2.m_Kmers[i])
				++n;
			}
		}
	if (N == 0)
		{
		asserta(n == 0);
		return 0.0;
		}
	double Fract = double(n)/N;
	return Fract;
	}

uint SyncmerIndex::GetSpaceToCount(map<uint, uint> &SpaceToCount) const
	{
	SpaceToCount.clear();

	uint N = GetSyncmerCount();
	for (uint i = 1; i < N; ++i)
		{
		uint Pos1 = m_SyncmerCoords[i-1];
		uint Pos2 = m_SyncmerCoords[i];
		asserta(Pos2 > Pos1);
		uint Space = Pos2 - Pos1 - 1;
		if (SpaceToCount.find(Space) == SpaceToCount.end())
			SpaceToCount[Space] = 1;
		else
			++(SpaceToCount[Space]);
		}
	map<uint, uint>::const_iterator p = SpaceToCount.end();
	--p;
	uint Max = p->first;
	return Max;
	}

byte *SyncmerIndex::GetDepthVec() const
	{
	byte *D = myalloc(byte, m_L);
	memset(D, 0, m_L);

	uint N = GetSyncmerCount();
	for (uint i = 0; i < N; ++i)
		{
		uint Pos = m_SyncmerCoords[i];
		for (uint j = 0; j < m_k; ++j)
			{
			if (D[Pos+j] < 255)
				++(D[Pos+j]);
			}
		}
	return D;
	}

uint SyncmerIndex::GetDepthToCount(vector<uint> &DepthToCount) const
	{
	DepthToCount.clear();
	DepthToCount.resize(256);
	byte *D = GetDepthVec();
	uint Max = 0;
	for (uint i = 0; i < m_L; ++i)
		{
		byte d = D[i];
		if (d > Max)
			Max = d;
		++(DepthToCount[d]);
		}
	myfree(D);
	return Max;
	}

double SyncmerIndex::GetCov1Fract() const
	{
	uint N = 0;
	uint n1 = 0;
	byte *D = GetDepthVec();
	for (uint i = 0; i < m_L; ++i)
		{
		byte d = D[i];
		if (d > 0)
			{
			++N;
			if (d == 1)
				++n1;
			}
		}
	if (N == 0)
		return 0.0;
	return double(n1)/double(N);
	}

double GetFractWindowsWithSyncmer(const SyncmerIndex &SI1,
  const SyncmerIndex &SI2, uint Size)
	{
	const uint L = SI1.m_L;
	const uint k = SI1.m_k;

	asserta(SI2.m_L == L);
	asserta(SI2.m_k == k);

	asserta(SI1.m_ST == SI2.m_ST);
	asserta(SI1.m_t == SI2.m_t || SI1.m_w == SI2.m_w);

	uint n = 0;
	const uint K = SI1.GetKmerCount();
	const uint WindowCount = L - Size + 1;
	const uint KmersPerWindow = Size - k + 1;
	for (uint WindowStart = 0; WindowStart < WindowCount; ++WindowStart)
		{
		for (uint i = 0; i < KmersPerWindow; ++i)
			{
			uint Pos = WindowStart + i;
			if (SI1.IsSyncmer(Pos) &&
			  (SI1.m_Kmers[Pos] == SI2.m_Kmers[Pos]))
				{
				++n;
				break;
				}
			}
		}

	double Fract = double(n)/WindowCount;
	return Fract;
	}

void SyncmerIndex::FromSeqDB(SyncmerType ST, const SeqDB &DB,
  uint k, uint t, uint w, double Fract)
	{
	uint SeqCount = DB.GetSeqCount();
	uint TotalL = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		uint L = DB.GetSeqLength(i);
		TotalL += 2*L + 3; // nuls
		}
	byte *DBSeq = myalloc(byte, TotalL);
	byte *DBSeqPtr = DBSeq;
	*DBSeqPtr++ = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const byte *Seq = DB.GetSeq(i);
		uint L = DB.GetSeqLength(i);
		memcpy(DBSeqPtr, Seq, L);
		DBSeqPtr += L;
		*DBSeqPtr++ = 0;
		byte *RCSeq = myalloc(byte, L);
		void RevCompSeq(const byte *Seq, unsigned L, byte *RCSeq);
		RevCompSeq(Seq, L, RCSeq);
		memcpy(DBSeqPtr, RCSeq, L);
		DBSeqPtr += L;
		*DBSeqPtr++ = 0;
		}

	Create(ST, k, t, w, DBSeq, TotalL);
	uint SlotCount = uint(TotalL/Fract);
	SetHashTable(SlotCount);
	}

//uint64 SyncmerIndex::GetShiftMask() const
//	{
//	uint64 ShiftMask = 0;
//	for (uint i = 0; i < 2u*m_k; ++i)
//		ShiftMask |= (uint64(1) << i);
//	return ShiftMask;
//	}
