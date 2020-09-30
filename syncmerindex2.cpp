#include "myutils.h"
#include "syncmerindex2.h"
#include "alpha.h"
#include "murmur.h"
#include <set>

void SyncmerIndex2::Create(const byte *Seq, uint L)
	{
	ClearVecs();

	m_Seq = Seq;
	m_L = L;

	SetKmers();
	SetSubmers();
	if (opt(validate))
		Validate();
	}

bool SyncmerIndex2::CalcIsSubmer(uint32 Coord) const
	{
	asserta(m_k <= 32);
	if (m_w != UINT_MAX)
		{
		bool IsMin = CalcIsMinimizer(Coord);
		return IsMin;
		}
	else if (m_s != UINT_MAX || m_d != 0)
		{
		bool IsSub = CalcIsSyncmer(Coord);
		return IsSub;
		}
	asserta(false);
	return false;
	}

uint64 SyncmerIndex2::WordToKmer(const byte *Seq) const
	{
	return ::WordToKmer(Seq, m_k);
	}

const byte *SyncmerIndex2::KmerToWord(uint64 Kmer, byte *Word) const
	{
	return ::KmerToWord(Kmer, m_k, Word);
	}

uint64 SyncmerIndex2::KmerToHash(uint64 Kmer) const
	{
	uint64 h = murmur64(Kmer);
	return h;
	}

void SyncmerIndex2::SetKmers()
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

bool SyncmerIndex2::IsSubmer(uint Coord) const
	{
	asserta(Coord < SIZE(m_CoordToSubmerIx));
	uint32 Ix = m_CoordToSubmerIx[Coord];
	return Ix != UINT32_MAX;
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
bool SyncmerIndex2::CalcIsMinimizer(uint32 Coord) const
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

uint SyncmerIndex2::GetUniqueSubmerCount() const
	{
	set<uint64> UniqueSubmers;
	uint SubmerCount = SIZE(m_SubmerCoords);
	for (uint i = 0; i < SubmerCount; ++i)
		{
		uint Pos = m_SubmerCoords[i];
		uint64 Submer = m_Kmers[Pos];
		UniqueSubmers.insert(Submer);
		}
	uint n = SIZE(UniqueSubmers);
	return n;
	}

void SyncmerIndex2::SetSubmers()
	{
	m_CoordToSubmerIx.clear();
	uint K = GetKmerCount();

	m_CoordToSubmerIx.resize(K, UINT32_MAX);
	uint Ix = 0;
	for (uint Coord = 0; Coord < K; ++Coord)
		{
		if (CalcIsSubmer(Coord))
			{
			m_CoordToSubmerIx[Coord] = Ix++;
			m_SubmerCoords.push_back(Coord);
			}
		}
	}

void SyncmerIndex2::LogSubmers() const
	{
	const uint M = SIZE(m_SubmerCoords);
	for (uint i = 0; i < M; ++i)
		{
		uint Pos = m_SubmerCoords[i];
		Log("%*.*s\n", m_k, m_k, m_Seq + Pos);
		}
	}

void SyncmerIndex2::Validate() const
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
			uint SubmerIx = m_CoordToSubmerIx[i];
			asserta(SubmerIx < SIZE(m_SubmerCoords));
			asserta(m_SubmerCoords[SubmerIx] == i);
			}
		else
			asserta(m_CoordToSubmerIx[i] == UINT32_MAX);
		}

	uint N = GetSubmerCount();
	for (uint i = 1; i < N; ++i)
		{
		uint Pos1 = m_SubmerCoords[i-1];
		uint Pos2 = m_SubmerCoords[i];
		asserta(Pos1 < m_L);
		asserta(Pos2 < m_L);
		asserta(Pos1 < Pos2);
		}
	}

void SyncmerIndex2::LogMaxDist() const
	{
	uint Max = 0;
	uint MaxPos1 = 0;
	uint MaxPos2 = 0;
	uint N = GetSubmerCount();
	for (uint i = 1; i < N; ++i)
		{
		uint Pos1 = m_SubmerCoords[i-1];
		uint Pos2 = m_SubmerCoords[i];
		asserta(Pos2 > Pos1);
		uint Space = Pos2 - Pos1;
		asserta(Space > 0);
		if (Space > Max)
			{
			Max = Space;
			MaxPos1 = Pos1;
			MaxPos2 = Pos2;
			}
		}

	for (uint Pos = MaxPos1; Pos <= MaxPos2; ++Pos)
		{
		uint64 Kmer = m_Kmers[Pos];
		char YN = ' ';
		if (IsSubmer(Pos))
			YN = '<';
		char YN2 = ' ';
		if (CalcIsSubmer(Pos))
			YN2 = '*';
		uint MinPos = GetMinSubkmerPos_Hash(Kmer, m_k, m_s);

//		Log("%5u", Pos);
		Log("%*.*s", m_k, m_k, m_Seq + Pos);
	
		byte Word[64];
		KmerToWord(Kmer, Word);
		Log("  %c", YN);
		Log("  %c", YN2);
		Log("  minpos=%u", MinPos);
		Log(" %*.*s", m_s, m_s, Word + MinPos);

		for (uint Pos = 0; Pos + m_s <= m_k; ++Pos)
			{
			uint64 SubKmer = GetSubkmer(Kmer, m_k, m_s, Pos);
			uint64 HashSubKmer = murmur64(SubKmer);
			Log(" %*.*s=%x", m_s, m_s, Word + Pos, HashSubKmer);
			}

		Log("\n");
		}
	}

void SyncmerIndex2::LogRange(uint32 Coord, uint32 n) const
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
		if (IsSubmer(i))
			YN = '<';
		Log("%5u", i);
		Log("  %*.*s", m_k, m_k, m_Seq + i);
		Log("  %s", KmerToWord(Kmer, Word));
		Log("  %c", YN);
		Log("\n");
		}
	}

double SyncmerIndex2::GetStride() const
	{
	uint K = GetKmerCount();
	asserta(SIZE(m_CoordToSubmerIx) == K);
	uint n = 0;
	for (uint i = 0; i < K; ++i)
		{
		if (m_CoordToSubmerIx[i] != UINT32_MAX)
			++n;
		}
	asserta(n > 0);
	double Stride = double(K)/double(n);
	return Stride;
	}

double SyncmerIndex2::GetFractKmersIndexed() const
	{
	uint n = GetIndexSize();
	uint K = GetKmerCount();
	asserta(n <= K);
	double Fract = double(n)/double(K);
	return Fract;
	}

double SyncmerIndex2::GetMeanDepth(double *ptrF0) const
	{
	vector<uint> DepthToCount;
	uint MaxDepth = GetDepthToCount(DepthToCount);
	asserta(SIZE(DepthToCount) == MaxDepth + 1);
	uint SumDepth = 0;
	for (uint d = 0; d <= MaxDepth; ++d)
		SumDepth += d*DepthToCount[d];
	uint n0 = DepthToCount[0];
	*ptrF0 = double(n0)/m_L;
	double MeanDepth = double(SumDepth)/m_L;
	return MeanDepth;
	}

double SyncmerIndex2::GetCompressionFactor() const
	{
	double f = GetFractKmersIndexed();
	asserta(f > 0);
	double c = 1.0/f;
	return c;
	}

double SyncmerIndex2::GetFractBasesCovered() const
	{
	uint n = GetBasesCovered();
	asserta(n <= m_L);
	double Fract = double(n)/double(m_L);
	return Fract;
	}

uint SyncmerIndex2::GetBasesCovered() const
	{
	uint K = GetSubmerCount();
	asserta(SIZE(m_SubmerCoords) == K);

	vector<bool> Cov(m_L, false);

	for (uint i = 0; i < K; ++i)
		{
		uint Pos = m_SubmerCoords[i];
		for (uint j = 0; j < m_k; ++j)
			Cov[Pos+j] = true;
		}

	uint n = 0;
	for (uint i = 0; i < m_L; ++i)
		if (Cov[i])
			++n;
	asserta(n <= m_L);

	return n;
	}

double GetConservedSubmerFract(const SyncmerIndex2 &SI1,
  const SyncmerIndex2 &SI2)
	{
	uint N = 0;
	uint n = 0;
	uint K = SI1.GetKmerCount();
	for (uint i = 0; i < K; ++i)
		{
		if (SI1.IsSubmer(i))
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

// d=1 means consecutive.
// SpaceToCount[0] always zero.
uint SyncmerIndex2::GetSpaceToCount(vector<uint> &SpaceToCount) const
	{
	SpaceToCount.resize(1);
	SpaceToCount[0] = 0;
	uint Max = 0;
	uint N = GetSubmerCount();
	for (uint i = 1; i < N; ++i)
		{
		uint Pos1 = m_SubmerCoords[i-1];
		uint Pos2 = m_SubmerCoords[i];
		asserta(Pos2 > Pos1);
		uint Space = Pos2 - Pos1;
		asserta(Space > 0);
		if (Space > Max)
			{
			Max = Space;
			SpaceToCount.resize(Max+1);
			}
		++(SpaceToCount[Space]);
		}
	asserta(SpaceToCount[0] == 0);
	return Max;
	}

uint *SyncmerIndex2::GetDepthVec() const
	{
	uint *D = myalloc(uint, m_L);
	memset(D, 0, m_L*sizeof(uint));

	uint N = GetSubmerCount();
	for (uint i = 0; i < N; ++i)
		{
		uint Pos = m_SubmerCoords[i];
		for (uint j = 0; j < m_k; ++j)
			{
			if (D[Pos+j] < UINT_MAX)
				++(D[Pos+j]);
			}
		}
	return D;
	}

uint SyncmerIndex2::GetDepthToCount(vector<uint> &DepthToCount) const
	{
	DepthToCount.clear();
	DepthToCount.resize(1);
	uint *D = GetDepthVec();
	uint Max = 0;
	for (uint i = 0; i < m_L; ++i)
		{
		byte d = D[i];
		if (d > Max)
			{
			Max = d;
			DepthToCount.resize(Max+1);
			}
		++(DepthToCount[d]);
		}
	myfree(D);
	return Max;
	}

double SyncmerIndex2::GetCov1Fract() const
	{
	uint N = 0;
	uint n1 = 0;
	uint *D = GetDepthVec();
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

bool SyncmerIndex2::CalcIsSyncmer(uint32 Coord) const
	{
	assert(Coord < SIZE(m_Hashes));
	if (m_s == UINT_MAX)
		{
		asserta(m_d != 0);
		uint64 h = m_Hashes[Coord];
		return h%m_d == 0;
		}

	asserta(m_s >= 2);
	if (!m_Open)
		asserta(m_s + 2 <= m_k);
	uint64 Kmer = m_Kmers[Coord];
	uint MinPos = GetMinSubkmerPos_Hash(Kmer, m_k, m_s);

	bool Yes = false;
	if (m_Open)
		Yes = (MinPos == m_k - m_s);
	else
		Yes = (MinPos == m_ds || MinPos == m_k - m_s);
	if (Yes && m_d != 0)
		{
		asserta(m_d > 1);
		uint64 Hash = m_Hashes[Coord];
		bool Yesd = (Hash%m_d == 0);
		return Yesd;
		}
	return Yes;
	}
