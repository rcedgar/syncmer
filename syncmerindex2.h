#ifndef syncmerindex_h
#define syncmerindex_h

#include <map>
#include "kmer.h"
#include "murmur.h"
#include "seqdb.h"

class SyncmerIndex2
	{
public:
	const char *m_Label = 0;
	const byte *m_Seq = 0;
	uint m_L = 0;
	uint m_k = UINT_MAX;		// k-mer length
	uint m_w = UINT_MAX;		// window size for minimizer
	uint m_s = UINT_MAX;		// substring size for syncmer
	uint m_ds = 0;				// offset from start for substring
	uint m_d = 0;				// downsample
	bool m_Open = false;

// Per-kmer vectors, length L + k - 1
	vector<uint64> m_Kmers;
	vector<uint64> m_Hashes;
	vector<uint32> m_CoordToSubmerIx;

// Per-syncmer vector
	vector<uint> m_SubmerCoords;

public:
	void Create(const byte *Seq, uint L);
	void ClearVecs()
		{
		m_Kmers.clear();
		m_Hashes.clear();
		m_CoordToSubmerIx.clear();
		m_SubmerCoords.clear();
		}

	bool CalcIsMinimizer(uint32 Coord) const;
	bool CalcIsSyncmer(uint32 Coord) const;
	void CopyParams(const SyncmerIndex2 &SI)
		{
		m_k = SI.m_k;
		m_w = SI.m_w;
		m_s = SI.m_s;
		m_ds = SI.m_ds;
		m_d = SI.m_d;
		m_Open = SI.m_Open;
		}

	const char *GetParamStr(string &Str, char cSep = ',') const
		{
		Ps(Str, "k=%u", m_k);
		if (m_s != UINT_MAX)
			Psa(Str, "%cs=%u", cSep, m_s);
		if (m_w != UINT_MAX)
			Psa(Str, "%cw=%u", cSep, m_w);
		if (m_ds > 0)
			Psa(Str, "%cds=%u", cSep, m_ds);
		if (m_d > 0)
			Psa(Str, "%cd=%u", cSep, m_d);
		if (m_Open)
			Psa(Str, "%copen=yes", cSep);
		return Str.c_str();
		}

	void SetSubmers();
	void LogSubmers() const;
	uint GetKmerCount() const { return m_L - m_k + 1; }
	uint GetSubmerCount() const { return SIZE(m_SubmerCoords); }
	uint GetWindowCount() const { return GetKmerCount() - m_w + 1; }
	uint GetIndexSize() const { return SIZE(m_SubmerCoords); }
	uint GetUniqueSubmerCount() const;
	uint64 WordToKmer(const byte *Seq) const;
	const byte *KmerToWord(uint64 Kmer, byte *Word) const;
	uint64 KmerToHash(uint64 Kmer) const;
	bool IsSubmer(uint Coord) const;
	void SetKmers();
	void Validate() const;
	double GetStride() const;
	uint GetBasesCovered() const;
	double GetFractBasesCovered() const;
	double GetFractKmersIndexed() const;
	double GetCompressionFactor() const;
	double GetMeanDepth(double *ptrF0) const;
	uint *GetDepthVec() const;
	uint GetSpaceToCount(vector<uint> &SpaceToCount) const;
	uint GetDepthToCount(vector<uint> &DepthToCount) const;
	double GetCov1Fract() const;
	void LogRange(uint32 Coord, uint32 n) const;
	void LogMaxDist() const;

	bool CalcIsSubmer(uint32 Coord) const;
	};

//double GetConservedSyncmerFract(const SyncmerIndex2 &SI1,
//  const SyncmerIndex2 &SI2);

#endif // syncmerindex_h
