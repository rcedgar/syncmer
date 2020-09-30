#ifndef syncmerindex_h
#define syncmerindex_h

#include <map>
#include "kmer.h"
#include "murmur.h"
#include "seqdb.h"

enum SyncmerType
	{
#define t(x)	ST_##x,
#include "stypes.h"
	};

SyncmerType StrToST(const string &s);
const char *STToStr(SyncmerType ST);

class SyncmerIndex
	{
public:
	SyncmerType m_ST = ST_None;
	const char *m_Label = 0;
	const byte *m_Seq = 0;
	uint m_L = 0;

	uint m_k = 0;
	uint m_w = 0;
	uint m_t = 0;

// Per-kmer vectors, length L + k - 1
	vector<uint64> m_Kmers;
	vector<uint64> m_Hashes;
	vector<uint32> m_CoordToSyncmerIx;

// Per-syncmer vector
	vector<uint> m_SyncmerCoords;

// Hash table
	uint m_SlotCount = 0;
	vector<uint32> m_HashTable;

// Temp data for building index
	byte *m_PlusCounts = 0;
	byte *m_MinusCounts = 0;

public:
	void Clear()
		{
		m_ST = ST_None;
		m_Label = 0;
		m_Seq = 0;
		m_L = 0;
		m_k = 0;
		m_w = 0;
		m_t = 0;
		m_Kmers.clear();
		m_Hashes.clear();
		m_CoordToSyncmerIx.clear();
		m_SyncmerCoords.clear();
		m_SlotCount = 0;
		m_HashTable.clear();
		if (m_PlusCounts != 0)
			{
			myfree(m_PlusCounts);
			m_PlusCounts = 0;
			}
		if (m_MinusCounts != 0)
			{
			myfree(m_MinusCounts);
			m_MinusCounts = 0;
			}
		}

public:
	void Create(SyncmerType ST, uint k, uint t, uint w,
	  const byte *Seq, uint L);
	void Create(SyncmerType ST, uint k, uint t, uint w,
	  const char *Label, const byte *Seq, uint L);
	void SetSyncmers();
	void LogSyncmers() const;
	uint GetKmerCount() const { return m_L - m_k + 1; }
	uint GetSyncmerCount() const { return SIZE(m_SyncmerCoords); }
	uint GetWindowCount() const { return GetKmerCount() - m_w + 1; }
	uint GetIndexSize() const { return SIZE(m_SyncmerCoords); }
	uint64 WordToKmer(const byte *Seq) const;
	const byte *KmerToWord(uint64 Kmer, byte *Word) const;
	uint64 KmerToHash(uint64 Kmer) const;
	bool IsSyncmer(uint Coord) const;
	void SetKmers();
	void Validate() const;
	double GetStride() const;
	uint GetBasesCovered() const;
	double GetFractBasesCovered() const;
	double GetFractKmersIndexed() const;
	byte *GetDepthVec() const;
	uint GetSpaceToCount(map<uint, uint> &SpaceToCount) const;
	uint GetDepthToCount(vector<uint> &DepthToCount) const;
	double GetCov1Fract() const;
	uint GetSubmerLength() const;
	void LogRange(uint32 Coord, uint32 n) const;

	bool CalcIsSyncmer(uint32 Coord) const;

	void FromSeqDB(SyncmerType ST, const SeqDB &DB,
	  uint k, uint t, uint w, double Fract);
	// uint64 GetShiftMask() const;
	uint32 GetPos(uint64 Slot) const;
	void SetHashTable(uint32 SlotCount);
	void CountPlus();

#define t(x)	bool Is##x(uint32 Coord) const;
#include "stypes.h"

public:
	static uint EstimateWindow(uint k, uint t);
	static uint EstimateStep(uint k, uint w);
	static bool ValidParams(SyncmerType ST, uint k, uint t, uint w);
	};

double GetConservedSyncmerFract(const SyncmerIndex &SI1,
  const SyncmerIndex &SI2);

double GetFractWindowsWithSyncmer(const SyncmerIndex &SI1,
  const SyncmerIndex &SI2, uint Size);

uint GetSubmerLength(uint k, uint t);

#endif // syncmerindex_h
