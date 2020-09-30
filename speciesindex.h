#ifndef speciesindex_h
#define speciesindex_h

#include "syncmerindex.h"
#include "seqdb.h"
#include <set>

class SpeciesIndex
	{
public:
	SyncmerType m_ST = ST_None;
	uint32 m_k = 0;
	uint32 m_t = 0;
	uint32 m_w = 0;
	uint32 m_Slots = 0;
	uint32 m_TotalSize = 0;
	uint32 m_SpeciesCount = 0;
	vector<string> m_SpeciesNames;
	vector<uint16> m_Sizes;
	vector<uint32> m_Offsets;
	vector<uint16> m_SpVec;
	SyncmerIndex m_SI;
	vector<uint> m_Counts;
	vector<uint> m_Order;

	const char *m_QLabel = 0;
	const byte *m_QSeq = 0;
	uint m_QL = 0;

public:
	void Clear()
		{
		m_ST = ST_None;
		m_k = 0;
		m_t = 0;
		m_w = 0;
		m_Slots = 0;
		m_TotalSize = 0;
		m_SpeciesNames.clear();
		m_Sizes.clear();
		m_Offsets.clear();
		m_SpVec.clear();
		m_SI.Clear();
		}

	void FromSeqDB(SyncmerType ST, uint k, uint t, uint w,
	  uint32 Slots, const SeqDB &DB);
	void ToFile(const string &FileName) const;
	void FromFile(const string &FileName);
	void Validate(const vector<set<uint64> > &HashesVec) const;
	void AppendUniqueSyncmerHashes(const byte *Seq, unsigned L,
	  set<uint64> &Hashes);

	void Search(const char *Label, const byte *Seq, uint L);
	void LogCounts() const;
	};

#endif // speciesindex_h
