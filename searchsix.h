#ifndef searchhashix_h
#define searchhashix_h

#include "syncmerindex.h"
#include "seqinfo.h"
#include "gobuff.h"
#include <list>

const float MISMATCH_SCORE = -3.0f;
const float GAP_OPEN_SCORE = -5.0f;
const float GAP_EXT_SCORE = -1.0f;
const float XDROP = 16.0f;
const float MIN_GSP_SCORE = 100.0f;
const uint MIN_GSP_LENGTH = 100;

struct SHGSP
	{
public:
	uint32 StartQ;
	uint32 EndQ;
	uint32 LoT;
	uint32 HiT;
	uint32 IdCount;
	float Score;

public:
	SHGSP()
		{
		StartQ = 0;
		EndQ = 0;
		LoT = 0;
		HiT = 0;
		IdCount = 0;
		Score = 0;
		}

public:
	uint GetLength() const { return HiT - LoT + 1; }
	bool GetPlus() const { return StartQ >= EndQ; }
	uint GetLoQ() const { return min(StartQ, EndQ); }
	uint GetHiQ() const { return max(StartQ, EndQ); }
	double GetPctId() const { return IdCount*100.0/GetLength(); }
	};

class SearchHashix
	{
public:
	const SyncmerIndex *m_SI;
	SyncmerIndex m_QSI;
	uint m_WordLength;
	uint m_Stride;
	SeqInfo *m_Query;
	vector<SHGSP *> m_GSPs;
	uint64 m_SumBases;
	double m_SumBasesPctId;
	double m_ANI;
	double m_AFQ;
	double m_AFT;

public:
	static FILE *m_fTab;

public:
	SearchHashix()
		{
		m_SI = 0;
		m_Query = 0;
		m_WordLength = 0;
		m_Stride = 0;
		m_ANI = 0.0;
		m_AFQ = 0.0;
		m_AFT = 0.0;
		m_SumBases = 0;
		m_SumBasesPctId = 0.0;
		}

	void SetHI(const SyncmerIndex *HI)
		{
		m_SI = HI;
		m_WordLength = HI->m_k;
		m_Stride = HI->m_t;
		}

	void Search(SeqInfo *Query);
	SHGSP *Extend(uint32 QPos, uint32 TPos);
	void ClearGSPs();
	void AppendGSP(SHGSP *GSP);
	void LogGSPs() const;
	void GetANI_AFQ(double &ANI, double &AFQ) const;
	};

#endif // searchhashix_h
