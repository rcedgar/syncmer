#ifndef spindex_h
#define spindex_h

#include "objmgr.h"
#include "seqinfo.h"

class Spindex
	{
public:
	uint m_k = 0;
	uint m_Slots = 0;
	byte *m_Counts = 0;
	SeqInfo *m_SI = 0;

public:
	void Init(uint k);
	void Clear()
		{
		m_k = 0;
		m_Slots = 0;

		myfree(m_Counts);
		m_Counts = 0;

		if (m_SI != 0)
			{
			ObjMgr::Down(m_SI);
			m_SI = 0;
			}
		}

	void AddFile_Pass1(const string &FileName);
	void AddSeq_Pass1(const byte *Seq, uint L);

	uint32 GetShiftMask() const
		{
		uint32 ShiftMask = 0;
		for (uint i = 0; i < 2u*m_k; ++i)
			ShiftMask |= (uint32(1) << i);
		return ShiftMask;
		}

	void LogStats() const;
	};

#endif // spindex_h
