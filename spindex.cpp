#include "myutils.h"
#include "spindex.h"
#include "alpha.h"
#include "fastaseqsource.h"

void Spindex::Init(uint k)
	{
	Clear();

	m_SI = ObjMgr::GetSeqInfo();

	m_k = k;
	m_Slots = myipow(4, k);
	m_Counts = myalloc(byte, m_Slots);
	memset_zero(m_Counts, m_Slots);
	}

void Spindex::AddFile_Pass1(const string &FileName)
	{
	FASTASeqSource SS;
	SS.Open(FileName);
	while (SS.GetNext(m_SI))
		{
		const byte *Seq = m_SI->m_Seq;
		const uint L = m_SI->m_L;
		AddSeq_Pass1(Seq, L);
		m_SI->RevCompInPlace();
		AddSeq_Pass1(Seq, L);
		}
	}

void Spindex::AddSeq_Pass1(const byte *Seq, uint L)
	{
	const uint32 ShiftMask = GetShiftMask();

	uint Count = 0;
	uint32 Word = 0;
	byte K = 0;
	for (uint SeqPos = 0; SeqPos < L - m_k + 1; ++SeqPos)
		{
		byte c = Seq[SeqPos];
		byte Letter = g_CharToLetterNucleo[c];
		if (Letter == INVALID_LETTER)
			{
			K = 0;
			Word = 0;
			continue;
			}
		if (K < m_k)
			++K;
		Word = ((Word << uint32(2)) | Letter) & ShiftMask;
		if (K == m_k)
			{
			assert(Word < m_Slots);
			byte Count = m_Counts[Word];
			if (Count < 0xff)
				m_Counts[Word] = Count + 1;
			}
		}
	}

void Spindex::LogStats() const
	{
	Log("k=%u\n", m_k);
	Log("%u slots\n", m_Slots);

	vector<uint> Counts(256);
	for (uint Word = 0; Word < m_Slots; ++Word)
		{
		byte Count = m_Counts[Word];
		++(Counts[Count]);
		}

	for (uint i = 0; i < 256; ++i)
		Log("%u\t%u\n", i, m_Counts[i]);
	}
