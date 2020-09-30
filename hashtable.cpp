#include "myutils.h"
#include "syncmerindex.h"
#include "alpha.h"

void SyncmerIndex::CountPlus()
	{
	asserta(m_PlusCounts == 0);
	m_PlusCounts = myalloc(byte, m_SlotCount);
	zero(m_PlusCounts, m_SlotCount);

	const uint K = GetKmerCount();
	for (uint SeqPos = 0; SeqPos < K; ++SeqPos)
		{
		if (IsSyncmer(SeqPos))
			{
			uint64 Kmer = m_Kmers[SeqPos];
			uint64 Hash = KmerToHash(Kmer);
			uint64 Slot = Hash%m_SlotCount;
			if (m_PlusCounts[Slot] < 255)
				++m_PlusCounts[Slot];
			}
		}
	}

void SyncmerIndex::SetHashTable(uint SlotCount)
	{
	m_SlotCount = SlotCount;
	m_HashTable.clear();
	m_HashTable.resize(m_SlotCount, UINT32_MAX);

	CountPlus();
	const uint K = GetKmerCount();
	for (uint SeqPos = 0; SeqPos < K; ++SeqPos)
		{
		if (IsSyncmer(SeqPos))
			{
			uint64 Kmer = m_Kmers[SeqPos];
			uint64 Hash = KmerToHash(Kmer);
			uint64 Slot = Hash%m_SlotCount;
			if (m_PlusCounts[Slot] == 1)
				{
				asserta(Slot < SIZE(m_HashTable));
				m_HashTable[Slot] = SeqPos;
				}
			}
		}
	}

uint32 SyncmerIndex::GetPos(uint64 Slot) const
	{
	asserta(Slot < SIZE(m_HashTable));
	uint32 Pos = m_HashTable[Slot];
	return Pos;
	}
