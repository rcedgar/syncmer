#include "myutils.h"
#include "speciesindex.h"
#include "sort.h"
#include <map>
#include <set>

void SpeciesIndex::Search(const char *Label, const byte *Seq, uint L)
	{
	//TraceCounts(Label, Seq, L);
	//return;
	m_QLabel = Label;
	m_QSeq = Seq;
	m_QL = L;

	m_Counts.clear();
	m_Counts.resize(m_SpeciesCount, 0);

	set<uint64> Hashes;
	AppendUniqueSyncmerHashes(Seq, L, Hashes);

	for (set<uint64>::const_iterator p = Hashes.begin();
	  p != Hashes.end(); ++p)
		{
		uint64 Hash = *p;
		uint32 Slot = uint32(Hash%m_Slots);
		uint Size = m_Offsets[Slot+1] - m_Offsets[Slot];
		uint Offset = m_Offsets[Slot];
		asserta(Offset + Size <= m_TotalSize);
		for (uint j = 0; j < Size; ++j)
			{
			uint16 SpIndex = m_SpVec[Offset+j];
			asserta(SpIndex < m_SpeciesCount);
			++m_Counts[SpIndex];
			}
		}
	m_Order.resize(m_SpeciesCount);
	QuickSortOrderDesc(m_Counts.data(), m_SpeciesCount, m_Order.data());
	}

void SpeciesIndex::LogCounts() const
	{
	asserta(SIZE(m_Counts) == m_SpeciesCount);
	asserta(SIZE(m_Order) == m_SpeciesCount);

	Log("%s[%u]", m_QLabel, m_QL);
	for (uint k = 0; k < min(m_SpeciesCount, 8u); ++k)
		{
		uint SpIndex = m_Order[k];
		asserta(SpIndex < m_SpeciesCount);
		uint Count = m_Counts[SpIndex];
		const string &SpeciesName = m_SpeciesNames[SpIndex];
		if (Count == 0)
			break;
		Log("  %s(%u)", SpeciesName.c_str(), Count);
		}
	Log("\n");
	}

void cmd_species_search()
	{
	const string &InputFileName = opt(species_search);
	const string &DBFileName = opt(db);
	const string &OutputFileName = opt(output);
	FILE *fOut = CreateStdioFile(OutputFileName);

	SpeciesIndex SPI;
	SPI.FromFile(DBFileName);

	SeqDB Input;
	Input.FromFasta(InputFileName);
	const uint SeqCount = Input.GetSeqCount();

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Searching");
		const char *Label = Input.GetLabel(SeqIndex);
		uint L = Input.GetSeqLength(SeqIndex);
		const byte *Seq = Input.GetSeq(SeqIndex);
		SPI.Search(Label, Seq, L);
		SPI.LogCounts();
		}
	}
