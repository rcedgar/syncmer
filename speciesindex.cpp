#include "myutils.h"
#include "speciesindex.h"
#include <map>

static const char *GetSpeciesNameFromLabel(const string &Label, string &SpeciesName)
	{
	vector<string> Fields;
	Split(Label, Fields, '.');
	SpeciesName = Fields[0];
	return SpeciesName.c_str();
	}

void SpeciesIndex::AppendUniqueSyncmerHashes(const byte *Seq, unsigned L,
  set<uint64> &Hashes)
	{
	m_SI.Create(m_ST, m_k, m_t, m_w, Seq, L);
	const uint SyncmerCount = SIZE(m_SI.m_SyncmerCoords);

	set<uint64> UniqueHashes;
	for (uint i = 0; i < SyncmerCount; ++i)
		{
		uint Coord = m_SI.m_SyncmerCoords[i];
		uint64 Hash = m_SI.m_Hashes[Coord];
		Hashes.insert(Hash);
		}
	}

void SpeciesIndex::Validate(const vector<set<uint64> > &HashesVec) const
	{
	const uint m_SpeciesCount = SIZE(m_SpeciesNames);
	asserta(SIZE(HashesVec) == m_SpeciesCount);
	asserta(SIZE(m_Sizes) == m_Slots);
	asserta(SIZE(m_Offsets) == m_Slots + 1);
	asserta(SIZE(m_SpVec) == m_TotalSize);
	for (uint SpeciesIndex = 0; SpeciesIndex < m_SpeciesCount; ++SpeciesIndex)
		{
		const set<uint64> &Hashes = HashesVec[SpeciesIndex];
		bool Found = false;
		for (set<uint64>::const_iterator p = Hashes.begin();
		  p != Hashes.end(); ++p)
			{
			uint64 Hash = *p;
			uint32 Slot = uint32(Hash%m_Slots);
			uint Size = m_Sizes[Slot];
			uint Offset = m_Offsets[Slot];
			asserta(Size + Offset <= m_TotalSize);
			for (uint j = 0; j < Size; ++j)
				{
				uint16 SpeciesIndex2 = m_SpVec[Offset+j];
				if (SpeciesIndex2 == SpeciesIndex)
					{
					Found = true;
					break;
					}
				}
			}
		if (!Found)
			Die("Not found");
		}
	}

void GetSpeciesToSeqIndexes(const SeqDB &DB,
  vector<string> &SpeciesNames,
  vector<vector<uint> > &SeqIndexesVec)
	{
	SpeciesNames.clear();
	SeqIndexesVec.clear();

	const uint SeqCount = DB.GetSeqCount();
	map<string, uint> SpeciesNameToIndex;
	vector<uint> EmptyVec;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string Label = string(DB.GetLabel(SeqIndex));
		string SpeciesName;
		GetSpeciesNameFromLabel(Label, SpeciesName);
		map<string, uint>::const_iterator p = 
		  SpeciesNameToIndex.find(SpeciesName);

		uint SpeciesIndex = UINT_MAX;
		if (p == SpeciesNameToIndex.end())
			{
			SpeciesIndex = SIZE(SpeciesNames);
			SeqIndexesVec.push_back(EmptyVec);
			SpeciesNames.push_back(SpeciesName);
			}
		else
			SpeciesIndex = p->second;
		SeqIndexesVec[SpeciesIndex].push_back(SeqIndex);
		}
	}

void SpeciesIndex::FromSeqDB(SyncmerType ST, uint k, uint t, uint w,
  uint32 Slots, const SeqDB &DB)
	{
	m_ST = ST;
	m_k = k;
	m_t = t;
	m_w = w;
	m_Slots = Slots;

	const uint SeqCount = DB.GetSeqCount();

	vector<vector<uint> > SeqIndexesVec;
	GetSpeciesToSeqIndexes(DB, m_SpeciesNames, SeqIndexesVec);
	m_SpeciesCount = SIZE(m_SpeciesNames);

	vector<set<uint64> > HashesVec(m_SpeciesCount);
	for (uint SpeciesIndex = 0; SpeciesIndex < m_SpeciesCount; ++SpeciesIndex)
		{
		const string &SpeciesName = m_SpeciesNames[SpeciesIndex];
		ProgressStep(SpeciesIndex, m_SpeciesCount, SpeciesName.c_str());
		const vector<uint> &SeqIndexes = SeqIndexesVec[SpeciesIndex];

		for (uint i = 0; i < SIZE(SeqIndexes); ++i)
			{
			uint SeqIndex = SeqIndexes[i];
			const byte *Seq = DB.GetSeq(SeqIndex);
			const unsigned L = DB.GetSeqLength(SeqIndex);

			AppendUniqueSyncmerHashes(Seq, L, HashesVec[SpeciesIndex]);
			}
		}

	Progress("Initialize counts\n");
	m_Sizes.resize(m_Slots, 0);
	vector<vector<bool> > SlotFoundVec(m_SpeciesCount);
	for (uint SpeciesIndex = 0; SpeciesIndex < m_SpeciesCount; ++SpeciesIndex)
		{
		ProgressStep(SpeciesIndex, m_SpeciesCount, "Counting");
		vector<bool> &SlotFound = SlotFoundVec[SpeciesIndex];
		SlotFound.resize(m_Slots, false);
		const set<uint64> &Hashes = HashesVec[SpeciesIndex];
		for (set<uint64>::const_iterator p = Hashes.begin();
		  p != Hashes.end(); ++p)
			{
			uint64 Hash = *p;
			uint32 Slot = uint32(Hash%m_Slots);
			SlotFound[Slot] = true;
			}
		for (uint32 Slot = 0; Slot < m_Slots; ++Slot)
			if (SlotFound[Slot])
				++(m_Sizes[Slot]);
		}

	ProgressLog("Count dist\n");
	vector<uint> CountDist(m_SpeciesCount);
	for (uint32 Slot = 0; Slot < m_Slots; ++Slot)
		{
		uint Count = m_Sizes[Slot];
		asserta(Count < m_SpeciesCount);
		++(CountDist[Count]);
		}

	for (uint i = 0; i < m_SpeciesCount; ++i)
		Log("%u  %10u\n", i, CountDist[i]);

	m_Offsets.reserve(m_Slots+1);
	m_Sizes.resize(m_Slots, 0);

	uint32 Offset = 0;
	for (uint32 Slot = 0; Slot < m_Slots; ++Slot)
		{
		uint16 Size = m_Sizes[Slot];
		m_Offsets.push_back(Offset);
		Offset += Size;
		}
	m_TotalSize = Offset;
	m_Offsets.push_back(m_TotalSize);

	m_SpVec.resize(m_TotalSize);
	vector<uint16> Sizes2(m_Slots);

	for (uint SpeciesIndex = 0; SpeciesIndex < m_SpeciesCount; ++SpeciesIndex)
		{
		ProgressStep(SpeciesIndex, m_SpeciesCount, "Building");
		uint16 Sp = uint16(SpeciesIndex);
		const vector<bool> &SlotFound = SlotFoundVec[SpeciesIndex];
		asserta(SIZE(SlotFound) == m_Slots);
		for (uint32 Slot = 0; Slot < m_Slots; ++Slot)
			{
			if (SlotFound[Slot])
				{
				asserta(Slot < SIZE(m_Offsets));
				uint32 Offset = m_Offsets[Slot];
				asserta(Slot < SIZE(Sizes2));
				uint32 Size2 = Sizes2[Slot];
				++(Sizes2[Slot]);
				uint32 Ix = Offset + Size2;
				asserta(Ix < SIZE(m_SpVec));
				m_SpVec[Ix] = Sp;
				}
			}
		}

	if (opt(validate))
		{
		ProgressLog("Validating...");
		Validate(HashesVec);
		ProgressLog(" ok\n");
		}

	Progress("Checking sizes...");
	for (uint32 Slot = 0; Slot < m_Slots; ++Slot)
		asserta(m_Sizes[Slot] <= m_Sizes[Slot]);
	Progress(" done.\n");

	ProgressLog("%u species\n", m_SpeciesCount);
	ProgressLog("%s slots\n", Int64ToStr(m_Slots));
	ProgressLog("%s total size\n", Int64ToStr(m_TotalSize));
	}

void SpeciesIndex::ToFile(const string &FileName) const
	{
	if (FileName.empty())
		return;

	FILE *fOut = CreateStdioFile(FileName);

	const uint32 m_SpeciesCount = SIZE(m_SpeciesNames);
	WriteStdioFile(fOut, &m_ST, sizeof(m_ST));
	WriteStdioFile(fOut, &m_k, sizeof(m_k));
	WriteStdioFile(fOut, &m_t, sizeof(m_t));
	WriteStdioFile(fOut, &m_w, sizeof(m_w));
	WriteStdioFile(fOut, &m_SpeciesCount, sizeof(m_SpeciesCount));
	WriteStdioFile(fOut, &m_Slots, sizeof(m_Slots));
	WriteStdioFile(fOut, &m_TotalSize, sizeof(m_TotalSize));

	for (uint i = 0; i < m_SpeciesCount; ++i)
		{
		const string &Name = m_SpeciesNames[i];
		uint32 n = SIZE(Name);
		WriteStdioFile(fOut, &n, sizeof(n));
		WriteStdioFile(fOut, Name.c_str(), n);
		}

	uint32 Bytes;
	//Bytes = Slots*sizeof(m_Sizes[0]);
	//ProgressLog("Writing sizes (%s bytes)\n", MemBytesToStr(Bytes));
	//WriteStdioFile(fOut, m_Sizes.data(), Bytes);

	Bytes = (m_Slots + 1)*sizeof(m_Offsets[0]);
	asserta(m_Offsets[m_Slots] == m_TotalSize);
	ProgressLog("Writing offsets (%s bytes)\n", MemBytesToStr(Bytes));
	WriteStdioFile(fOut, m_Offsets.data(), Bytes);

	Bytes = m_TotalSize*sizeof(m_SpVec[0]);
	ProgressLog("Writing spvecs (%s bytes)\n", MemBytesToStr(Bytes));
	WriteStdioFile(fOut, m_SpVec.data(), Bytes);

	Progress("Done.\n");

	CloseStdioFile(fOut);
	}

void SpeciesIndex::FromFile(const string &FileName)
	{
	Clear();

	FILE *f = OpenStdioFile(FileName);

	ReadStdioFile(f, &m_ST, sizeof(m_ST));
	ReadStdioFile(f, &m_k, sizeof(m_k));
	ReadStdioFile(f, &m_t, sizeof(m_t));
	ReadStdioFile(f, &m_w, sizeof(m_w));

	ReadStdioFile(f, &m_SpeciesCount, sizeof(m_SpeciesCount));
	ReadStdioFile(f, &m_Slots, sizeof(m_Slots));
	ReadStdioFile(f, &m_TotalSize, sizeof(m_TotalSize));
	ProgressLog("%u species\n", m_SpeciesCount);
	ProgressLog("%s slots\n", Int64ToStr(m_Slots));
	ProgressLog("%s total size\n", Int64ToStr(m_TotalSize));
	for (uint i = 0; i < m_SpeciesCount; ++i)
		{
		uint32 n;
		ReadStdioFile(f, &n, sizeof(n));
		string Name;
		for (uint j = 0; j < n; ++j)
			{
			char c;
			ReadStdioFile(f, &c, 1);
			Name += c;
			}
		// Log("[%2u] %s\n", i, Name.c_str());
		m_SpeciesNames.push_back(Name);
		}

	m_Sizes.resize(m_Slots);
	m_Offsets.resize(m_Slots+1);
	m_SpVec.resize(m_TotalSize);

	uint32 Bytes;
	//Bytes = m_Slots*sizeof(m_Sizes[0]);
	//Progress("Reading sizes (%s bytes)\n", MemBytesToStr(Bytes));
	//ReadStdioFile(f, m_Sizes.data(), Bytes);

	Bytes = (m_Slots + 1)*sizeof(m_Offsets[0]);
	Progress("Reading offsets (%s bytes)\n", MemBytesToStr(Bytes));
	ReadStdioFile(f, m_Offsets.data(), Bytes);

	Bytes = m_TotalSize*sizeof(m_SpVec[0]);
	Progress("Reading spvecs (%s bytes)\n", MemBytesToStr(Bytes));
	ReadStdioFile(f, m_SpVec.data(), Bytes);

	Progress("Done.\n");

	CloseStdioFile(f);
	}

void cmd_species_index()
	{
	const string &InputFileName = opt(species_index);
	const string &OutputFileName = opt(output);

	SyncmerType ST = ST_Syncmer1;
	uint32 k = 14;
	uint32 t = 50;
	uint32 w = 0;
	uint32 Slots = 0;
	double LoadFactor = 0.6;

	if (optset_st)
		ST = StrToST(opt(st));
	if (optset_k)
		k = opt(k);
	if (optset_t)
		t = opt(t);
	if (optset_w)
		w = opt(w);
	if (optset_slots)
		Slots = opt(slots);
	if (optset_load_factor)
		LoadFactor = opt(load_factor);

	SeqDB DB;
	DB.FromFasta(InputFileName);
	uint64 DBSize = DB.GetLetterCount();

	if (optset_slots)
		Slots = opt(slots);
	else
		{
		uint64 EstimatedSyncmerCount = DBSize/t;
		int64 iMinSlotCount = int64(EstimatedSyncmerCount/LoadFactor);
		uint64 GetPrime(uint64 n);
		uint64 Slots64 = GetPrime(iMinSlotCount);
		asserta(Slots64 < UINT32_MAX/2);
		Slots = uint32(Slots64);
		}

	ProgressLog("\n");
	ProgressLog("    ST  %s\n", STToStr(ST));
	ProgressLog("     k  %u\n", k);
	ProgressLog("     t  %u\n", t);
	ProgressLog("     w  %u\n", w);
	ProgressLog("    DB  %u (%s)\n", DBSize, MemBytesToStr(DBSize));
	ProgressLog(" Slots  %u (%s)\n", Slots, MemBytesToStr(Slots));
	ProgressLog("  Load  %.2f\n", LoadFactor);
	ProgressLog("\n");

	SpeciesIndex SPI;
	SPI.FromSeqDB(ST, k, t, w, Slots, DB);
	SPI.ToFile(OutputFileName);
	}
