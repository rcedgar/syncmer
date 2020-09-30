#include "myutils.h"
#include "syncmerindex.h"
#include "seqdb.h"
#include "fastaseqsource.h"
#include "seqinfo.h"
#include "objmgr.h"
#include "sort.h"

void GetSpeciesToSeqIndexes(const SeqDB &DB,
  vector<string> &SpeciesNames,
  vector<vector<uint> > &SeqIndexesVec);

static void Align(const SeqDB &DB, const SyncmerIndex &QSI,
  const SyncmerIndex &TSI)
	{
	const uint SyncmerCount = SIZE(QSI.m_SyncmerCoords);
	vector<uint> Diags;
	for (uint i = 0; i < SyncmerCount; ++i)
		{
		uint32 Coord = QSI.m_SyncmerCoords[i];
		asserta(Coord < QSI.m_L);
		uint64 Hash = QSI.m_Hashes[Coord];
		uint32 Slot = (Hash%TSI.m_SlotCount);
		uint32 TargetCoord = TSI.m_HashTable[Slot];
		if (TargetCoord != UINT32_MAX)
			{
			uint Base = TSI.m_L + TargetCoord;
			asserta(Coord <= Base);
			uint Diag = Base - Coord;
			Diags.push_back(Diag);
			}
		}
	}

void cmd_species_align()
	{
	const string &QueryFileName = opt(species_align);
	const string &DBFileName = opt(db);

	SeqDB DB;
	DB.FromFasta(DBFileName);

	FASTASeqSource FSS;
	FSS.Open(QueryFileName);

	SeqInfo *Query = ObjMgr::GetSeqInfo();

	SyncmerType ST = ST_Syncmer4;
	uint32 k = 14;
	uint32 t = 9;
	uint32 w = 0;
	double LoadFactor = 0.6;

	if (optset_st)
		ST = StrToST(opt(st));
	if (optset_k)
		k = opt(k);
	if (optset_t)
		t = opt(t);
	if (optset_w)
		w = opt(w);
	if (optset_load_factor)
		LoadFactor = opt(load_factor);

	ProgressLog("\n");
	ProgressLog("    ST  %s\n", STToStr(ST));
	ProgressLog("     k  %u\n", k);
	ProgressLog("     t  %u\n", t);
	ProgressLog("     w  %u\n", w);
	ProgressLog("  Load  %.2f\n", LoadFactor);
	ProgressLog("\n");

	vector<string> SpeciesNames;
	vector<vector<uint> > SeqIndexesVec;
	vector<SyncmerIndex *> SIs;
	GetSpeciesToSeqIndexes(DB, SpeciesNames, SeqIndexesVec);

	const uint SpeciesCount = SIZE(SpeciesNames);
	for (uint SpeciesIndex = 0; SpeciesIndex < SpeciesCount; ++SpeciesIndex)
		{
		ProgressStep(SpeciesIndex, SpeciesCount, "Building SIs %s",
		  SpeciesNames[SpeciesIndex].c_str());

		SyncmerIndex *SI = new SyncmerIndex;
		const vector<uint> &SeqIndexes = SeqIndexesVec[SpeciesIndex];
		const unsigned n = SIZE(SeqIndexes);
		uint SumL = 0;
		if (n == 1)
			{
			const byte *Seq = DB.GetSeq(SeqIndexes[0]);
			SumL = DB.GetSeqLength(SeqIndexes[0]);
			SI->Create(ST, k, t, w, Seq, SumL);
			}
		else
			{
			for (uint i = 0; i < n; ++i)
				SumL += DB.GetSeqLength(SeqIndexes[i]);
			byte *CatSeq = myalloc(byte, SumL);
			uint Offset = 0;
			for (uint i = 0; i < n; ++i)
				{
				const byte *Seq = DB.GetSeq(SeqIndexes[i]);
				const uint L = DB.GetSeqLength(SeqIndexes[i]);
				memcpy(CatSeq + Offset, Seq, L);
				Offset += L;
				}
			asserta(Offset == SumL);
			SI->Create(ST, k, t, w, CatSeq, SumL);
			}
		uint SlotCount = uint(SumL/LoadFactor);
		SI->SetHashTable(SlotCount);
		SIs.push_back(SI);
		}

	while (FSS.GetNext(Query))
		{
		SyncmerIndex QSI;
		const char *QLabel = Query->m_Label;
		const byte *Q = Query->m_Seq;
		uint QL = Query->m_L;
		QSI.Create(ST, k, t, w, QLabel, Q, QL);
		//uint SlotCount = uint(QL/LoadFactor);
		for (uint SpeciesIndex = 0; SpeciesIndex < SpeciesCount; ++SpeciesIndex)
			{
			const SyncmerIndex &TSI = *SIs[SpeciesIndex];
			Align(DB, QSI, TSI);
			}
		}
	}
