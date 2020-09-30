#include "myutils.h"
#include "searchsix.h"
#include "seqdb.h"
#include "objmgr.h"
#include "alpha.h"

FILE *SearchHashix::m_fTab;

#include "myutils.h"
#include "searchsix.h"
#include "seqdb.h"
#include "objmgr.h"
#include "alpha.h"

void SearchHashix::AppendGSP(SHGSP *GSP)
	{
	m_GSPs.push_back(GSP);
	}

SHGSP *SearchHashix::Extend(uint32 SeedPosQ, uint32 SeedPosT)
	{
	const byte *Q = m_Query->m_Seq;
	const byte *T = m_SI->m_Seq;
	const unsigned QL = m_Query->m_L;
	const unsigned TL = m_SI->m_L;

	const unsigned w = m_WordLength;

	uint IdCount = 0;
	uint BestIdCount = 0;

	float Score = 0.0f;
	float BestScore = 0.0f;

// Extend right
	uint QPos = SeedPosQ;
	uint TPos = SeedPosT;
	uint EndPosQ = QPos;
	for (;;)
		{
		byte t = T[TPos++];
		if (t == 0)
			break;

		if (QPos >= QL)
			break;
		byte q = Q[QPos];

		if (q == t)
			{
			++IdCount;
			Score += 1.0f;
			if (Score > BestScore)
				{
				BestIdCount = IdCount;
				BestScore = Score;
				EndPosQ = QPos;
				}
			}
		else
			{
			Score += MISMATCH_SCORE;
			float Drop = (BestScore - Score);
			if (Drop > XDROP)
				break;
			}
		++QPos;
		}

// Extend left
	QPos = SeedPosQ;
	TPos = SeedPosT;
	uint StartPosQ = SeedPosQ;
	for (;;)
		{
		// Not needed coz nuls
		//if (TPos == 0)
		//	break;
		byte t = T[--TPos];
		if (t == 0)
			break;

		if (QPos == 0)
			break;
		byte q = Q[QPos];
		if (q == t)
			{
			++IdCount;
			Score += 1.0f;
			if (Score > BestScore)
				{
				BestIdCount = IdCount;
				BestScore = Score;
				StartPosQ = QPos;
				}
			}
		else
			{
			Score += MISMATCH_SCORE;
			float Drop = (BestScore - Score);
			if (Drop > XDROP)
				break;
			}
		--QPos;
		}

	if (BestScore < MIN_GSP_SCORE)
		return 0;

	assert(EndPosQ > StartPosQ);
	assert(StartPosQ <= SeedPosQ);
	uint GSPLength = EndPosQ - StartPosQ + 1;
	if (GSPLength < MIN_GSP_LENGTH)
		return 0;

	SHGSP *GSP = new SHGSP;
	uint StartPosT = SeedPosT - (SeedPosQ - StartPosQ);
	uint EndPosT = StartPosT + GSPLength - 1;
	asserta(EndPosT > StartPosT);

	GSP->StartQ = StartPosQ;
	GSP->EndQ = EndPosQ;
	GSP->LoT = StartPosT;
	GSP->HiT = EndPosT;
	GSP->IdCount = BestIdCount;
	GSP->Score = BestScore;

	return GSP;
	}

void SearchHashix::ClearGSPs()
	{
	for (vector<SHGSP *>::const_iterator p = m_GSPs.begin();
	  p != m_GSPs.end(); ++p)
		delete *p;
	m_GSPs.clear();
	}

void SearchHashix::LogGSPs() const
	{
	Log("\n");
	Log("%u GPS\n", SIZE(m_GSPs));
	uint64 TotalBases = 0;
	double TotalBasesPctId = 0;
	for (vector<SHGSP *>::const_iterator p = m_GSPs.begin();
	  p != m_GSPs.end(); ++p)
		{
		const SHGSP &GSP = **p;

		uint Len = GSP.GetLength();
		double PctId = GSP.GetPctId();

		TotalBases += Len;
		TotalBasesPctId += Len*PctId;

		bool Plus = (GSP.HiT < m_SI->m_L/2);
		Log("%10u", GSP.GetLoQ());
		Log("  %10u", GSP.GetHiQ());
		Log("  %10u", GSP.LoT);
		Log("  %c", pom(Plus));
		Log("  %10u", GSP.HiT);
		Log("  %8u", Len);
		Log("  %6.1f", PctId);
		Log("\n");
		}
	Log("\n");
	Log("ANI %.1f\n", TotalBasesPctId/TotalBases);
	Log("AFQ %.2f\n", double(TotalBases)/m_Query->m_L);
	}

void SearchHashix::GetANI_AFQ(double &ANI, double &AFQ) const
	{
	uint64 TotalBases = 0;
	double TotalBasesPctId = 0;
	for (vector<SHGSP *>::const_iterator p = m_GSPs.begin();
	  p != m_GSPs.end(); ++p)
		{
		const SHGSP &GSP = **p;

		uint Len = GSP.GetLength();
		double PctId = GSP.GetPctId();

		TotalBases += Len;
		TotalBasesPctId += Len*PctId;
		}
	if (TotalBases == 0)
		{
		ANI = 0;
		AFQ = 0;
		}
	else
		{
		ANI = TotalBasesPctId/TotalBases;
		AFQ = double(TotalBases)/m_Query->m_L;
		}
	}

void SearchHashix::Search(SeqInfo *Query)
	{
	ClearGSPs();
	m_Query = Query;

	const byte *QSeq = Query->m_Seq;
	const uint QL = Query->m_L;
	SyncmerType ST = m_SI->m_ST;
	uint k = m_SI->m_k;
	uint t = m_SI->m_t;
	uint w = m_SI->m_w;

	m_QSI.Clear();
	m_QSI.Create(ST, k, t, w, QSeq, QL);

	const uint K = m_QSI.GetKmerCount();
	for (uint QPos = 0; QPos < K; ++QPos)
		{
		if (m_QSI.IsSyncmer(QPos))
			{
			uint64 Kmer = m_QSI.m_Kmers[QPos];
			uint64 Hash = m_QSI.KmerToHash(Kmer);
			uint64 Slot = Hash%m_SI->m_SlotCount;
			uint32 TPos = m_SI->GetPos(Slot);
			if (TPos != UINT32_MAX)
				{
				SHGSP *GSP = Extend(QPos, TPos);
				if (GSP != 0)
					{
					AppendGSP(GSP);
					QPos += GSP->GetLength();
					goto Next;
					}
				}
			}
	Next:;
		}
	}

void cmd_searchsix()
	{
	const string &QueryFileName = opt(searchsix);
	const string &DBFileName = opt(db);
	const string &STSStr = opt(sts);
	vector<string> STNames;
	Split(STSStr, STNames, '+');

	asserta(optset_k);
	asserta(optset_t || optset_w);
	uint k = opt(k);
	uint t = (optset_t ? opt(t) : UINT_MAX);
	uint w = (optset_w ? opt(w) : UINT_MAX);

	vector<SyncmerType> STs;
	for (uint i = 0; i < SIZE(STNames); ++i)
		{
		SyncmerType ST = StrToST(STNames[i]);
		STs.push_back(ST);
		}

	SeqDB DB;
	DB.FromFasta(DBFileName);

	SeqDB QDB;
	QDB.FromFasta(QueryFileName);

	ProgressLog("@ANI");
	for (uint i = 0; i < SIZE(STNames); ++i)
		{
		SyncmerType ST = STs[i];
		SyncmerIndex SI;
		SI.FromSeqDB(ST, DB, k, t, w, 0.3);

		SearchHashix SH;
		SH.SetHI(&SI);

		SeqInfo *Query = ObjMgr::GetSeqInfo();
		const unsigned QuerySeqCount = QDB.GetSeqCount();
		uint64 TotalBases= 0;
		double TotalANIBases = 0.0;
		double TotalAFQBases = 0.0;
		for (unsigned QuerySeqIndex = 0; QuerySeqIndex < QuerySeqCount; ++QuerySeqIndex)
			{
			QDB.GetSI(QuerySeqIndex, *Query);
			SH.Search(Query);
			double ANI, AFQ;
			SH.GetANI_AFQ(ANI, AFQ);

			uint QL = Query->m_L;
			TotalBases += QL;
			TotalANIBases += ANI*QL;
			TotalAFQBases += AFQ*QL;
			}

		const char *Name = STNames[i].c_str();
		double ANI = TotalANIBases/TotalBases;
		double AFQ = TotalAFQBases/TotalBases;
		double IndexSizeKilo = SI.GetIndexSize()/1000.0;
		ProgressLog("	%s.ANI=%.1f	%s.AFQ=%.4f %s.Ix=%.1f",
		  Name, ANI, Name, AFQ, Name, IndexSizeKilo);
		}
	ProgressLog("\n");
	}
