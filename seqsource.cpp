#include "myutils.h"
#include "fastaseqsource.h"
#include "fastqseqsource.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "filetype.h"
#include "seqdb.h"

#define TIME_LOCKS	0

#if TIME_LOCKS
#include "getticks.h"
static TICKS g_tLocks;
static TICKS g_tUnLocks;
#endif

SeqSource::SeqSource()
	{
	m_SI = ObjMgr::GetSeqInfo();
	m_SeqCount = 0;
	m_DoGetLock = true;
	}

SeqSource::~SeqSource()
	{
	ObjMgr::Down(m_SI);
	}

bool SeqSource::GetNext(SeqInfo *SI)
	{
	if (m_DoGetLock)
		{
#if	TIME_LOCKS
		TICKS t1 = GetClockTicks();
#endif
		LOCK_CLASS();
#if	TIME_LOCKS
		TICKS t2 = GetClockTicks();
		g_tLocks += (t2 - t1);
#endif
		}
	bool Ok = GetNextLo(SI);
	if (m_DoGetLock)
		{
#if	TIME_LOCKS
		TICKS t1 = GetClockTicks();
#endif
		UNLOCK_CLASS();
#if	TIME_LOCKS
		TICKS t2 = GetClockTicks();
		g_tUnLocks += (t2 - t1);
#endif
		}

	if (!Ok)
		{
#if	TIME_LOCKS
		Log("SeqSource locks %.3e, unlocks %.3e\n", double(g_tLocks), double(g_tUnLocks));
#endif
		return false;
		}

	++m_SeqCount;
	return true;
	}
