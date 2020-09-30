#include "myutils.h"

#define TRACE	0

#if TRACE
static const char *WordToStr(uint64 Word, uint n)
	{
	static char Tmp[128];
	for (uint i = 0; i < n; ++i)
		{
//		uint BitsToShift = 2*(n - i - 1);
		uint BitsToShift = 2*i;
		const uint64 ShiftedWord = (Word >> BitsToShift);
		uint c = "ACGT"[ShiftedWord & 0x3];
		Tmp[i] = c;
		}
	Tmp[n] = 0;
	return Tmp;
	}
#endif

static uint64 GetSub(uint64 Word, uint s, uint Pos)
	{
	uint64 Mask = 0;
	for (uint i = 0; i < 2*s; ++i)
		Mask |= (uint64(1) << i);

	uint64 Shifted = (Word >> (2*Pos));
	uint64 Sub = (Shifted & Mask);
	return Sub;
	}

static uint GetMinPos(uint64 Word, uint k, uint s)
	{
	for (uint i = 0; i < k - s + 1; ++i)
		{
		}
	}

//static void Test(uint64 Word, uint k, uint s, uint Pos)
//	{
//	uint64 Sub = GetSub(Word, s, Pos);
//	Log("Word %08" PRIx64 " %s", Word, WordToStr(Word, k));
//	Log("  Sub Pos=%u %" PRIx64 " %s", Pos, Sub, WordToStr(Sub, s));
//	Log("\n");
//	}

void cmd_dict()
	{
	//Test(9, 4, 2, 0);
	//Test(9, 4, 2, 1);
	//Test(9, 4, 2, 2);
	//return;

	const string &OutputFileName = opt(dict);
	asserta(optset_k);
	asserta(optset_s);

	const uint k = opt(k);
	const uint s = opt(s);

	const uint HiPos = k - s + 1;
	const uint64 HiWord = myipow64(4, k);
	uint64 Mask = 0;
	for (uint64 Bit = 0; Bit < 2*s; ++Bit)
		Mask |= (uint64(1) << Bit);

#if TRACE
	uint64 FromWord = 0;
	uint64 ToWord = FromWord + 100;
#else
	uint64 FromWord = 0;
	uint64 ToWord = HiWord;
#endif
	uint64 DictSize = 0;
	for (uint64 Word = FromWord; Word < ToWord; ++Word)
		{
		uint64 MinSub = 0;
		uint MinPos = 0;
		for (uint Pos = 0; Pos < HiPos; ++Pos)
			{
			uint64 Sub = (Word >> (2*Pos)) & Mask;
			if (Pos == 0 || Sub < MinSub)
				{
				MinPos = Pos;
				MinSub = Sub;
				}
			}
		if (MinPos == 0 || MinPos == k - s)
			++DictSize;

#if TRACE
		{
		Log("Word %08" PRIx64 " %s", Word, WordToStr(Word, k));
		Log("  Sub Pos=%u %" PRIx64 " %s", MinPos, MinSub, WordToStr(MinSub, s));
		Log("\n");
		}
#endif
		}


	uint w = k - s;
	uint L = k + w - 1;
	double c = double(HiWord)/DictSize;

	ProgressLog("k=%u", k);
	ProgressLog(" s=%u", s);
	ProgressLog(" w=%u", w);
	ProgressLog(" L=%u", L);
	ProgressLog(" c=%.2f", c);
	ProgressLog(" (w+1)/2=%.2f", (w + 1)/2.0);
	ProgressLog(" density=%.3g", 1.0/c);
	ProgressLog(" size=%.3g", double(DictSize));
	ProgressLog("\n");
	}
