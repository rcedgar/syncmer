#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

void TestSubmers();
void TestMinimizers();
void TestDist();

const char *Uint64ToBinaryStr(uint64 i, string &s)
	{
	s.clear();
	while (i != 0)
		{
		uint bit = (i & 1);
		s.push_back(bit ? '1' : '0');
		i >>= 1;
		}
	reverse(s.begin(), s.end());
	return s.c_str();
	}

uint64 rot64(uint64 i, uint64 n)
	{
	return (i << n) | (i >> (64-n));
	}

void cmd_test()
	{
	opt(test);
	TestSubmers();
	}
