#include "myutils.h"
#include "cigar.h"

void PathToCIGAR(const char *Path, string &CIGAR)
	{
	char LastC = *Path;
	unsigned n = 1;
	char Tmp[32];
	for (unsigned i = 1; ; ++i)
		{
		char c = Path[i];
		if (c == 0)
			break;
		if (c == LastC)
			{
			++n;
			continue;
			}
		else
			{
			assert(n > 0);
			if (LastC == 'D')
				LastC = 'I';
			else if (LastC == 'I')
				LastC = 'D';
			sprintf(Tmp, "%u%c", n, LastC);
			CIGAR += string(Tmp);
			LastC = c;
			n = 1;
			}
		}
	if (n > 0)
		{
		if (LastC == 'D')
			LastC = 'I';
		else if (LastC == 'I')
			LastC = 'D';
		sprintf(Tmp, "%u%c", n, LastC);
		CIGAR += string(Tmp);
		}
	}

void CIGARGetOps(const string &CIGAR, vector<char> &Ops,
  vector<unsigned> &Lengths)
	{
	Ops.clear();
	Lengths.clear();
	if (CIGAR.empty())
		return;

	unsigned L = SIZE(CIGAR);
	unsigned n = 0;
	for (unsigned i = 0; i < L; ++i)
		{
		char c = CIGAR[i];
		if (isdigit(c))
			n = n*10 + (c - '0');
		else if (isupper(c) || c == '=')
			{
			if (n == 0)
				Die("Operation '%c' has zero length in CIGAR '%s'", c, CIGAR.c_str());
			Ops.push_back(c);
			Lengths.push_back(n);
			n = 0;
			}
		else
			Die("Invalid char '%c' in CIGAR '%s'", c, CIGAR.c_str());
		}
	if (n > 0)
		Die("Missing operation at end of CIGAR '%s'", CIGAR.c_str());
	}

unsigned CIGARToQL(const string &CIGAR)
	{
	vector<char> Ops;
	vector<unsigned> Lengths;
	CIGARGetOps(CIGAR, Ops, Lengths);
	const unsigned N = SIZE(Ops);
	asserta(SIZE(Lengths) == N);
	unsigned QL = 0;
	for (unsigned i = 0; i < N; ++i)
		{
		char Op = Ops[i];
		switch (Op)
			{
		case 'M':
		case 'I':
			QL += Lengths[i];
			break;

		case 'D':
			break;

		default:
			Die("Unsupported op '%c' in CIGAR '%s'", Op, CIGAR.c_str());
			}
		}
	return QL;
	}
