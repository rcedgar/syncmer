#include "myutils.h"
#include "seqdb.h"
#include "linereader.h"
#include "alpha.h"
#include "cigar.h"

static void WriteAln(FILE *f, bool Plus,
  const string &QLabel, const byte *Q, uint QLo, uint QHi, uint QL,
  const string &RLabel, const byte *R, uint RLo, uint RHi, uint RL,
  const vector<char> &Ops, const vector<uint> &Lengths)
	{
	if (f == 0)
		return;

	const uint OpCount = SIZE(Ops);
	asserta(SIZE(Lengths) == OpCount);

	string QRow;
	string RRow;
	int QPos = int(QLo);
	int RPos = int(Plus ? RLo : RHi);
	for (uint i = 0; i < OpCount; ++i)
		{
		char Op = Ops[i];
		uint n = Lengths[i];
		switch (Op)
			{
		case 'M':
			{
			for (uint j = 0; j < n; ++j)
				{
				asserta(QPos < int(QL));
				asserta(RPos >= 0 && RPos < int(RL));
				QRow += Q[QPos++];
				if (Plus)
					RRow += R[RPos++];
				else
					{
					byte r = R[RPos];
					RRow += g_CharToCompChar[r];
					--RPos;
					}
				}
			break;
			}
		case 'D':
			{
			for (uint j = 0; j < n; ++j)
				{
				asserta(RPos >= 0 && RPos < int(RL));
				QRow += '-';
				if (Plus)
					RRow += R[RPos++];
				else
					{
					byte r = R[RPos];
					RRow += g_CharToCompChar[r];
					--RPos;
					}
				}
			break;
			}
		case 'I':
			{
			for (uint j = 0; j < n; ++j)
				{
				asserta(QPos < int(QL));
				QRow += Q[QPos++];
				RRow += '-';
				}
			break;
			}

		default:
			Die("Unknown op '%c'", Op);
			}
		}
	fprintf(f, "\n");
	fprintf(f, ">%s  Q:%u-%u(%u), R:%u-%u(%u) %c\n",
	  QLabel.c_str(),
	  QLo, QHi, QHi - QLo + 1,
	  RLo, RHi, RHi - RLo + 1,
	  pom(Plus));
	//fprintf(f, "Q: %s\n", QRow.c_str());
	//fprintf(f, "R: %s\n", RRow.c_str());

	const uint ColCount = SIZE(QRow);
	const uint COLS = 100;
	for (uint Col = 0; Col < ColCount; Col += COLS)
		{
		uint k = COLS;
		if (Col + k > ColCount)
			k = ColCount - Col;
		if (Col > 0)
			fprintf(f, "\n");
		fprintf(f, "Q: %*.*s\n", k, k, QRow.c_str() + Col);
		fprintf(f, "R: %*.*s\n", k, k, RRow.c_str() + Col);
		}
	}

void cmd_paf2features()
	{
	const string &PAFFileName = opt(paf2aln);
	const string &QueryFileName = opt(input);
	const string &DBFileName = opt(db);
	const string &OutputFileName = opt(output);

	FILE *fOut = CreateStdioFile(OutputFileName);

	SeqDB Query;
	SeqDB DB;

	Query.FromFastx(QueryFileName);
	DB.FromFasta(DBFileName);

	vector<string> Fields;
	t_LineBuff LB;
	LineReader LR;
	LR.Open(PAFFileName);
	vector<char> Ops;
	vector<uint> Lengths;
	while (LR.ReadLine(LB))
		{
		const string Line = string(LB.Data);
		Split(Line, Fields, '\t');

		const string &QLabel = Fields[0];
		uint QSeqIndex = Query.GetSeqIndex(QLabel);
		uint QL = StrToUint(Fields[1]);
		uint QLo = StrToUint(Fields[2]);
		uint QHi = StrToUint(Fields[3]) - 1;
		asserta(QLo < QHi);
		asserta(QHi < QL);

		const string &RLabel = Fields[5];
		uint RSeqIndex = DB.GetSeqIndex(RLabel);
		uint RL = StrToUint(Fields[6]);
		uint RLo = StrToUint(Fields[7]);
		uint RHi = StrToUint(Fields[8]) - 1;
		asserta(RLo < RHi);
		asserta(QHi < RL);

		bool Plus;
		const string &Strand = Fields[4];
		if (Strand == "+")
			Plus = true;
		else if (Strand == "-")
			Plus = false;
		else
			Die("Bad strand '%s'", Strand.c_str());

	//  s2:i:0  de:f:0.1493     rl:i:30 cg:Z:13M2I33M2D33M1D10M1...
		string CIGAR;
		for (uint i = 12; i < SIZE(Fields); ++i)
			{
			const string &Field = Fields[i];
			if (StartsWith(Field, "cg:Z:"))
				{
				CIGAR = Field.substr(5, string::npos);
				break;
				}
			}
		if (CIGAR.empty())
			Die("Missing CIGAR");

		CIGARGetOps(CIGAR, Ops, Lengths);
		const uint OpCount = SIZE(Ops);
		asserta(SIZE(Lengths) == OpCount);

		const byte *Q = Query.GetSeq(QSeqIndex);
		const byte *R = DB.GetSeq(RSeqIndex);

		uint QL2 = Query.GetSeqLength(QSeqIndex);
		uint RL2 = DB.GetSeqLength(RSeqIndex);
		asserta(QL2 == QL);
		asserta(RL2 == RL);

		WriteAln(g_fLog, Plus,
		  QLabel, Q, QLo, QHi, QL,
		  RLabel, R, RLo, RHi, RL,
		  Ops, Lengths);
		}

	CloseStdioFile(fOut);
	}
