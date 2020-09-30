#include "myutils.h"
#include "spindex.h"

void cmd_make_spindex()
	{
// Text file with one FASTA filename per line
	const string &ListFileName = opt(make_spindex);
	string InputDir = (optset_inputdir ? opt(inputdir) : ".");
	if (!EndsWith(InputDir, "/") && !EndsWith(InputDir, "\\"))
		InputDir += "/";

	Spindex SI;
	SI.Init(13);

	string FastaFileName;
	vector<string> FastaFileNames;
	FILE *f = OpenStdioFile(ListFileName);
	while (ReadLineStdioFile(f, FastaFileName))
		FastaFileNames.push_back(FastaFileName);
	CloseStdioFile(f);

	const uint N = SIZE(FastaFileNames);
	for (uint i = 0; i < N; ++i)
		{
		const string &FastaFileName = FastaFileNames[i];
		ProgressStep(i, N, "Pass 1 %s", FastaFileName.c_str());
		SI.AddFile_Pass1(InputDir + FastaFileName);
		}

	SI.LogStats();
	}
