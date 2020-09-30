#include "myutils.h"
#include "filetype.h"
#include "alpha.h"

bool FastaFileIsNucleo(FILE *f)
	{
	unsigned SampleSize = 1024;
	uintB CurrPos = GetStdioFilePosB(f);
	uintB FileSize = GetStdioFileSizeB(f);

	SetStdioFilePos64(f, 0);
	byte lastc = '\n';
	bool InLabel = false;
	unsigned NucleoCount = 0;
	unsigned LetterCount = 0;
	unsigned UpperCount = 0;
	for (uint64 Pos = 0; Pos < FileSize; ++Pos)
		{
		byte c;
		ReadStdioFile(f, &c, 1);
		if (c == '\r')
			continue;
		if (c == '>' && (lastc == '\n'))
			InLabel = true;
		else if (InLabel && (c == '\n'))
			InLabel = false;
		else if (!InLabel && isalpha(c))
			{
			if (isupper(c))
				++UpperCount;
			++LetterCount;
			if (g_IsNucleoChar[c])
				++NucleoCount;
			if (LetterCount >= SampleSize)
				break;
			}
		lastc = c;
		}

	bool IsNucleo = (LetterCount > 0 && double(NucleoCount)/double(LetterCount) > 0.9);
	return IsNucleo;
	}

FILE_TYPE GetFileType(const string &FileName, bool *ptrNucleo)
	{
	FILE_TYPE Type = FT_Unknown;
	*ptrNucleo = false;

	FILE *f = OpenStdioFile(FileName);
	uintB FileSize = GetStdioFileSizeB(f);
	if (FileSize == 0)
		Die("Empty file %s", FileName.c_str());

	byte b;
	ReadStdioFile(f, &b, 1);

	if (b == '>')
		{
		Type = FT_FASTA;
		*ptrNucleo = FastaFileIsNucleo(f);
		}
	else if (b == '@')
		{
		Type = FT_FASTQ;
		*ptrNucleo = true;
		}
	CloseStdioFile(f);

	if (Type == FT_Unknown)
		Die("Unknown file format %s", FileName.c_str());

	return Type;
	}
