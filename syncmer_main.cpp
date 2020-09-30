#include "myutils.h"
#include "seqdb.h"

int main(int argc, char **argv)
	{
#ifdef _MSC_VER
	_putenv("TZ=");
#endif
	setbuf(stdout, 0);
	setbuf(stderr, 0);

	MyCmdLine(argc, argv);

	if (!opt(quiet))
		{
		PrintProgramInfo(stdout);
		PrintCopyright(stdout);
		}

	SetLogFileName(opt(log));
	LogProgramInfoAndCmdLine();
	ResetRand(1);

	if (optset_seqlength)
		{
		extern uint BENCHL;
		BENCHL = opt(seqlength);
		}

	CMD Cmd = GetCmd();
	switch (Cmd)
		{
#define A(x)	case CMD_##x: { void cmd_##x(); opt(x); cmd_##x(); break; }
#include "cmds.h"
	default:
		asserta(false);
		}

	CheckUsedOpts(opt_log_used_opts);

	LogElapsedTimeAndRAM();
	return 0;
	}
