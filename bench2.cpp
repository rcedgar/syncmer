#include "myutils.h"
#include "syncmerindex.h"
#include "randseq.h"

void InitBench2(uint WindowSize);
void Bench2(uint k, uint t, uint w, const string &Name1, const string &Name2);

void cmd_bench2()
	{
	const string Names = opt(bench2);
	vector<string> Fields;
	Split(Names, Fields, '+');
	asserta(SIZE(Fields) == 2);
	asserta(optset_k && optset_t);
	uint k = opt(k);
	uint t = opt(t);
	uint w = (optset_w ? opt(w) : SyncmerIndex::EstimateWindow(k, t));
	uint WindowSize = 64;
	if (optset_window)
		WindowSize = opt(window);
	Log("k=%u, t=%u, w=%u, Window %u\n", k, t, w, WindowSize);
	const string &Name1 = Fields[0];
	const string &Name2 = Fields[1];
	InitBench2(WindowSize);
	Bench2(k, t, w, Name1, Name2);
	}
