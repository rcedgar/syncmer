#ifndef randseq_h
#define randseq_h

extern uint BENCHL;

void MakeRandSeq(byte *Seq, uint L);
void MutateSeq(const byte *InputSeq, uint L, uint PctId, byte *MutatedSeq);

#endif // randseq_h
