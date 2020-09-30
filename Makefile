CPP = ccache g++
CPPOPTS = -fopenmp -msse -mfpmath=sse -O3 -DNDEBUG -c

CC = ccache gcc
CCOPTS = -fopenmp -msse -mfpmath=sse -O3 -DNDEBUG -c

LNK = g++
LNKOPTS = -O3 -fopenmp -pthread -lpthread -static

HDRS = \
  alpha.h \
  bestmatch1.h \
  cigar.h \
  cmds.h \
  countsort.h \
  fastaseqsource.h \
  fastq.h \
  fastqseqsource.h \
  fileseqsource.h \
  filetype.h \
  gobuff.h \
  jenkinshash.h \
  kmer.h \
  linereader.h \
  lockobj.h \
  lockobjs.h \
  murmur.h \
  myopts.h \
  myutils.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  omplock.h \
  primes.h \
  quarts.h \
  randseq.h \
  searchsix.h \
  seqdb.h \
  seqinfo.h \
  seqsource.h \
  sort.h \
  speciesindex.h \
  spindex.h \
  stypes.h \
  syncmerindex.h \
  syncmerindex2.h \

OBJS = \
  o/alpha.o \
  o/bench.o \
  o/bench1.o \
  o/bench2.o \
  o/cigar.o \
  o/compress.o \
  o/cov1fract.o \
  o/depthdist.o \
  o/dict.o \
  o/equiv.o \
  o/fastaseqsource.o \
  o/fastq.o \
  o/fastqseqsource.o \
  o/fileseqsource.o \
  o/filetype.o \
  o/kmer.o \
  o/linereader.o \
  o/lockobj.o \
  o/makespindex.o \
  o/mum2breaks.o \
  o/newbench.o \
  o/objmgr.o \
  o/paf2aln.o \
  o/paf2features.o \
  o/pair.o \
  o/prime.o \
  o/quarts.o \
  o/randseq.o \
  o/searchsix.o \
  o/seqdbfromfasta.o \
  o/hashtable.o \
  o/spacing.o \
  o/spacings.o \
  o/speciesalign.o \
  o/speciesindex.o \
  o/speciessearch.o \
  o/spindex.o \
  o/syncmerindex.o \
  o/syncmerindex2.o \
  o/syncmer_main.o \
  o/myutils.o \
  o/seqdb.o \
  o/seqinfo.o \
  o/seqsource.o \
  o/test.o \
  o/testdist.o \
  o/testminimizers.o \
  o/testsubmers.o \
  o/wga.o \

syncmer : o/ $(OBJS)
	$(LNK) $(LNKOPTS) $(OBJS) -o o/syncmer
	strip -d o/syncmer

o/ :
	mkdir -p o/

o/alpha.o : alpha.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/alpha.o alpha.cpp

o/bench.o : bench.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/bench.o bench.cpp

o/bench1.o : bench1.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/bench1.o bench1.cpp

o/bench2.o : bench2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/bench2.o bench2.cpp

o/cigar.o : cigar.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/cigar.o cigar.cpp

o/compress.o : compress.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/compress.o compress.cpp

o/cov1fract.o : cov1fract.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/cov1fract.o cov1fract.cpp

o/depthdist.o : depthdist.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/depthdist.o depthdist.cpp

o/dict.o : dict.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/dict.o dict.cpp

o/equiv.o : equiv.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/equiv.o equiv.cpp

o/fastaseqsource.o : fastaseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fastaseqsource.o fastaseqsource.cpp

o/fastq.o : fastq.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fastq.o fastq.cpp

o/fastqseqsource.o : fastqseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fastqseqsource.o fastqseqsource.cpp

o/fileseqsource.o : fileseqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/fileseqsource.o fileseqsource.cpp

o/filetype.o : filetype.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/filetype.o filetype.cpp

o/kmer.o : kmer.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/kmer.o kmer.cpp

o/linereader.o : linereader.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/linereader.o linereader.cpp

o/lockobj.o : lockobj.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/lockobj.o lockobj.cpp

o/makespindex.o : makespindex.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/makespindex.o makespindex.cpp

o/mum2breaks.o : mum2breaks.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/mum2breaks.o mum2breaks.cpp

o/newbench.o : newbench.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/newbench.o newbench.cpp

o/objmgr.o : objmgr.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/objmgr.o objmgr.cpp

o/paf2aln.o : paf2aln.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/paf2aln.o paf2aln.cpp

o/paf2features.o : paf2features.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/paf2features.o paf2features.cpp

o/pair.o : pair.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/pair.o pair.cpp

o/prime.o : prime.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/prime.o prime.cpp

o/quarts.o : quarts.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/quarts.o quarts.cpp

o/randseq.o : randseq.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/randseq.o randseq.cpp

o/searchsix.o : searchsix.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/searchsix.o searchsix.cpp

o/seqdbfromfasta.o : seqdbfromfasta.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqdbfromfasta.o seqdbfromfasta.cpp

o/hashtable.o : hashtable.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/hashtable.o hashtable.cpp

o/spacing.o : spacing.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/spacing.o spacing.cpp

o/spacings.o : spacings.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/spacings.o spacings.cpp

o/speciesalign.o : speciesalign.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/speciesalign.o speciesalign.cpp

o/speciesindex.o : speciesindex.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/speciesindex.o speciesindex.cpp

o/speciessearch.o : speciessearch.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/speciessearch.o speciessearch.cpp

o/spindex.o : spindex.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/spindex.o spindex.cpp

o/syncmerindex.o : syncmerindex.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/syncmerindex.o syncmerindex.cpp

o/syncmerindex2.o : syncmerindex2.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/syncmerindex2.o syncmerindex2.cpp

o/syncmer_main.o : syncmer_main.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/syncmer_main.o syncmer_main.cpp

o/myutils.o : myutils.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/myutils.o myutils.cpp

o/seqdb.o : seqdb.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqdb.o seqdb.cpp

o/seqinfo.o : seqinfo.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqinfo.o seqinfo.cpp

o/seqsource.o : seqsource.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/seqsource.o seqsource.cpp

o/test.o : test.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/test.o test.cpp

o/testdist.o : testdist.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/testdist.o testdist.cpp

o/testminimizers.o : testminimizers.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/testminimizers.o testminimizers.cpp

o/testsubmers.o : testsubmers.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/testsubmers.o testsubmers.cpp

o/wga.o : wga.cpp $(HDRS)
	$(CPP) $(CPPOPTS) -o o/wga.o wga.cpp
