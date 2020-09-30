#ifndef cigar_h
#define cigar_h

void PathToCIGAR(const char *Path, string &CIGAR);
unsigned CIGARToQL(const string &CIGAR);
void CIGARGetOps(const string &CIGAR, vector<char> &Ops,
  vector<uint> &Lengths);

#endif // cigar_h
