#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstring>
#include <exception>
#include <htslib/sam.h>

using namespace std;

class Prune {
private:
	string bamfile;
	string table;
	string refSeq;
	unordered_map<string, long long> read_idx;
	unordered_map<int, unordered_map <int, vector<long long>>> pairdb;
	unordered_map<int, long> ctgdb;
	unordered_map<string, int> ctgidxdb;
	unordered_map<int, string> sctgdb;
	unordered_map<long long, long> removedb;

	bool Split(string source, string delim, vector<string>&target);
public:
	Prune();
	Prune(string bamfile, string table, string refSeq);
	~Prune();
	void SetParameter(string bamfile, string table, string refSeq);
	bool GeneratePairsAndCtgs();
	bool GenerateRemovedb();
	int CreatePrunedBam();
};
