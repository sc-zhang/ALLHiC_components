#pragma once

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <cstring>
#include <htslib/sam.h>

using namespace std;

class Prune {
private:
	string bamfile;
	string table;
	unordered_map<int, unordered_map <int, long long>> pairdb;
	unordered_map<int, long> ctgdb;
	unordered_map<string, int> ctgidxdb;
	unordered_map<int, string> sctgdb;
	unordered_map<int, unordered_set <int>> allretaindb;

	bool Split(string source, string delim, vector<string>&target);
public:
	Prune();
	Prune(string bamfile, string table);
	~Prune();
	void SetParameter(string bamfile, string table);
	bool GeneratePairsAndCtgs();
	bool GenerateRemovedb();
	long long CreatePrunedBam();
};
