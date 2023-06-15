#pragma once
#ifndef __PRUNE_H__
#define __PRUNE_H__
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <cstring>

class Prune {
private:
	std::string bamfile;
	std::string table;
	std::unordered_map<int, std::unordered_map <int, long long>> pairdb;
	std::unordered_map<int, long> ctgdb;
	std::unordered_map<std::string, int> ctgidxdb;
	std::unordered_map<int, std::string> sctgdb;
	std::unordered_map<int, std::unordered_set <int>> allremovedb;

	bool Split(std::string source, std::string delim, std::vector<std::string>&target);
public:
	Prune();
	Prune(std::string bamfile, std::string table);
	~Prune();
	void SetParameter(std::string bamfile, std::string table);
	bool GeneratePairsAndCtgs();
	bool GenerateRemovedb();
	long long CreatePrunedBam();
};

#endif