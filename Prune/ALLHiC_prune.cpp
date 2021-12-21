#include <iostream>
#include <string>
#include <cstring>
#include <ctime>
#include "Prune.h"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc != 7) {
		cout << "************************************************************************\n";
		cout << "    Usage: "<<argv[0]<<" -i Allele.ctg.table -b sorted.bam -r draft.asm.fasta\n";
		cout << "      -h : help and usage.\n";
		cout << "      -i : Allele.ctg.table\n";
		cout << "      -b : sorted.bam\n";
		cout << "      -r : draft.asm.fasta\n";
		cout << "************************************************************************\n";
	}
	else {
		string bamfile;
		string table;
		string refSeq;
		clock_t startt, endt;
		startt = clock();
		for (long i = 1; i < 7; i += 2) {
			if (strcmp(argv[i], "-i") == 0) {
				table = argv[i + 1];
				continue;
			}
			if (strcmp(argv[i], "-b") == 0) {
				bamfile = argv[i + 1];
				continue;
			}
			if (strcmp(argv[i], "-r") == 0) {
				refSeq = argv[i + 1];
				continue;
			}
		}
		Prune prune;
		prune.SetParameter(bamfile, table, refSeq);
		cout<<"Getting contig pairs"<<endl; 
		prune.GeneratePairsAndCtgs();
		cout<<"Generating remove reads"<<endl;
		prune.GenerateRemovedb();
		cout<<"Creating prunned bam file"<<endl;
		int rmcnt = 0;
		rmcnt = prune.CreatePrunedBam();
		cout<<"Removed "<<rmcnt<<" reads"<<endl;
		
		endt = clock();
		cout << "use time: " << (endt - startt) / CLOCKS_PER_SEC << "s\n";
	}
	return 0;
}
