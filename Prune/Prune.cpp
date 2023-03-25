#include "Prune.h"
#include <htslib/sam.h>


Prune::Prune() {
	this->bamfile = "";
	this->table = "";
}

Prune::Prune(std::string bamfile, std::string table) {
	this->bamfile = bamfile;
	this->table = table;
}

Prune::~Prune(){}

//Split string by delimiter
bool Prune::Split(std::string source, std::string delim, std::vector<std::string>&target) {
	target.clear();
	char *p;
	p = strtok(const_cast<char*>(source.c_str()), delim.c_str());
	if (!p) {
		return false;
	}
	while (p) {
		target.push_back(p);
		p = strtok(NULL, delim.c_str());
	}
	return true;
}


void Prune::SetParameter(std::string bamfile, std::string table) {
	this->bamfile = bamfile;
	this->table = table;
}


//Read bamfiles and read them by samtools, then create pairdbs and ctgdbs;
bool Prune::GeneratePairsAndCtgs() {
	if (bamfile == "" || table == "") {
		return false;
	}
	else {
		int ctg1, ctg2;
		std::string sctg1, sctg2;
		bam1_t *rec = bam_init1();
		htsFile *inbam = hts_open(bamfile.c_str(), "rb");
		sam_hdr_t *hdr = sam_hdr_read(inbam);
		int res;

		while((res = sam_read1(inbam, hdr, rec))>=0){
			ctg1 = rec->core.tid;
			ctg2 = rec->core.mtid;
			if(ctg2==-1){
				continue;
			}
			sctg1 = hdr->target_name[ctg1];
			sctg2 = hdr->target_name[ctg2];

			ctgidxdb[sctg1] = ctg1;
			ctgidxdb[sctg2] = ctg2;
			
			sctgdb[ctg1] = sctg1;
			sctgdb[ctg2] = sctg2;
			
			if(ctg1==ctg2){
				continue;
			}
			if(sctg1.compare(sctg2)>=0){
				int tmp = ctg1;
				ctg1 = ctg2;
				ctg2 = tmp;
			}
			pairdb[ctg1][ctg2]++;
			ctgdb[ctg1]++;
			ctgdb[ctg2]++;
		}
		hts_close(inbam);
		delete rec;
		delete hdr;
	}
	return true;
}

//Create removedb_Allele.txt, removedb_nonBest.txt and log.txt;
bool Prune::GenerateRemovedb() {
	std::ifstream fin;
	std::unordered_set<std::string>tempset;
	std::unordered_map<int, int> retaindb;
	std::unordered_map<int, int> numdb;
	std::unordered_map<int, std::unordered_set <int>> removedb;
	std::vector<std::string>data;
	std::string temp;
	std::string sctg1, sctg2;
	int ctg1, ctg2;
	long long num_r;
	
	fin.open(table);
	if (fin) {
		while (getline(fin, temp)) {
			tempset.clear();
			Split(temp, "\t", data);
			if (data.size() <= 3) {
				continue;
			}
			for (long i = 2; i < data.size() - 1; i++) {
				sctg1 = data[i];
				for (long j = i + 1; j < data.size(); j++) {
					sctg2 = data[j];
					ctg1 = ctgidxdb[sctg1];
					ctg2 = ctgidxdb[sctg2];
					if(sctg1.compare(sctg2)>=0){
						int tmp = ctg1;
						ctg1 = ctg2;
						ctg2 = tmp;
					}
					removedb[ctg1].insert(ctg2);
				}
			}
			retaindb.clear();
			numdb.clear();
			for (long i = 2; i < data.size(); i++) {
				sctg1 = data[i];
				ctg1 = ctgidxdb[sctg1];
				for(std::unordered_map<int, long>::iterator iter=ctgdb.begin(); iter!=ctgdb.end(); iter++){
					ctg2 = iter->first;
					sctg2 = sctgdb[ctg2];
					int nctg1=ctg1, nctg2=ctg2;
					if(sctg1.compare(sctg2)>=0){
						nctg1 = ctg2;
						nctg2 = ctg1;
					}
					if(removedb.count(nctg1) && removedb[nctg1].count(nctg2)){
						continue;
					}
					if(pairdb.count(nctg1)==0 || pairdb[nctg1].count(nctg2)==0){
						continue;
					}
					num_r = pairdb[nctg1][nctg2];
					if(retaindb.count(ctg2)==0){
						retaindb[ctg2] = ctg1;
						numdb[ctg2] = num_r;
					}else{
						if(num_r>numdb[ctg2]){
							retaindb[ctg2] = ctg1;
							numdb[ctg2] = num_r;
						}
					}
				}
			}
			for(std::unordered_map<int, int>::iterator iter=retaindb.begin(); iter!=retaindb.end(); iter++){
				ctg1 = iter->first;
				ctg2 = iter->second;
				sctg1 = sctgdb[ctg1];
				sctg2 = sctgdb[ctg2];
				if(sctg1.compare(sctg2)>=0){
					int tmp = ctg1;
					ctg1 = ctg2;
					ctg2 = tmp;
				}
				allretaindb[ctg1].insert(ctg2);
			}
		}
	}else{
		return false;
	}
	return true;
}

//Directly to create prunning.bam through pipe with samtools
long long Prune::CreatePrunedBam() {
	int ctg1, ctg2;
	std::string sctg1, sctg2;
	std::string outbam = "prunning.bam";
	htsFile *in = hts_open(bamfile.c_str(), "rb");
	htsFile *out = hts_open(outbam.c_str(), "wb");
	sam_hdr_t *hdr = sam_hdr_read(in);
	bam1_t *rec = bam_init1();
	int res;
	long long rmcnt = 0;

	if(sam_hdr_write(out, hdr)<0){
		hts_close(in);
		hts_close(out);
		delete hdr;
		delete rec;
		return -1;
	}
	while((res = sam_read1(in, hdr, rec))>=0){
		ctg1 = rec->core.tid;
		ctg2 = rec->core.mtid;
		sctg1 = sctgdb[ctg1];
		sctg2 = sctgdb[ctg2];
		if(sctg1.compare(sctg2)>=0){
			int tmp = ctg1;
			ctg1 = ctg2;
			ctg2 = tmp;
		}
		if(allretaindb.count(ctg1)==0 || allretaindb[ctg1].count(ctg2)==0 || rec->core.mtid==-1){
			rmcnt++;
			continue;
		}
		if(sam_write1(out, hdr, rec)<0){
			hts_close(in);
			hts_close(out);
			delete hdr;
			delete rec;
			return -1;
		}
	}
	if(hts_close(in)>=0&&hts_close(out)>=0){
		delete hdr;
		delete rec;
		return rmcnt;
	}

	return -1;
}
