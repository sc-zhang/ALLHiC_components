#!/usr/bin/env python
import sys
import os
import argparse
import gc
import time


def time_print(str):
    print("\033[32m%s\033[0m %s"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), str))


def get_opts():
	group = argparse.ArgumentParser()
	group.add_argument("-r", "--reference", help="Chromosome fasta", required=True)
	group.add_argument("-c", "--contig", help="Contig fasta", required=True)
	group.add_argument("-p", "--ploidy", help="Ploidy of genome", required=True, type=int)
	group.add_argument("-o", "--output", help="Output allele table", required=True)
	group.add_argument("-w", "--win", help="Size of window, defalut: 20k", default="20k")
	group.add_argument("-s", "--step", help="Size of step, default: 10k", default="10k")
	group.add_argument("-d", "--dir", help="Work folder, default: wrk_dir", default="wrk_dir")
	group.add_argument("-t", "--threads", help="Number of thread, default: 1", type=int, default=1)
	return group.parse_args()


def read_fasta(in_fa):
	fa_db = {}
	with open(in_fa, 'r') as fin:
		for line in fin:
			if line[0] == '>':
				id = line.strip().split()[0][1:]
				fa_db[id] = []
			else:
				fa_db[id].append(line.strip())
	
	for id in fa_db:
		fa_db[id] = ''.join(fa_db[id])
	return fa_db


def gen_sub_seq(fa_db, win_size, step_size, wrk_dir):
	with open(os.path.join(wrk_dir, 'sub_chr.fa'), 'w') as fout:
		for id in sorted(fa_db):
			if 'tig' in id or 'ctg' in id:
				continue
			for i in range(0, len(fa_db[id])-win_size+1, step_size):
				sub_id = "%s-%d"%(id, i+1)
				sub_seq = fa_db[id][i: i+win_size]
				fout.write(">%s\n%s\n"%(sub_id, sub_seq))
	return os.path.join(wrk_dir, 'sub_chr.fa')


def gen_allele_table(ref_fa, ctg_fa, allele_table, ploidy, win_size, step_size, wrk_dir, threads):
	if not os.path.exists(wrk_dir):
		os.mkdir(wrk_dir)
	time_print("Loading reference genome")
	fa_db = read_fasta(ref_fa)
	
	time_print("Generating sequences")
	sub_chr_fn = gen_sub_seq(fa_db, win_size, step_size, wrk_dir)

	del fa_db
	gc.collect()
		
	time_print("Mapping")
	paf_fn = os.path.join(wrk_dir, "mapping.paf")
	cmd = "minimap2 -k19 -w19 -t%s %s %s > %s"%(threads, sub_chr_fn, ctg_fa, paf_fn)
	os.system(cmd)
	
	time_print("Generating allele table")
	map_db = {}
	map_len_db = {}
	with open(paf_fn, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			ctg = data[0]
			tsp = int(data[2])
			tep = int(data[3])
			chrn = data[5]
			chrn_base = chrn.split('-')[0]
			map_len = abs(tsp-tep)+1
			if ctg not in map_db:
				map_db[ctg] = []
			if ctg not in map_len_db:
				map_len_db[ctg] = {}
			if chrn_base not in map_len_db[ctg]:
				map_len_db[ctg][chrn_base] = 0
			map_len_db[ctg][chrn_base] += map_len
			map_db[ctg].append([map_len, chrn])
	
	new_map_db = {}
	for ctg in map_db:
		best_match = ''
		max_len = 0
		for chrn_base in map_len_db[ctg]:
			if map_len_db[ctg][chrn_base] > max_len:
				best_match = chrn_base
				max_len = map_len_db[ctg][chrn_base]
		for map_len, chrn in map_db[ctg]:
			chrn_base = chrn.split('-')[0]
			if chrn_base == best_match:
				if chrn not in new_map_db:
					new_map_db[chrn] = []
				new_map_db[chrn].append([map_len, ctg])
	
	allele_db = {}
	tmp_list = []
	for chrn in new_map_db:
		allele_db[chrn] = []
		tmp_list = sorted(new_map_db[chrn], reverse=True)
		for i in range(0, ploidy):
			if i >= len(tmp_list):
				break
			allele_db[chrn].append(tmp_list[i][1])
	
	tmp_list = []
	for chrn in allele_db:
		id, idx = chrn.split('-')
		idx = int(idx)
		tmp_list.append([id, idx, sorted(list(set(allele_db[chrn])))])
	time_print("Generating success")

	time_print("Writing allele table")
	with open(allele_table, 'w') as fout:
		for id, idx, allele_list in sorted(tmp_list):
			fout.write("%s\t%d\t%s\n"%(id, idx, '\t'.join(allele_list)))
	time_print("Writing success")

	del allele_db, tmp_list
	gc.collect()
	
	time_print("Finished")


if __name__ == "__main__":
	opts = get_opts()
	ref_fa = opts.reference
	ctg_fa = opts.contig
	ploidy = opts.ploidy
	win_size = opts.win
	step_size = opts.step
	threads = opts.threads
	allele_table = opts.output
	wrk_dir = opts.dir
	win_size = int(win_size.lower().replace('m', '000000').replace('k', '000'))
	step_size = int(step_size.lower().replace('m', '000000').replace('k', '000'))

	gen_allele_table(ref_fa, ctg_fa, allele_table, ploidy, win_size, step_size, wrk_dir, threads)

