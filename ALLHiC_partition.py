#!/usr/bin/env python
import argparse
import os

import pysam
import numpy as np
import time


class UnionFind():
	def __init__(self, size):
		self.__f = [i for i in range(0, size)]
	
	
	def find(self, x):
		if(self.__f[x] == x):
			return x
		self.__f[x] = self.find(self.__f[x])
		return self.__f[x]
	

	def union(self, x, y):
		fx = self.find(x)
		fy = self.find(y)
		if fx != fy:
			self.__f[fy] = fx


def getOpts():
	groups = argparse.ArgumentParser()
	groups.add_argument('-r', '--ref', help="Contig level assembly fasta", required=True)
	groups.add_argument('-b', '--bam', help="Prunned bam file", required=True)
	groups.add_argument('-d', '--bed', help='dup.bed', required=True)
	groups.add_argument('-a', '--anchors', help='anchors file with dup.mono.anchors', required=True)
	groups.add_argument('-p', '--poly', help="Ploid count of polyploid", type=int, required=True)
	groups.add_argument('-e', '--exclude', help="A list file contains exclude contigs for partition, default=\"\"", default="")
	groups.add_argument('-o', '--out', help="Output directory, default=workdir", default="workdir")
	return groups.parse_args()


def getSignal(inBam, seqCount, seqList, qryDB, excludeDB):
	seqIdx = {}
	for i in range(0, seqCount):
		seqIdx[seqList[i]] = i

	seqMat = [[0 for i in range(0, seqCount)] for j in range(0, seqCount)]
	with pysam.AlignmentFile(inBam, 'rb') as fin:
		for line in fin:
			ctg1 = line.reference_name
			ctg2 = line.next_reference_name
			pos1 = line.reference_start+1
			pos2 = line.next_reference_start+1
			if pos1 == -1 or pos2 == -1 or ctg1 == ctg2 or ctg1 in excludeDB or ctg2 in excludeDB:
				continue
			idx1 = seqIdx[ctg1]
			idx2 = seqIdx[ctg2]
			if idx1>idx2:
				idx1, idx2 = idx2, idx1
			seqMat[idx1][idx2] += 1
	
	sigList = []
	for idx1 in range(0, seqCount-1):
		for idx2 in range(idx1+1, seqCount):
			if seqMat[idx1][idx2] >= 10:
				ctg1 = seqList[idx1]
				ctg2 = seqList[idx2]
				if (ctg1 not in qryDB) or (ctg2 not in qryDB) or (len(qryDB[ctg1])+len(qryDB[ctg2]))==0:
					ovlp = 0.0
				else:
					ovlpCount = len(qryDB[ctg1].intersection(qryDB[ctg2]))
					ovlp = ovlpCount*2.0/(len(qryDB[ctg1])+len(qryDB[ctg2]))
				sigList.append([idx1, idx2, seqMat[idx1][idx2], ovlp])
	return sigList


def checkLongestGroups(lengthList, polyCount):
	#avgL = np.average(lengthList[: polyCount])
	minL = min(lengthList[: polyCount])
	maxL = max(lengthList[: polyCount])
	print("\tMax %d groups, %s"%(polyCount, ','.join(map(str, lengthList[: polyCount]))))
	if maxL<=minL*3: #avgL*1.5>maxL and avgL*0.5<minL:
		return False
	else:
		return True


def allHiCPartition(refFasta, inBam, bed, anchors, polyCount, exclude, outDir):
	# Get full file path
	refFasta = os.path.abspath(refFasta)
	inBam = os.path.abspath(inBam)
	bed = os.path.abspath(bed)
	anchors = os.path.abspath(anchors)

	if exclude != "":
		exclude = os.path.abspath(exclude)

	outDir = os.path.abspath(outDir)
	if not os.path.exists(outDir):
		os.mkdir(outDir)
	
	# Enter work directory
	os.chdir(outDir)
	print("Loading fasta")
	excludeDB = {}
	if exclude != "":
		with open(exclude, 'r') as fin:
			for line in fin:
				excludeDB[line.strip()] = 1
	faDB = {}
	with open(refFasta, 'r') as fin:
		id = ""
		seq = ""
		for line in fin:
			if line[0] == '>':
				if seq != "" and id not in excludeDB:
					faDB[id] = seq
				id = line.strip()[1:]
				seq = ""
			else:
				seq += line.strip()
	if id not in excludeDB:
		faDB[id] = seq

	# Get overlap
	print("Loading anchors")
	anchorsDB = {}
	with open(anchors, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == '#':
				continue
			data = line.strip().split()
			anchorsDB[data[0]] = data[1]

	qryDB = {}
	with open(bed, 'r') as fin:
		for line in fin:
			data = line.strip().split()
			tig = data[0]
			gene = data[3]
			if tig not in qryDB:
				qryDB[tig] = set()
			if gene not in anchorsDB:
				continue
			qryDB[tig].add(anchorsDB[gene])

	# Get signals
	print("Getting signals")
	seqCount  = len(faDB)
	seqList = sorted(faDB)
	seqLen = []
	for i in range(0, seqCount):
		seqLen.append(len(faDB[seqList[i]]))

	sigList = getSignal(inBam, seqCount, seqList, qryDB, excludeDB)
	
	# Save signal list
	print("Saving signal list")
	with open("signal.txt", 'w') as fout:
		for idx1, idx2, signal, ovlp in sigList:
			fout.write("%s\t%s\t%d\t%f\n"%(seqList[idx1], seqList[idx2], signal, ovlp))
	
	# Initial UnionFind
	sigList = sorted(sigList, key=lambda x: (-x[3], x[2]))
	sigCount = len(sigList)
	print("Generating Union find")
	uf = UnionFind(seqCount)
	for idx1, idx2, signal, ovlp in sigList:
		uf.union(idx1, idx2)
	
	# Get current groups
	currentGroupCount = 0
	for idx in range(0, seqCount):
		if uf.find(idx) == idx:
			currentGroupCount += 1
	
	print("\tInitial group count: %d, edge count: %d"%(currentGroupCount, sigCount))

	groupDB = {}
	for idx in range(0, seqCount):
		gid = uf.find(idx)
		if gid not in groupDB:
			groupDB[gid] = []
		groupDB[gid].append(idx)	
	
	lengthList = []
	for gid in groupDB:
		curLen = 0
		for idx in groupDB[gid]:
			curLen += seqLen[idx]
		lengthList.append(curLen)
	lengthList = sorted(lengthList, reverse=True)

	i = 1
	while sigList[i][3] > 0:
		i += 1

	print("\tRemoved: %d edges while contigs were overlaped"%i)

	while currentGroupCount < polyCount or checkLongestGroups(lengthList, polyCount):
		sig = sigList[i][2]
		while sigList[i][2] == sig:
			i += 1
		
		uf = UnionFind(seqCount)
		for idx in range(i, sigCount):
			idx1, idx2, signal, ovlp = sigList[idx]
			uf.union(idx1, idx2)
		
		lengthList = []
		currentGroupCount = 0
		for idx in range(0, seqCount):
			if uf.find(idx) == idx:
				currentGroupCount += 1

		groupDB = {}
		for idx in range(0, seqCount):
			gid = uf.find(idx)
			if gid not in groupDB:
				groupDB[gid] = []
			groupDB[gid].append(idx)	
		
		lengthList = []
		for gid in groupDB:
			curLen = 0
			for idx in groupDB[gid]:
				curLen += seqLen[idx]
			lengthList.append(curLen)
		lengthList = sorted(lengthList, reverse=True)

		print("\tCurrent group count: %d, removed edge count: %d"%(currentGroupCount, i))
		i += 1


	with open("remove.list", "w") as fout:
		for idx in range(0, i):
			idx1, idx2, signal, ovlp = sigList[idx]
			fout.write("Remove %d: %s, %s, %d, %f\n"%(idx+1, seqList[idx1], seqList[idx2], signal, ovlp))
	

	checkLongestGroups(lengthList, polyCount)
	lengthDB = {}
	for gid in groupDB:
		curLen = 0
		for idx in groupDB[gid]:
			curLen += seqLen[idx]
		lengthDB[gid] = curLen

	groupList = []
	for gid in groupDB:
		groupList.append([gid, lengthDB[gid]])

	groupList = sorted(groupList, key=lambda x: -x[1])

	print("Writing group list")
	with open("group.txt", "w") as fout:
		for i in range(0, len(groupList)):
			idx = groupList[i][0]
			fout.write("group%d\t"%(i+1))
			tmp = []
			for subIdx in sorted(groupDB[idx]):
				tmp.append(seqList[subIdx])
			fout.write("%s\n"%'\t'.join(tmp))

	print("Finished")


if __name__ == "__main__":
	opts = getOpts()
	refFasta = opts.ref
	inBam = opts.bam
	bed = opts.bed
	anchors = opts.anchors
	polyCount = opts.poly
	exclude = opts.exclude
	outDir = opts.out
	allHiCPartition(refFasta, inBam, bed, anchors, polyCount, exclude, outDir)
