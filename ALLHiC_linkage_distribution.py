#!/usr/bin/env python
import sys
import os
import pysam
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_linkage_dist(in_bam, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    print("Getting linkages between contigs")
    link_db = {}
    with pysam.AlignmentFile(in_bam, 'rb') as fin:
        for line in fin:
            for line in fin:
                ctg1 = line.reference_name
                ctg2 = line.next_reference_name
                pos1 = line.reference_start+1
                pos2 = line.next_reference_start+1
                if pos1 == -1 or pos2 == -1 or ctg1 == ctg2:
                    continue
                if ctg1 not in link_db:
                    link_db[ctg1] = {}
                if ctg2 not in link_db[ctg1]:
                    link_db[ctg1][ctg2] = 0
                if ctg2 not in link_db:
                    link_db[ctg2] = {}
                if ctg1 not in link_db[ctg2]:
                    link_db[ctg2][ctg1] = 0
                link_db[ctg1][ctg2] += 1
                link_db[ctg2][ctg1] += 1

    
    print("Writing linkage distribution")
    link_list = []
    for ctg in link_db:
        sig = 0
        #for ctg2 in link_db[ctg]:
        #    sig += link_db[ctg][ctg2]
        sig = len(link_db[ctg])
        link_list.append([ctg, sig])
    
    link_list = sorted(link_list, key=lambda x: -x[1])
    dist_db = {}
    with open(os.path.join(out_dir, 'linkages.txt'), 'w') as fout:
        for ctg, links in link_list:
            fout.write("%s\t%d\n"%(ctg, links))
            links = int(links/10)
            if links not in dist_db:
                dist_db[links] = 0
            dist_db[links] += 1
    
    x_vals = []
    y_vals = []
    for links in sorted(dist_db):
        x_vals.append(links)
        y_vals.append(dist_db[links])
    

    print("Drawing distributions")

    plt.figure(figsize=(10, 8), dpi=100)
    plt.plot(x_vals, y_vals)
    plt.savefig(os.path.join(out_dir, "dist.pdf"), bbox_inches="tight")

    print("Finished")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python %s <in_bam> <out_dir>"%sys.argv[0])
    else:
        in_bam, out_dir = sys.argv[1:]
        get_linkage_dist(in_bam, out_dir)
