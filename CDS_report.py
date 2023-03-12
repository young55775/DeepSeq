import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse


# show low cov range
## example: a = [1,2,3,4,5,6,10,11,12,13,15,16,19], return [(1, 6), (10, 13), (15, 16)]
def stat2region(lst):
    start = lst[0]
    ran = []
    for i in range(len(lst) - 1):
        end = lst[i]
        if lst[i + 1] != lst[i] + 1:
            end = lst[i]
            ran.append((start, end))
            start = lst[i + 1]
        else:
            continue
    ran.append((lst[-1], lst[-1]))
    return ran


def ran_by_chrom(df):
    lst = df['pos'].to_list()
    ran = stat2region(lst)
    chrom = [df['chr'].to_list()[0]] * len(ran)
    return pd.DataFrame({'chr': chrom, 'range': ran})


def recognition(pos, anno):
    lst1 = anno[anno['pos1'] <= pos].index
    lst2 = anno[anno['pos2'] >= pos].index
    return [n for n in lst2 if n in lst1]


# tup1 exon range; tup2 seq range
def intersect(tup1, tup2):
    a1 = tup1[0]
    a2 = tup1[1]
    b1 = tup2[0]
    b2 = tup2[1]
    # judge if there is an intersection
    if a2 < b1 or a1 > b2:
        return False
    else:
        lst = [a1, a2, b1, b2]
        lst.sort()
    return lst[1], lst[2]


def search(data, anno, chromosome):
    pos = list(zip(anno['pos1'].to_list(), anno['pos2'].to_list()))
    exon = anno['exon'].to_list()
    gene = anno['gene'].to_list()
    # start
    out = []
    for i in range(len(pos)):
        if intersect(pos[i], data):
            a = intersect(pos[i], data)
            out.append("{}\t{}\texon{}\t({},{})\tlength={}\n".format(chromosome, gene[i], exon[i], a[0], a[1],
                                                                     (a[1] - a[0]) + 1))
    return out


def main(ref, cov, outpath, cutoff):
    # read file
    depth = pd.read_csv(str(cov), sep='\t', names=['chr', 'pos', 'cov'])
    info = pd.read_csv(str(ref), sep='\t', names=['gene', 'exon', 'chr', 'pos1', 'pos2'])
    depth_filter = depth.drop(depth[depth['cov'] > cutoff].index)
    del depth
    # del < 300
    # depth_clear = depth.drop(depth[depth['cov'] < 300].index)
    # print(depth_clear)
    # del depth
    # mean = depth_clear['cov'].mean()
    # std = depth_clear['cov'].std()
    # depth_filter = depth_clear.drop(depth_clear[depth_clear['cov'] > mean - std].index)
    # print(depth_filter)
    # del depth_clear
    if depth_filter.size != 0:
        # group each chromosome
        result = depth_filter.groupby('chr').apply(ran_by_chrom)
        m = list(zip(result['chr'].to_list(), result['range'].to_list()))
        dt = {}
        for i in m:
            if i[0] not in dt.keys():
                dt[i[0]] = []
            dt[i[0]].append(i[1])
        # annotation
        # chrom_dict = {'NC_003279.8': 'I', 'NC_003280.10': 'II', 'NC_003281.10': 'III', 'NC_003282.8': 'IV',
        #               'NC_003283.11': 'V',
        #               'NC_003284.9': 'X'}
        # chrom = {v: k for k, v in chrom_dict.items()}

        out = []
        for k, v in dt.items():
            anno = info[info['chr'] == k]
            for i in tqdm(v, desc='chromosome {}'.format(k)):
                out.append(search(i, anno, k))
        # output
        with open(str(outpath) + '/anno_on_CDS.txt', 'w+') as f:
            for i in out:
                if i:
                    for j in i:
                        f.write(j)
    else:
        print('Nice seq quality!!!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", help='exon_range file path')
    parser.add_argument("--cov", help="file of coverage per base")
    parser.add_argument("--out", help="output path")
    parser.add_argument("--cut", help="set cutoff, default 1000", default=1000)
    args = parser.parse_args()
    if not args.ref:
        raise Exception("No reference file!")
    elif not args.cov:
        raise Exception("No coverage file!")
    elif not args.out:
        raise Exception("Not output path")
    else:
        main(args.ref, args.cov, args.out, int(args.cut))
