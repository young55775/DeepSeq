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


def search(data, anno, chromosome):
    beg = []
    end = []

    pos = list(zip(anno['pos1'].to_list(), anno['pos2'].to_list()))
    # start
    start = data[0]
    ind = recognition(start, anno)
    if ind:
        for i in ind:
            beg.append(chromosome + '\t' + str(start) + '\t' + str(anno.loc[i]['gene']) + '\t' + "exon{}".format(
                anno.loc[i]['exon']) + '\t' + str(data[0] - int(anno.loc[i]['pos1'])))
    else:
        beg.append(chromosome + '\t' + str(start) + '\t' + 'Non-coding region')
    # end
    fin = data[1]
    ind = recognition(fin, anno)
    if ind:
        for i in ind:
            end.append(chromosome + '\t' + str(fin) + '\t' + str(anno.loc[i]['gene']) + '\t' + "exon{}".format(
                anno.loc[i]['exon']) + '\t' + str(data[1] - int(anno.loc[i]['pos1'])))
    else:
        end.append(chromosome + '\t' + str(fin) + '\t' + 'Non-coding region')
    return beg, end


def main(ref, cov, outpath):
    # read file
    depth = pd.read_csv(str(cov), sep='\t', names=['chr', 'pos', 'cov'])
    info = pd.read_csv(str(ref), sep='\t', names=['gene', 'exon', 'chr', 'pos1', 'pos2'])
    # del < 300
    depth_clear = depth.drop(depth[depth['cov'] < 300].index)
    del depth
    mean = depth_clear['cov'].mean()
    std = depth_clear['cov'].std()
    depth_filter = depth_clear.drop(depth_clear[depth_clear['cov'] > mean - std].index)
    del depth_clear
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
    with open(str(outpath) + '/anno.txt', 'w+') as f:
        for i in out:
            for j in i[0]:
                f.write('{} | '.format(j))
            f.write('->')
            for j in i[1]:
                f.write('{} | '.format(j))
            f.write('\n')

    # output filter
    with open(str(outpath) + '/anno_filter.txt', 'w+') as f:
        for i in out:
            if i[0][0].split('\t')[-1] != 'Non-coding region' or i[1][0].split('\t')[-1] != 'Non-coding region':
                for j in i[0]:
                    f.write('{} | '.format(j))
                f.write('->')
                for j in i[1]:
                    f.write('{} | '.format(j))
                f.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", help='exon_range file path')
    parser.add_argument("--cov", help="file of coverage per base")
    parser.add_argument("--out", help="output path")
    args = parser.parse_args()
    if not args.ref:
        raise Exception("No reference file!")
    elif not args.cov:
        raise Exception("No coverage file!")
    elif not args.out:
        raise Exception("Not output path")
    else:
        main(args.ref, args.cov, args.out)
