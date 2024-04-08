#!usr/bin/python
########
####   Developed by xuchu ,used for analysis subcellular RNA 
#########

import argparse
import os
import configparser

def mapping(fastq):
    name = fastq.split('/')[-1].split('_')[0]
    os.system("nohup trim_galore {0} --fastqc --quality 20 --phred33 --length 15".format(fastq))
    os.system("nohup clumpify.sh qin=auto in={0}_R1_trimmed.fq.gz out={0}_dedup.fq dedupe=t repair=f lowcomplexity=f rcomp=f deleteinput=f deletetemp=t ".format(name))
    os.system('kallisto quant -i /gpfs/hulab/lulu/data/xuchu/software/genomeindex/kallisto/Human.hg38.gencodev44.index -o ./{0} --single -l 30 -s 10 -t 20 -g /gpfs/hulab/lulu/data/xuchu/software/genomeindex/human/gencode.v44.annotation.gtf ./{0}_dedup.fq'.format(name))
    os.system('kallisto quant -i /gpfs/hulab/lulu/data/xuchu/software/genomeindex/kallisto/tRNA/hg38_tRNA -o ./{0} --single -l 50  -s 20 -t 20 ./{0}_dedup.fq'.format(name))
    return name
	
def main():
    parser = argparse.ArgumentParser(description='Help Manual ==> design for RNAseq'
									 ,epilog     ='design by XuChu')
    parser.add_argument("-T","-treat",required=True, nargs = '*', help="treat file dir, example: -T /d1 /d2 /d3")
    parser.add_argument("-C","-control", nargs = '*', help= "control file dir, example: -C /di1 /di2 /di3")
    args = parser.parse_args()
    print(args)
    R_treat = args.T
    R_control = args.C
    print(R_treat)
    print(R_control)
    samplenames = []##get all sample name
    condition = []
    for i in R_treat:
        r_out = mapping(i)
        samplenames.append(r_out)
        condition.append('Treat')
    for i in R_control:
        r_con = mapping(i)
        samplenames.append(r_con)
        condition.append('Control')
    title = samplenames[0]
    sams = [i+'Aligned.sortedByCoord.out.bam' for i in samplenames]
    sps = ' '.join(sams)
    condition = zip(sams,condition)
    condition = ['\t'.join(i) for i in condition]
    print(condition)
    os.system('featureCounts -T 5 -t exon -g gene_id -a /gpfs/hulab/lulu/data/xuchu/software/genomeindex/gencode.v44.annotation.gtf -o out_counts.txt {1}'.format(title,sps))
#/gpfs/hulab/lulu/reference/human_ref/hg38.refGene.gtf -o out_counts.txt {1}'.format(title,sps))
    with open('out_counts.txt') as handle:
        out = []
        handle.readline()
        for line in handle:
            line = line.strip().split('\t')
            name = line[0]
            value = line[6:]
            mat = '\t'.join([name]+value)
            out.append(mat)
    with open('rowcount_forDESeq.txt'.format(title),'w') as b:
        b.write('\n'.join(out))
    with open('sample.sheet'.format(title),'w') as c:
        c.write('\n'.join(condition))
if __name__ == '__main__':
	main()

