#!usr/bin/python
#########
####   Developed by xuchu 2020 ,used for analysis subcellular RNA 
#########

import argparse
import os
import configparser
import subprocess

def mapping_paired(input_dirs,output_dir):
    f = os.listdir(input_dirs)
    for i in f:
        if i.endswith('_1.fq.gz'):
            R1 = i
    for i in f:
        if i.endswith('_2.fq.gz'):
            R2 = i
    print(R1)
    os.chdir(output_dir)
    name_R1 = R1.split('.')[0]
    name_R2 = R2.split('.')[0]
    name = R1.split('_')[0]
    cmd_trim = 'nohup trim_galore --paired --quality 20 -o {0} {1}/{2} {1}/{3}'.format(output_dir,input_dirs,R1,R2)
    cmd_clumpify = 'nohup clumpify.sh qin=auto in={0}_val_1.fq.gz in2={1}_val_2.fq.gz out={0}_R1.dedup.fastq out2={1}_R2.dedup.fastq dedupe=t repair=f lowcomplexity=f rcomp=f deleteinput=f deletetemp=t'.format(name_R1,name_R2)
    cmd_starmapping = 'nohup /gpfs/hulab/apps/miniconda3/py310/envs/xuchu/bin/STAR --runThreadN 40 --quantMode TranscriptomeSAM --genomeDir /gpfs/hulab/lulu/data/xuchu/software/genomeindex/human/STAR_index/Normal_RNASeq/hg38 --readFilesIn {0}_R1.dedup.fastq {1}_R2.dedup.fastq --outFileNamePrefix ./{2} --outFilterMultimapNmax 20 --outSAMtype BAM Unsorted'.format(name_R1,name_R2,name)
    cmd_RSEM = 'nohup rsem-calculate-expression --paired-end -no-bam-output --alignments -p 20 {0}Aligned.toTranscriptome.out.bam /gpfs/hulab/lulu/data/xuchu/software/genomeindex/human/RSEM_index/Normal_RNASeq/hg38/hg38 ./RSEM_{0}'.format(name) 
    
   # subprocess.call(cmd_trim,shell=True)
    subprocess.call(cmd_clumpify,shell=True)
    subprocess.call(cmd_starmapping,shell=True)
    subprocess.call(cmd_RSEM,shell=True)
    #break
def main():
    parser = argparse.ArgumentParser(description='Help Manual ==> design for RNAseq'
									 ,epilog     ='design by XuChu')
    parser.add_argument("-T","-treat",required=True, nargs = '*', help="treat file dir, example: -T /d1 /d2 /d3")
    parser.add_argument("-C","-control", required=True,nargs = '*', help= "control file dir, example: -C /di1 /di2 /di3")
    parser.add_argument("-D","-dir", required=True, help = "output dir")
    args = parser.parse_args()
    print(args)
    output_dir = args.D
    R_treat = args.T
    R_control = args.C
    print(R_treat)
    print(R_control)
    #os.system('conda activate xuchu')
    samplenames = []##get all sample name
    condition = []
    for i in R_treat:
        r_out = mapping_paired(i,output_dir)
        samplenames.append(r_out)
        condition.append('Treat')
    for i in R_control:
        r_con = mapping_paired(i,output_dir)
        samplenames.append(r_con)
        condition.append('Control')
    title = samplenames[0]
    sams = [i+'Aligned.out.bam' for i in samplenames]
    sps = ' '.join(sams)
    condition = zip(sams,condition)
    condition = ['\t'.join(i) for i in condition]
    print(condition)
    os.system('featureCounts -T 5 -t exon -g gene_id -a /gpfs/hulab/lulu/reference/human_ref/hg38.refGene.gtf -p -o out_counts.txt {1}'.format(title,sps))
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

