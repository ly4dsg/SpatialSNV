import pysam 
from pysam import VariantFile
from multiprocess import Pool
import pandas as pd 
import os

from .utils import LoadChromFromBAM,logger_config,DirCheck

def fillter(read):
    m = 1
    if 'weak_evidence' in read.filter or 'germline' in read.filter or 'strand_bias' in read.filter or 'slippage' in read.filter or 'contamination' in read.filter:
        m = 0
    return m


def CallBack_process(bam,only_autosome,sample,out_dir,tmp_dir,vcf,barcode_tag,umi_tag,qmap,stereo,threads = 1):
    DirCheck(out_dir,create=True)
    DirCheck(tmp_dir,create=True)

    chrom_list = LoadChromFromBAM(bam,only_autosome)
    logger = logger_config(f'CallBack')

    if stereo:
        callback = StereoBack
        logger.debug(f'Call back on stereo-seq')
    else:
        callback = DropseqBack
        logger.debug(f'Call back on drop-seq like')

    logger.info(f'Save stat in {tmp_dir}')

    if threads == 1:
        for chrom in chrom_list:
            callback(bam,vcf,sample,tmp_dir,chrom,barcode_tag,umi_tag,qmap,logger)
    else:
        split_pool = Pool(threads)
        for chrom in chrom_list:
            split_pool.apply_async(callback, args=(bam,vcf,sample,tmp_dir,chrom,barcode_tag,umi_tag,qmap,logger))
        split_pool.close()
        split_pool.join()
    
    



def StereoBack(bam,vcf,sample,outdir,chrom,barcode_tag,umi_tag,qmap,logger):
    logger.info(f'Start {chrom} calling...')
    inbam = pysam.AlignmentFile(bam,"rb")
    vcf_input = VariantFile(vcf,"r")
    result_name = os.path.join(outdir,f'{sample}.{chrom}.stat')
    result = open(result_name,"w")
    total = 0
    for read in vcf_input.fetch(chrom,reopen=True):
        total += 1
        if all(len(allele) == 1 for allele in read.alleles) and fillter(read):
            vcf_pos = read.pos
            for i in inbam.pileup(chrom,vcf_pos-1,vcf_pos,truncate=True,min_mapping_quality=qmap):
                for pileupread in i.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_qcfail and not pileupread.alignment.is_duplicate:
                        read_base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if read_base in read.alts:
                            barcode = str(pileupread.alignment.get_tag('Cx'))+'_'+str(pileupread.alignment.get_tag('Cy'))
                            umi = str(pileupread.alignment.get_tag(umi_tag))
                            snv_name = "%s_%s:%s>%s"%(chrom,vcf_pos,read.ref,read_base)
                            result.write("%s,%s,%s\n"%(snv_name,barcode,umi))


    if total == 0:
        logger.warning(f'Not Found Any SNV in {chrom}')
    else:
        logger.debug(f'>> Finish {chrom} calling ::: Process {total} SNVs')

def DropseqBack(bam,vcf,sample,outdir,chrom,barcode_tag,umi_tag,qmap,logger):
    logger.info(f'Start {chrom} calling...')
    inbam = pysam.AlignmentFile(bam,"rb")
    vcf_input = VariantFile(vcf,"r")
    result_name = os.path.join(outdir,f'{sample}.{chrom}.stat')
    result = open(result_name,"w")

    total = 0
    for read in vcf_input.fetch(chrom,reopen=True):
        total += 1
        if all(len(allele) == 1 for allele in read.alleles) and fillter(read):
            vcf_pos = read.pos
            for i in inbam.pileup(chrom,vcf_pos-1,vcf_pos,truncate=True,min_mapping_quality=qmap):
                for pileupread in i.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_qcfail and not pileupread.alignment.is_duplicate and pileupread.alignment.has_tag(barcode_tag) and pileupread.alignment.has_tag(umi_tag):
                        read_base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if read_base in read.alts:
                            barcode = str(pileupread.alignment.get_tag(barcode_tag))
                            umi = str(pileupread.alignment.get_tag(umi_tag))
                            snv_name = "%s_%s:%s>%s"%(chrom,vcf_pos,read.ref,read_base)
                            result.write("%s,%s,%s\n"%(snv_name,barcode,umi))

    if total == 0:
        logger.warning(f'Not Found Any SNV in {chrom}')
    else:
        logger.debug(f'>> Finish {chrom} calling ::: Process {total} SNVs')