import pysam 
from pysam import VariantFile
from multiprocess import Pool
import pandas as pd 
import os
import glob
import shutil
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
    alt_result_name = os.path.join(outdir,f'{sample}.{chrom}.alt.stat')
    ref_result_name = os.path.join(outdir,f'{sample}.{chrom}.ref.stat')
    alt_result = open(alt_result_name,"w")
    ref_result = open(ref_result_name,"w")
    total = 0
    for read in vcf_input.fetch(chrom,reopen=True):
        total += 1
        if all(len(allele) == 1 for allele in read.alleles) and fillter(read):
            vcf_pos = read.pos
            for i in inbam.pileup(chrom,vcf_pos-1,vcf_pos,truncate=True,min_mapping_quality=qmap):
                for pileupread in i.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_qcfail and not pileupread.alignment.is_duplicate:
                        read_base = pileupread.alignment.query_sequence[pileupread.query_position]
                        barcode = str(pileupread.alignment.get_tag('Cx'))+'_'+str(pileupread.alignment.get_tag('Cy'))
                        umi = str(pileupread.alignment.get_tag(umi_tag))
                        if read_base in read.alts:
                            snv_name = "%s_%s:%s>%s"%(chrom,vcf_pos,read.ref,read_base)
                            alt_result.write("%s,%s,%s\n"%(snv_name,barcode,umi))
                        if read_base == read.ref:
                            snv_name = "%s_%s:%s"%(chrom,vcf_pos,read.ref)
                            ref_result.write("%s,%s,%s\n"%(snv_name,barcode,umi))


    if total == 0:
        logger.warning(f'Not Found Any SNV in {chrom}')
    else:
        logger.debug(f'>> Finish {chrom} calling ::: Process {total} SNVs')

def DropseqBack(bam,vcf,sample,outdir,chrom,barcode_tag,umi_tag,qmap,logger):
    logger.info(f'Start {chrom} calling...')
    inbam = pysam.AlignmentFile(bam,"rb")
    vcf_input = VariantFile(vcf,"r")
    alt_result_name = os.path.join(outdir,f'{sample}.{chrom}.alt.stat')
    ref_result_name = os.path.join(outdir,f'{sample}.{chrom}.ref.stat')
    alt_result = open(alt_result_name,"w")
    ref_result = open(ref_result_name,"w")
    total = 0
    for read in vcf_input.fetch(chrom,reopen=True):
        total += 1
        if all(len(allele) == 1 for allele in read.alleles) and fillter(read):
            vcf_pos = read.pos
            for i in inbam.pileup(chrom,vcf_pos-1,vcf_pos,truncate=True,min_mapping_quality=qmap):
                for pileupread in i.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_qcfail and not pileupread.alignment.is_duplicate and pileupread.alignment.has_tag(barcode_tag) and pileupread.alignment.has_tag(umi_tag):
                        read_base = pileupread.alignment.query_sequence[pileupread.query_position]
                        barcode = str(pileupread.alignment.get_tag(barcode_tag))
                        umi = str(pileupread.alignment.get_tag(umi_tag))
                        if read_base in read.alts:
                            snv_name = "%s_%s:%s>%s"%(chrom,vcf_pos,read.ref,read_base)
                            alt_result.write("%s,%s,%s\n"%(snv_name,barcode,umi))
                        if read_base == read.ref:
                            snv_name = "%s_%s:%s"%(chrom,vcf_pos,read.ref)
                            ref_result.write("%s,%s,%s\n"%(snv_name,barcode,umi))

    if total == 0:
        logger.warning(f'Not Found Any SNV in {chrom}')
    else:
        logger.debug(f'>> Finish {chrom} calling ::: Process {total} SNVs')





def BulidResult(sample,sample_tmp_path,sample_matrix_path,x_offset,y_offset,bin_size,removetmp):
    # alt and supported reads
    alt_vcf_list = glob.glob(f"{sample_tmp_path}/*.alt.stat")
    df = pd.DataFrame(columns=['snv','barcode','umi'])
    for i in alt_vcf_list:
        if os.path.getsize(i):
            tmp=pd.read_csv(i,header=None)
            tmp.columns=['snv','barcode','umi']
            df = pd.concat([df,tmp])
    df = df.drop_duplicates(['snv','barcode','umi'])
    df['num'] = 1
    pi = pd.pivot_table(df,index=['barcode','snv'],values='num',aggfunc='count').reset_index()
    pi['num'] = 1
    deep = pd.pivot_table(df,index=['barcode','snv'],values='num',aggfunc='count').reset_index()# how many read suppose this snv
#     deep_name = os.path.join(out_dir,"SNV.depth")
#     stat_name = os.path.join(out_dir,"SNV.stat")
#     pi.to_csv(stat_name,index=None)
#     deep.to_csv(deep_name,index=None)
    if x_offset != None and bin_size != 1:
        pi['x'] = pi['barcode'].map(lambda x : (int(x.split('_')[0]) - x_offset)//bin_size*bin_size)
        pi['y'] = pi['barcode'].map(lambda x : (int(x.split('_')[1]) - y_offset)//bin_size*bin_size)
        pi['location'] = pi['x'].map(str)+'_'+pi['y'].map(str)
        pi = pd.pivot_table(pi,index=['location','snv'],values='num',aggfunc='sum')
        pi = pi.reset_index()
        pi = pi.rename({'location':'barcode'},axis=1)
    else:
        pi = pd.pivot_table(pi,index=['barcode','snv'],values='num',aggfunc='sum')
        pi = pi.reset_index()

    if x_offset != None and bin_size != 1:
        deep['x'] = deep['barcode'].map(lambda x : (int(x.split('_')[0]) - x_offset)//bin_size*bin_size)
        deep['y'] = deep['barcode'].map(lambda x : (int(x.split('_')[1]) - y_offset)//bin_size*bin_size)
        deep['location'] = deep['x'].map(str)+'_'+deep['y'].map(str)
        deep = pd.pivot_table(deep,index=['location','snv'],values='num',aggfunc='sum')
        deep = deep.reset_index()
        deep = deep.rename({'location':'barcode'},axis=1)
    else:
        deep = pd.pivot_table(pi,index=['barcode','snv'],values='num',aggfunc='sum')
        deep = deep.reset_index()
    
    rawsample_matrix_path = sample_matrix_path
    sample_matrix_path = os.path.join(rawsample_matrix_path,f'{sample}_alt')
    get_mtx(pi,sample_matrix_path)

    sample_matrix_path_depth = os.path.join(rawsample_matrix_path,f'{sample}_depth')
    get_mtx(deep,sample_matrix_path_depth)

    ref_vcf_list = glob.glob(f"{sample_tmp_path}/*.ref.stat")
    df = pd.DataFrame(columns=['snv','barcode','umi'])
    for i in ref_vcf_list:
        if os.path.getsize(i):
            tmp=pd.read_csv(i,header=None)
            tmp.columns=['snv','barcode','umi']
            df = pd.concat([df,tmp])
    df = df.drop_duplicates(['snv','barcode','umi'])
    df['num'] = 1
    pi = pd.pivot_table(df,index=['barcode','snv'],values='num',aggfunc='count').reset_index()
    pi['num'] = 1
    if x_offset != None and bin_size != 1:
        pi['x'] = pi['barcode'].map(lambda x : (int(x.split('_')[0]) - x_offset)//bin_size*bin_size)
        pi['y'] = pi['barcode'].map(lambda x : (int(x.split('_')[1]) - y_offset)//bin_size*bin_size)
        pi['location'] = pi['x'].map(str)+'_'+pi['y'].map(str)
        pi = pd.pivot_table(pi,index=['location','snv'],values='num',aggfunc='sum')
        pi = pi.reset_index()
        pi = pi.rename({'location':'barcode'},axis=1)
    else:
        pi = pd.pivot_table(pi,index=['barcode','snv'],values='num',aggfunc='sum')
        pi = pi.reset_index()
    
    sample_matrix_path_ref = os.path.join(rawsample_matrix_path,f'{sample}_ref')
    get_mtx(pi,sample_matrix_path_ref)


    if removetmp:
        for i in (alt_vcf_list + ref_vcf_list):
            os.remove(i)


def get_mtx(pi,outdir_mtx):
    if not os.path.exists(outdir_mtx):
        os.makedirs(outdir_mtx)
    barcodes = pi['barcode'].drop_duplicates(keep='first').reset_index()
    barcodes_dict = dict(zip(barcodes['barcode'], barcodes.index + 1))
    features = pi['snv'].drop_duplicates(keep='first').reset_index()
    features_dict = dict(zip(features['snv'], features.index + 1))
    pi['barcode'] = pi['barcode'].apply(lambda x: barcodes_dict[x])
    pi['snv'] = pi['snv'].apply(lambda x: features_dict[x])
    hd = pd.DataFrame([['%%MatrixMarket matrix coordinate integer general'],\
                                   ['%'], \
                                   [' '.join([str(len(features)), str(len(barcodes)), str(len(pi))])]])
    hd.to_csv(outdir_mtx+'/matrix.mtx.gz', compression='gzip', sep = '\t', index=False, header=False)
    pi[['snv', 'barcode', 'num']].to_csv(outdir_mtx+'/matrix.mtx.gz', compression='gzip', mode='a+',\
                                                                                                 sep=' ', header=False, index=False)
    barcodes['barcode'].to_csv(outdir_mtx+'/barcodes.tsv.gz', compression='gzip', sep='\t', header=False, index=False)
    features['snv_id'] = features['snv']
    features['feature_type'] = 'Gene Expression'
    features[['snv', 'snv_id', 'feature_type']].to_csv(outdir_mtx+'/features.tsv.gz', compression='gzip', sep='\t', header=False, index=False)