import pysam
import os


from multiprocess import Pool

from .utils import LoadChromFromBAM,logger_config,DirCheck


def SplitChromBAM_process(bam,only_autosome,sampleid,out_dir,threads = 1):#chrom,bam,output_dir,sample_name):
    DirCheck(out_dir,create=True)
    chrom_list = LoadChromFromBAM(bam,only_autosome)
    logger = logger_config(f'SplitChromBAM')
    if threads == 1:
        for chrom in chrom_list:
            split_bam(bam,chrom,sampleid,out_dir,logger)
    else:
        split_pool = Pool(threads)
        for chrom in chrom_list:
            split_pool.apply_async(split_bam, args=(bam,chrom,sampleid,out_dir,logger))
        split_pool.close()
        split_pool.join()


def split_bam(bam,chrom,sample_name,output_dir,logger):
    inputfile = pysam.AlignmentFile(bam, "rb")
    tmp_bam = os.path.join(output_dir, '%s.%s.bam'%(sample_name,chrom))
    output = pysam.AlignmentFile(tmp_bam, 'wb', template = inputfile)
    m = 0
    for read in inputfile.fetch(str(chrom), multiple_iterators=True):
        output.write(read)
        m+=1
        if m>1000000:
            break
    logger.info(f'{chrom} finish')
    inputfile.close()
    output.close()
