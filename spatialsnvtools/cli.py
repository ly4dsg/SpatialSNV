#!/usr/bin/env python3
import click
from .SplitChromBAM import SplitChromBAM_process
from .PerpareBAMforCalling import PerpareBAMforCalling_process
from .SNVCalling import SNVCalling_process
from .CallBack import CallBack_process,BulidResult
from .utils import BulidConfig


@click.group()
def main():
    pass

@main.command(name="SplitChromBAM",no_args_is_help=True)
@click.option("-b","--bam", type = click.Path(exists=True, file_okay=True, dir_okay=False), required = True, help = "")
@click.option("-s","--sample", type = str, required = True, help = "Sample ID ")
@click.option("-o","--outdir", type = str, required = True,  help = "Output dir")
@click.option("-@","--threads", type = int, default = 10, help = "Sets the number of threads")
@click.option('--only_autosome', is_flag=True,default=False,help="Only analysis autosome")
def SplitChromBAM(
    bam: str,
    sample: str,
    outdir: str,
    threads: int = 10,
    only_autosome: bool = False,
):
    SplitChromBAM_process(bam,only_autosome,sample,outdir,threads = threads)

@main.command(name="PerpareBAMforCalling",no_args_is_help=True)
@click.option("-b","--bam", type = click.Path(exists=True, file_okay=True, dir_okay=False), required = True, help = "")
@click.option("-o","--outdir", type = str, required = True,  help = "Output dir")
@click.option("-s","--proxy", type = str, required = True, help = "Sample ID ")
@click.option("-f","--fasta", type = click.Path(exists=True, file_okay=True, dir_okay=False), required = True, help = "Refence fasta")
@click.option("-d","--dbsnp", type = click.Path(exists=True, file_okay=True, dir_okay=False), required = True, help = "dbsnp")
@click.option("-c","--barcode", type = str, default = 'CR', help = "Cell Barcode in Bam file (ex. CR for 10X Genomics)")
@click.option("-u","--umi", type = str, default ='UR', help = "UMI(Molecular Barcodes) in Bam file (ex. UR for 10X Genomics)")
@click.option('--stereo', is_flag=True,default=False, help="Please check your data is stereo(barcode is Cx and Cy) or not")
@click.option('--gem',type = str, help="gem match your raw stereo bam")
@click.option("-x",'--xsetoff',type = int, help="gem_x + x_offset = bam_x")
@click.option("-y",'--ysetoff',type = int, help="gem_y + y_offset = bam_y")
@click.option('--tmpdir',type=str,help="Specify a temporary directory")
@click.option("-@","--threads", type = int, default = 10, help = "Sets the number of threads")
@click.option('--samtools',type=str,help="Specify samtools path")
@click.option('--picard',type=str,help="Specify picard.jar path")
@click.option('--gatk',type=str,help="Specify GATK path")
@click.option('--removetmp', is_flag=True,default=False, help="Remove all tmp file")
def PerpareBAMforCalling(
    bam: str,
    outdir: str,
    proxy: str,
    fasta: str,
    dbsnp: str,
    barcode: str = 'CB',
    umi: str = 'UB',
    stereo: bool = True,
    gem: str = None,
    xsetoff: int = 0,
    ysetoff: int = 0,
    tmpdir: str = None,
    samtools: str = None,
    gatk: str = None,
    picard: str = None,
    threads: int = 10,
    removetmp: bool = False,
    ):
    config = BulidConfig(samtools = samtools, gatk = gatk, picard = picard, threads= threads)
    PerpareBAMforCalling_process(bam,outdir,proxy,barcode,umi,stereo,gem,xsetoff,ysetoff,config,fasta,dbsnp,tmpdir,removetmp)


@main.command(name="SNVCalling",no_args_is_help=True)
@click.option("-s","--sample", multiple=True,required=True, help="Bam file with inde(ex: -b example.bam1 -b example.bam2)")
@click.option("-b","--bam", type=click.Path(exists=True, file_okay=True, dir_okay=False), multiple=True,required=True, help="Bam file with inde(ex: -b example.bam1 -b example.bam2)")
@click.option("-o","--outvcf", type = str, required = True,  help = "Output VCF")
@click.option("-f","--fasta", type = click.Path(exists=True, file_okay=True, dir_okay=False), required = True, help = "Refence fasta")
@click.option("--pon",type=click.Path(exists=True, file_okay=True, dir_okay=False),required = True,help="PON")
@click.option("--germline",type=click.Path(exists=True, file_okay=True, dir_okay=False),required = True, help="Germline source")
@click.option('--gatk',type=str,help="Specify GATK path")
@click.option("-L","--chrom",type = str, required = False,  help = "Specific CHROM")
def SNVCalling(
    sample: str,
    bam: str,
    fasta: str,
    outvcf: str,
    pon: str,
    germline: str,
    gatk: str,
    chrom: str = None,
    ):
    if len(sample) != len(bam):
        raise Exception(f"The lenth of BAM not match with SAMPLE")
    config = BulidConfig(gatk = gatk,threads= 1)
    SNVCalling_process(sample,bam,fasta,outvcf,pon,germline,config,chrom)



@main.command(name="CallBack",no_args_is_help=True)
@click.option("-b","--bam", type=click.Path(exists=True, file_okay=True, dir_okay=False), required=True, help="Bam file with inde(ex: -b example.bam1 -b example.bam2)")
@click.option("-v","--vcf", type=click.Path(exists=True, file_okay=True, dir_okay=False), required = True,  help = "Output VCF")
@click.option("-s","--sample",required=True, help="Bam file with inde(ex: -b example.bam1 -b example.bam2)")
@click.option("-o","--outdir", type = str, required = True,  help = "Output dir")
@click.option('--tmpdir',type=str,help="Specify a temporary directory")
@click.option('--only_autosome', is_flag=True,default=False,help="Only analysis autosome")
@click.option('--stereo', is_flag=True,default=False, help="Please check your data is stereo(barcode is Cx and Cy) or not")
@click.option("-c","--barcode", type = str, default = 'CB', help = "Cell Barcode in Bam file (ex. CR for 10X Genomics)")
@click.option("-u","--umi", type = str, default ='UB', help = "UMI(Molecular Barcodes) in Bam file (ex. UR for 10X Genomics)")
@click.option("-@","--threads", type = int, default = 10, help = "Sets the number of threads")
@click.option("--qmap", type = int, default = 10, help = "Sets the number of qmap")
@click.option("-x",'--xsetoff',default = None,type = int, help="gem_x + x_offset = bam_x")
@click.option("-y",'--ysetoff',default = None,type = int, help="gem_y + y_offset = bam_y")
@click.option('--binsize',type = int,default = None, help="binsize")
@click.option('--removetmp', is_flag=True,default=False, help="Remove all tmp file")
def CallBack(
    bam: str,
    vcf: str,
    sample: str,
    outdir: str,
    stereo: bool = True,
    only_autosome: bool = True,
    barcode: str = 'CB',
    umi: str = 'UB',
    tmpdir: str = None,
    qmap: int = 10,
    threads: int = 10,
    xsetoff: int = None,
    ysetoff: int = None,
    binsize : int = None,
    removetmp: bool = False,
):
    CallBack_process(bam,only_autosome,sample,outdir,tmpdir,vcf,barcode,umi,qmap,stereo,threads = threads)
    BulidResult(sample = sample ,sample_tmp_path = tmpdir,sample_matrix_path = outdir,x_offset = xsetoff,y_offset = ysetoff,bin_size = binsize,removetmp = removetmp)