import pysam
import os 
import shutil
import pandas as pd
import numpy as np

from .utils import logger_config,DirCheck,RunCMD
##########################################
####  BarcodeUMIBinding Part 
##########################################
def Code_barcode(barcode):
    code = {1:'AA',2:'AT',3:'AC',4:'AG',5:'TA',6:'TT',7:'TC',8:'TG',9:'CA',0:'CC'}
    barcode_list = list()
    for m in barcode:
        barcode_list.append(code[int(m)])
    return ''.join(barcode_list)

def Decoding_umi(umi):
    umi_dict = {0:'A',1:'C',2:'G',3:'T'}
    umi_num = np.base_repr(int(umi,16),base=4)
    umi_list = list()
    for i in umi_num:
        umi_list.append(umi_dict[int(i)])
    return ''.join(umi_list)

def LoadLocation(gem,xoff,yoff):
    is_gzip = gem.endswith('.gz')
    if is_gzip:
        df = pd.read_csv(gem,sep='\t',compression='gzip', comment='#')
    else:
        df = pd.read_csv(gem,sep='\t', comment='#')
    df = df.drop_duplicates(['x','y'])
    df['bam_x'] = df['x'] + xoff
    df['bam_y'] = df['y'] + yoff
    df['location'] = df['bam_x'].map(str)+'-'+df['bam_y'].map(str)
    uniq_coo = pd.DataFrame(df['location'])
    uniq_coo['num'] = 1
    location_dict = dict(zip(uniq_coo['location'],uniq_coo['num']))
    return location_dict

def BarcodeUMIBinding(bam,outbam,stereo = False,barcode_tag = "CR",umi_tag = 'UR',location_dict = None):

    inputfile = pysam.AlignmentFile(bam, "rb")
    output = pysam.AlignmentFile(outbam, 'wb', template = inputfile)

    if not stereo:
        for read in inputfile:
            _barcode = read.has_tag(barcode_tag)
            _umi = read.has_tag(umi_tag)
            if _barcode and _umi:
                barcode = read.get_tag(barcode_tag)
                umi = read.get_tag(umi_tag)
                BI = barcode + "-" + umi
                read.set_tag('LY', BI)
                output.write(read)
    else:
        for read in inputfile:
            x = str(read.get_tag('Cx'))
            y = str(read.get_tag('Cy'))
            umi = str(read.get_tag('UR'))

            if location_dict != None:
                location = x + '-' + y
                intissue = location_dict.get(location,False)
                if intissue:
                    new_barcode = Code_barcode(x) + '-' + Code_barcode(y)
                    new_umi = Decoding_umi(umi)
                    encode_barcode = new_barcode + '-' + new_umi
                    read.set_tag('LY', encode_barcode)
                    read.set_tag('CB', new_barcode)
                    read.set_tag('UM', new_umi)
                    output.write(read)
    output.close()
    inputfile.close()

##########################################
####  GATK process Part
##########################################
def GATKprocess(config,outdir,tmpdir,proxy,fasta,dbsnp):
    logger = logger_config("GATKprocess")
    Markduplicate_cmd = config.picard_cmd("MarkDuplicates", tmpdir)
    Markduplicate_cmd.extend(
            [   
                "-I",
                f"{tmpdir}/{proxy}.tag.bam",
                "-O",
                f"{tmpdir}/{proxy}.dup.bam",
                "-M",
                f"{tmpdir}/{proxy}.txt",
                "--TMP_DIR",
                f"{tmpdir}",
                "--REMOVE_SEQUENCING_DUPLICATES",
                "true",
                "--BARCODE_TAG",
                "LY",
            ]
        )
    RunCMD(Markduplicate_cmd, "MarkDuplicates",f"{proxy}",logger)

    Addsm_cmd = config.picard_cmd("AddOrReplaceReadGroups", tmpdir)
    Addsm_cmd.extend(
        [
            "-I",
            f"{tmpdir}/{proxy}.dup.bam",
            "-SM",
            "tumor",
            "-PL",
            "iDrop",
            "-LB",
            "lib1",
            "-PU",
            "unit1",
            "-ID",
            "3",
            "-O",
            f"{tmpdir}/{proxy}.sm.bam",
        ]
    )
    RunCMD(Addsm_cmd,"AddOrReplaceReadGroups",f"{proxy}",logger)

    SplitNCigarReads_cmd = config.gatk_cmd("SplitNCigarReads",tmpdir)
    SplitNCigarReads_cmd.extend(
        [
            "-R",
            fasta,
            "-I",
            f"{tmpdir}/{proxy}.sm.bam",
            "-O",
            f"{tmpdir}/{proxy}.nc.bam",
            "--max-reads-in-memory",
            300000,
        ]
    )
    RunCMD(SplitNCigarReads_cmd,"SplitNCigarReads",f"{proxy}",logger)

    BaseRecalibrator_cmd = config.gatk_cmd("BaseRecalibrator")
    BaseRecalibrator_cmd.extend(
        [
            "-R",
            fasta,
            "-I",
            f"{tmpdir}/{proxy}.nc.bam",
            "-O",
            f"{tmpdir}/{proxy}.base",
            "--use-original-qualities",
            "--known-sites",
            dbsnp,
        ]
    )
    RunCMD(BaseRecalibrator_cmd,"BaseRecalibrator",f"{proxy}",logger)

    ApplyBQSR_cmd = config.gatk_cmd("ApplyBQSR")
    ApplyBQSR_cmd.extend(
        [
            "-R",
            fasta,
            "-I",
            f"{tmpdir}/{proxy}.nc.bam",
            "-O",
            f"{outdir}/{proxy}.rdfcall.bam",
            "-bqsr",
            f"{tmpdir}/{proxy}.base",
            "--add-output-sam-program-record",
            "--use-original-qualities",
        ]
    )
    RunCMD(ApplyBQSR_cmd,"ApplyBQSR",f"{proxy}",logger)


def PerpareBAMforCalling_process(bam,outdir,proxy,barcode_tag,umi_tag,stereo,gem,xsetoff,ysetoff,config,fasta,dbsnp,tmpdir,removetmp):


    logger = logger_config("PerpareBAMforCalling")
    if tmpdir == None:
        tmpdir = os.path.join(outdir,"tmp")
    
    logger.info(f"Tmp dir: {tmpdir}")
    DirCheck(tmpdir,create = True)
    
    if stereo:
        LocationDict = LoadLocation(gem,xsetoff,ysetoff)
    else:
        LocationDict = None

    TagBAM = f"{tmpdir}/{proxy}.tag.bam"
    BarcodeUMIBinding(
        bam = bam,
        outbam = TagBAM,
        stereo = stereo,
        barcode_tag = barcode_tag,
        umi_tag = umi_tag,
        location_dict = LocationDict
        )
    logger.info(f"BarcodeUMIBinding: Finish")

    GATKprocess(
        config = config,
        outdir = outdir,
        tmpdir = tmpdir,
        proxy = proxy,
        fasta = fasta,
        dbsnp = dbsnp
        )

    if removetmp:
        logger.info(f"Remove Tmp Dir: {tmpdir}")
        shutil.rmtree(tmpdir)
