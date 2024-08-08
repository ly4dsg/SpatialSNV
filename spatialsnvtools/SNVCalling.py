

from .utils import logger_config,DirCheck,RunCMD


def SNVCalling_process(sample,bam,fasta,outvcf,pon,germline,config,chrom = None):
    logger = logger_config("SNVCalling")
    mutect2_cmd = config.gatk_cmd("Mutect2")

    for single_sample,single_bam in zip(sample,bam):
        logger.debug(f"Match BAM Path: [{single_sample}]-{single_bam}")
        mutect2_cmd.extend(
            [
                "-I",
                single_bam
            ]
        )
    raw_vcf = outvcf[:-7] + '.raw.vcf.gz'
    mutect2_cmd.extend(
        [
            "-R",
            fasta,
            "--minimum-mapping-quality",
            10,
            "-tumor",
            "tumor",
            "--native-pair-hmm-threads",
            "10",
            "--output",
            raw_vcf,
            "--germline-resource",
            f"{germline}",
            "--panel-of-normals",
            f"{pon}",
        ]
    )
    if chrom != None:
        mutect2_cmd.extend(
            [
                "-L",
                f"{chrom}",]
        )
    sample_str = ','.join(sample)
    RunCMD(mutect2_cmd,"Mutect2","",logger)

    filter_cmd = config.gatk_cmd("FilterMutectCalls")
    filter_cmd.extend(
        [
            "-V",
            raw_vcf,
            "-R",
            fasta,
            "-O",
            outvcf,
            "--min-median-base-quality",
            10
        ]
    )
    RunCMD(filter_cmd,"FilterMutectCalls","",logger)