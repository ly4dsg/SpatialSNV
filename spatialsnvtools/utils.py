import pysam
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from subprocess import run,PIPE
import subprocess


from .log import logger_config

def BulidConfig(samtools = None,gatk = None,picard = None,threads = 1):
    logger = logger_config("CheckSoftware")
    if samtools != None:
        samtools = check_software('samtools',samtools,logger = logger)
    if picard != None:
        picard = check_software('picard',picard,logger = logger)
    gatk = check_software('gatk',gatk,logger = logger)
    config = Config.set_path(picard = picard, gatk = gatk, samtools = samtools, threads = threads )
    return config

def check_software(software_name, software_path,logger = None):
    if not software_path:
        try:
            output = subprocess.check_output(["which", f"{software_name}"], stderr=subprocess.STDOUT)
            output = output.decode().rstrip('\n')
            if logger != None:
                logger.info(f"Check {software_name} : Find in {output}")
            return output
        except subprocess.CalledProcessError as e:
            if logger != None:
                logger.debug(f"{software_name} may not be installed, please use install {software_name}\n Or set --{software_name} PATH_TO_{software_name}")
    else:
        if os.path.exists(software_path):
            if logger != None:
                logger.info(f"Check {software_name} : Find in {software_path}")
            return software_path
        else:
            raise FileNotFoundError(f"{software_name} not exist")

@dataclass
class Config:
    picard: Path
    gatk: Path
    samtools: Path
    threads: int
    def set_path(picard,gatk,samtools,threads):
        return Config(picard = picard,
                     gatk = gatk,
                     samtools = samtools,
                     threads = threads)
    
    def picard_cmd(self, command: str, tmp_dir: Path, mem: str = "100g"):
        return [
            "java",
            f"-Djava.io.tmp_dir={tmp_dir}",
            f"-Xmx{mem}",
            "-XX:+UseParallelGC",
            "-XX:GCTimeLimit=20",
            "-XX:GCHeapFreeLimit=10",
            "-jar",
            self.picard,
            command,
            "--TMP_DIR",
            tmp_dir,
            "--VALIDATION_STRINGENCY",
            "SILENT",
            "--VERBOSITY",
            "WARNING",
            "--QUIET",
            "true",
        ]
    def gatk_cmd(self,command: str, tmp_dir: str = None,mem: str = "100g"):
        if tmp_dir != None:
            return [
                self.gatk,
                "--java-options",
                f"-Xmx{mem}",
                command,
                "--tmp-dir",
                tmp_dir,
            ]
        else:
            return [
                self.gatk,
                "--java-options",
                f"-Xmx{mem}",
                command,
            ]
    def samtools_cmd(self, command: str):
        return [
            self.samtools,
            command,
            "-@",
            self.threads,
        ]

def RunCMD(cmd: list[Any],name: str,id: str,logger):
    cmd = [str(arg) for arg in cmd]
    proc = run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise Exception(f"Error running {name}:\n\t{proc.stderr}")
    else:
        logger.info(f"{id}-{name}: completed")

def check_bam_index(bam):
    bai_file = f"{bam}.bai"
    logger = logger_config("BAMIndexCheck")
    if not os.path.exists(bai_file):
        raise Exception(f"{bam} without index")
    logger.info(f"Find BAM Index: {bai_file}")
    

def LoadChromFromBAM(bam, only_autosome = True):
    logger = logger_config("LoadChromFromBAM")
    check_bam_index(bam)
    inputfile = pysam.AlignmentFile(bam,"rb")
    chrom_tuple = inputfile.references
    inputfile.close()
    if only_autosome:
        chrom_tuple_tmp = [chrom for chrom in chrom_tuple if chrom.startswith('chr') and 'Un' not in chrom ]
        if len(chrom_tuple_tmp) == 0 :
            chrom_tuple_tmp = [chrom for chrom in chrom_tuple if "K" not in chrom and "G" not in chrom and "L" not in chrom and "J" not in chrom and "V" not in chrom and "Q" not in chrom and "N" not in chrom and "h" not in chrom and '.' not in chrom and 'Z' not in chrom]
        chrom_str = "[" + ','.join(chrom_tuple_tmp) + "]"
        logger.info(f"Find CHROM: {chrom_str}")
        return chrom_tuple_tmp
    else:
        chrom_str = "[" + ','.join(chrom_tuple) + "]"
        logger.info(f"Find CHROM: {chrom_str}")
        return chrom_tuple

def DirCheck(path, create=False):
    if os.path.isdir(path):
        return True
    else:
        if create:
            try:
                os.makedirs(path)
                return True
            except OSError:
                raise OSError(f"Failed to create {path}")
        return False