# SpatialSNV

A novel method for calling and analyzing SNVs from spatial transcriptomics data.

---

We divided the process of calling mutations from spatial transcriptomics data into two parts: **Data Preprocessing** and **Data Analysis**.  

All analyzed jupyter notebooks are saved in the `article` folder

## Install

To install `spatialsnv`, use pip:

```bash
pip install spatialsnv
```
We recommend using Python version 3.10.14. You also need to install the following tools: **samtools**,**gatk**,**picard**

## Data Preprocessing

### Splitting BAM File by Chromosome for Speed Up (Optional)

```bash
spatialsnvtools SplitChromBAM -b demo.bam –s demo –o demo_split -@ 10 –only_autosome

```

**Options:**
- `-b, --bam FILE`  
  BAM file that needs to be split by chromosome **[required]**
- `-s, --sample TEXT`  
  Sample ID **[required]**
- `-o, --outdir TEXT`  
  Output directory for the split BAM files **[required]**
- `-@, --threads INTEGER`  
  Sets the number of threads
- `--only_autosome`  
  Only analyze autosomes
- `--help`  
  Show this message and exit.


### Mutation Calling Data Preprocessing

> 10x or drop-seq
```bash
spatialsnvtools PerpareBAMforCalling -b demo.bam -o process_out -s demo -c 'CR' -u 'UR' -@ 10 --fasta GRCh38.p12.genome.fa --dbsnp dbsnp.chr9.hg38.vcf.gz --removetmp --picard $pathtopicard/picard.jar --gatk $pathtogatk/gatk --samtools $pathtosamtools/samtools
```
> stereo-seq
```bash
spatialsnvtools PerpareBAMforCalling -b demo.stereo.bam -o stereo_process_out -s stereo_demo --stereo --gem demo.gem.gz -x 0 -y 0 -@ 10 --fasta GRCh38.p12.genome.fa --dbsnp dbsnp.chr9.hg38.vcf.gz --removetmp --picard $pathtopicard/picard.jar --gatk $pathtogatk/gatk --samtools $pathtosamtools/samtools
```

**Options:**
- `-b, --bam FILE`  
  BAM file that needs preprocessing **[required]**
- `-o, --outdir TEXT`  
  Output directory for the preprocessing results **[required]**
- `-s, --proxy TEXT`  
  Sample ID **[required]**
- `-f, --fasta FILE`  
  Reference FASTA used for mutation calling **[required]**
- `-d, --dbsnp FILE`  
  dbSNP file used for BQSR **[required]**
- `-c, --barcode TEXT`  
  Cell Barcode in BAM file (e.g., `CR` for 10X Genomics)
- `-u, --umi TEXT`  
  UMI (Molecular Barcodes) in BAM file (e.g., `UR` for 10X Genomics)
- `--stereo`  
  Ensure that your data is stereo (barcode is `Cx` and `Cy`)
- `--gem TEXT`  
  GEM file matching your raw stereo BAM
- `-x, --xsetoff INTEGER`  
  `gem_x + x_offset = bam_x`
- `-y, --ysetoff INTEGER`  
  `gem_y + y_offset = bam_y`
- `--tmpdir TEXT`  
  Specify a temporary directory
- `-@, --threads INTEGER`  
  Sets the number of threads
- `--samtools TEXT`  
  Specify the path to `samtools`, if not specified, automatically detected
- `--picard TEXT`  
  Specify the path to `picard.jar`, if not specified, automatically detected
- `--gatk TEXT`  
  Specify the path to `GATK`, if not specified, automatically detected
- `--removetmp`  
  Remove all temporary files
- `--help`  
  Show this message and exit.
  
### SNV Calling on Preprocessed BAM Files
```bash
spatialsnvtools SNVCalling -b demo.processed.bam -s demo -o demo.vcf.gz -f GRCh38.p12.genome.fa --pon 1000g_pon.hg38.vcf.gz --germline af-only-gnomad.hg38.vcf.gz
```


**Options:**
- `-s, --sample TEXT`  
  Sample ID (e.g., `-b example.bam1 -b example.bam2`) **[required]**
- `-b, --bam FILE`  
  BAM file(s) with index (e.g., `-b example.bam1 -b example.bam2`) **[required]**
- `-o, --outvcf TEXT`  
  Output VCF file **[required]**
- `-f, --fasta FILE`  
  Reference FASTA file **[required]**
- `--pon FILE`  
  Panel of Normals (PON) file
- `--germline FILE`  
  Germline source file 
- `--gatk TEXT`  
  Specify the path to `GATK`
- `-L, --chrom TEXT`  
  Specify chromosome(s) to analyze
- `--help`  
  Show this message and exit.

### Traceback SNVs to Spatial Transcriptomics

> 10x or drop-seq
```bash
spatialsnvtools CallBack --bam demo.processed.bam --vcf demo.vcf.gz -o demo_matrix -s demo --tmpdir demo_tmp --only_autosome -c "CB"  -u "UB" -@ 1
```
> stereo-seq
```bash
spatialsnvtools CallBack --bam demo.stereo.bam --vcf demo.stereo.vcf.gz -o demo_stereo_matrix -s demo_stereo --tmpdir demo_stereo_tmp --stereo -x 0 -y 0 --binsize 100 --only_autosome -@ 1 --umi UM --removetmp
```

**Options:**
- `-b, --bam FILE`  
  BAM file(s) with index (e.g., `-b example.bam1 -b example.bam2`) **[required]**
- `-v, --vcf FILE`  
  VCF file for SNV data **[required]**
- `-s, --sample TEXT`  
  Sample ID **[required]**
- `-o, --outdir TEXT`  
  Output directory **[required]**
- `--tmpdir TEXT`  
  Specify a temporary directory
- `--only_autosome`  
  Only analyze autosomes
- `--stereo`  
  Ensure that your data is stereo (barcode is `Cx` and `Cy`)
- `-c, --barcode TEXT`  
  Cell Barcode in BAM file (e.g., `CR` for 10X Genomics)
- `-u, --umi TEXT`  
  UMI (Molecular Barcodes) in BAM file (e.g., `UR` for 10X Genomics)
- `-@, --threads INTEGER`  
  Sets the number of threads
- `--qmap INTEGER`  
  Sets the number of qmap
- `-x, --xsetoff INTEGER`  
  `gem_x + x_offset = bam_x`
- `-y, --ysetoff INTEGER`  
  `gem_y + y_offset = bam_y`
- `--binsize INTEGER`  
  Set the bin size
- `--removetmp`  
  Remove all temporary files
- `--help`  
  Show this message and exit.


## Data Analysis
To be determined
