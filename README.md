# SpatialSNV

A novel method for calling and analyzing SNVs from spatial transcriptomics data.

---

We divided the process of calling mutations from spatial transcriptomics data into two parts: **Data Preprocessing** and **Data Analysis**.

## Install
To be determined

## Data Preprocessing
To be determined

### Splitting BAM File by Chromosome for Speed Up (Optional)

```bash
spatialsnvtools SplitChromBAM -b demo.bam –s demo –o demo_split -@ 10 –only_autosome

```

**Usage:** `spatialsnvtools SplitChromBAM [OPTIONS]`

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
  

## Data Analysis
To be determined
