import GEOparse
import os
njobs = 12

# wang 2022 (Diaz group)
outdir = "../data/wang2022_NC"
gse = GEOparse.get_GEO(geo="GSE174554", destdir=outdir, how = 'full')
gse.download_supplementary_files(directory = outdir, download_sra = False, nproc = njobs)

# park 2024
## 1. snRNA-seq
geo="GSE189650"
outdir = "../data/park2024_PO/RNA"
os.makedirs(outdir, exist_ok = True)
gse = GEOparse.get_GEO(geo=geo, destdir=outdir, how = 'full')
gse.download_supplementary_files(directory = outdir, download_sra = False, nproc = njobs)

## 2. scATAC-seq
geo="GSE157910"
outdir = "../data/park2024_PO/ATAC" # sn vs sc?
os.makedirs(outdir, exist_ok = True)
gse = GEOparse.get_GEO(geo=geo, destdir=outdir, how = 'full')
gse.download_supplementary_files(directory = outdir, download_sra = False, nproc = njobs)
