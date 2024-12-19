import GEOparse

# wang 2022 (Diaz group)
outdir = "/home/jiehoonk/project/crisprscreen/data/wang2022_NC"
gse = GEOparse.get_GEO(geo="GSE174554", destdir=outdir, how = 'full')
gse.download_supplementary_files(directory = outdir, download_sra = False, nproc = 64)

## check metadata

