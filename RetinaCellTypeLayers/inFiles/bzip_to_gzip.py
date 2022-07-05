import glob
import os


DEL_FASTQ = True

bz2_files = sorted(glob.glob("*.bz2"))

for file in bz2_files:
    print("Decompressing file %s" % file)
    os.system("bzip2 -dk %s" % file)

    print("Compressing file...")
    os.system("gzip -c %s > %s" % (
        file.split(".bz2")[0],
        file.split(".bz2")[0] + ".gz"
        )
    )

    if DEL_FASTQ:
        os.system("rm %s" % file.split(".bz2")[0])
