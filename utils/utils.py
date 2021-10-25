import os
import glob

def globBAMfiles(**kwargs):
    if kwargs.get("type"):
        filelist = glob.glob("*_" + kwargs.get("type") + ".bam")
    else:
        filelist = glob.glob("*[!dedup][!sorted][!featurecounts].bam")
    return filelist

def getBAMname(file, **kwargs):
    if kwargs.get("type"):
        name = file.split("_" + kwargs.get("type") + ".bam")[0]
    else:
         name = file.split(".bam")[0]
    return name

def sortBAM(file):
    fname = getBAMname(file)
    print("Sorting File: %s" % (file))
    os.system('samtools sort %s -o %s' % \
             (file , fname + '_sorted.bam')
             )

    sortedfname = fname + '_sorted.bam'
    return sortedfname

def indexBAM(file, **kwargs):
    fname = getBAMname(file, type=kwargs.get("type"))
    print("Indexing File: %s" %file)
    os.system('samtools index %s' %file)

def umiDedup(file):
    fname = getBAMname(file, type="sorted")
    umistr = 'umi_tools dedup --per-gene --gene-tag=XT --assigned-status-tag=XS -I %s -S %s --log %s' % \
             (file, fname + "_dedup.bam", fname + '_deduplog.log')
    print('\n' + 'Running: ' + umistr + '\n')
    os.system(umistr)
