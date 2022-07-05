import os
import glob

COND_MAP = {'TAAGGCGA': 'LS8A', #'N701',
            'CGTACTAG': 'LS8B', #'N702',
            'AGGCAGAA': 'LS8C', #'N703',
            'TCCTGAGC': 'LS8D', #'N704',
            'GGACTCCT': 'LS8E', #'N705',
            'TAGGCATG': 'LS8F', #'N706',
            'CTCTCTAC': 'LS8H'} #'N707'}
            
def parse_index():

    for key in COND_MAP:
        pattern = "^@.*%s" % key
        file = "Undetermined_S0_L001_R1.fastq.gz" # Lane 1
        outfile = "%s_iParsed_L001_R1.fastq.gz" % COND_MAP.get(key)

        cmd = "zgrep -A3 \"%s\" %s | grep -v \"^--\" | gzip -c > %s" % (
            pattern, file, outfile
            )

        print("Running this cmd: %s" % cmd)
        os.system(cmd)

    for key in COND_MAP:
        pattern = "^@.*%s" % key
        file = "Undetermined_S0_L002_R1.fastq.gz" # Lane 2
        outfile = "%s_iParsed_L002_R1.fastq.gz" % COND_MAP.get(key)

        cmd = "zgrep -A3 \"%s\" %s | grep -v \"^--\" | gzip -c > %s" % (
            pattern, file, outfile
            )

        print("Running this cmd: %s" % cmd)
        os.system(cmd)

if __name__ == '__main__':
    parse_index()