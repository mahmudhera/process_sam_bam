import pandas as pd
import time
import subprocess
import pysam

def main():
    genome_ids = range(2, 19)
    start_time = time.time()
    query_start_range = 1000
    query_end_range = 2000

    for genome_id in genome_ids:
        samfilename = f'camisim_sample/bam/Genome{genome_id}.0.sam'
        bamfilename = f'camisim_sample/bam/Genome{genome_id}.0.bam'
        subprocess.call(['samtools', 'view', bamfilename, '-o', samfilename])

        # Read the SAM file as a pandas DataFrame
        df = pd.read_csv(samfilename, sep='\t', header=None)
        df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
        df = df[df['TLEN'] > 0]
        print(df)

        df = df[df['POS'] >= query_start_range]
        df = df[df['PNEXT'] <= query_end_range]
        num_rows = df.shape[0]
        print(f"Number of rows: {num_rows}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time using approach 1: {elapsed_time} seconds")

    """
    start_time = time.time()
    for genome_id in genome_ids:
        bamfilename = f'camisim_sample/bam/Genome{genome_id}.0.bam'
        bamfile = pysam.AlignmentFile(bamfilename, "rb")
    """
        

if __name__ == "__main__":
    main()
