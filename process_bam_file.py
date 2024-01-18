import pandas as pd
import time
import subprocess
import pysam

def main():
    #genome_ids = range(2, 19)
    genome_ids = [8]
    start_time = time.time()
    query_start_range = 1000
    query_end_range = 2000

    genome_id_to_ref_name = {}

    for genome_id in genome_ids:
        samfilename = f'camisim_sample/bam/Genome{genome_id}.0.sam'
        bamfilename = f'camisim_sample/bam/Genome{genome_id}.0.bam'
        subprocess.call(['samtools', 'view', bamfilename, '-o', samfilename])

        # Read the SAM file as a pandas DataFrame
        df = pd.read_csv(samfilename, sep='\t', header=None)
        df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
        
        mask1 = (df['POS'] >= query_start_range) & (df['POS'] <= query_end_range)
        mask2 = (df['POS']+df['TLEN'] >= query_start_range) & (df['POS']+df['TLEN'] <= query_end_range)
        df = df[mask1 | mask2]
        num_rows = df.shape[0]
        print(f"Number of rows: {num_rows}")
        print(df)

        genome_id_to_ref_name[genome_id] = list(df['RNAME'].unique())

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time using approach 1: {elapsed_time} seconds")

    
    start_time = time.time()
    for genome_id in genome_ids:
        bamfilename = f'camisim_sample/bam/Genome{genome_id}.0.bam'
        bamfile = pysam.AlignmentFile(bamfilename, "rb")

        ref_names = genome_id_to_ref_name[genome_id]
        num_matches = 0
        for ref_name in ref_names[:1]:
            list_of_matches = list(bamfile.fetch(ref_name, query_start_range, query_end_range))
            num_matches += len(list_of_matches)
            print(f'Reference name: {ref_name}')
            for match in list_of_matches:
                print(str(match))

            coverage_list = [match.get_tag('RC') for match in list_of_matches]
            print(f'Coverage list: {coverage_list}')

            coverage_list_using_pileupcolumn = [pileupcolumn.n for pileupcolumn in bamfile.pileup(ref_name, query_start_range, query_end_range)]
            print(coverage_list_using_pileupcolumn)

        print(f"Number of matches: {num_matches}")
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time using approach 2: {elapsed_time} seconds")


if __name__ == "__main__":
    main()
