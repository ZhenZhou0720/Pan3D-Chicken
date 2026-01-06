#!/usr/bin/env python3

import argparse
import pandas as pd
import os

 
def get_args():
    parser = argparse.ArgumentParser(description='Parse and classify structural variants')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file path')
    parser.add_argument('-te', '--te', required=True, help='TE bed file path')
    parser.add_argument('-trf', '--trf', required=True, help='TRF gff file path')
    parser.add_argument('-indel', '--indel', required=True, help='InDel fragment file path')
    parser.add_argument('-delins', '--delins', required=True, help='SV.INS_DEL bed file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    return parser.parse_args()


def run_blast(chr_id, pos1, pos2):
    reffa = 'GGswu.v23.fa'
    samtools = '/public/software/samtools/samtools-1.12/bin/samtools'
    makeblastdb = '/public/software/blast/ncbi-blast-2.11.0+/bin/makeblastdb'
    blastn = '/public/software/blast/ncbi-blast-2.11.0+/bin/blastn'
    
    left_pos = max(0, int(pos1) - 600)
    right_pos = int(pos2) + 600
    
    os.system(f'{samtools} faidx {reffa} {chr_id}:{left_pos}-{pos1} >left.fa 2>/dev/null')
    os.system(f'{samtools} faidx {reffa} {chr_id}:{pos2}-{right_pos} >right.fa 2>/dev/null')
    
    os.system(f'{makeblastdb} -in left.fa -dbtype nucl 1>/dev/null 2>&1')
    os.system(f'{blastn} -db left.fa -query right.fa -out bltoutfile -evalue 1 -outfmt 6 1>/dev/null 2>&1')
    
    homo = 0
    if os.path.getsize('bltoutfile') > 0:
        with open('bltoutfile', 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                homo = int(parts[3]) - int(parts[4]) - int(parts[5])
    
    os.system('cat bltoutfile >>all.bltoutfile 2>/dev/null')
    os.system('rm -f left.fa* right.fa* bltoutfile')
    return homo


def check_overlap(snp_frag, frag_df):
    frag_df = frag_df.rename(columns={0:'chr', 1:'start', 2:'end'})
    frag_df[['start', 'end']] = frag_df[['start', 'end']].astype(int)
    
    chr_filtered = frag_df[frag_df['chr'] == snp_frag[0]]
    overlap = chr_filtered[(chr_filtered['start'] <= snp_frag[1]) & (chr_filtered['end'] >= snp_frag[1])]
    
    ins_len = overlap['end'].iloc[0] - overlap['start'].iloc[0] if not overlap.empty else 0
    return not overlap.empty, ins_len


def main():
    args = get_args()
    
    vcf = pd.read_csv(args.input, comment='#', sep='\t', header=None).iloc[:, :9]
    vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    
    te_df = pd.read_csv(args.te, sep='\t', header=None)
    trf_df = pd.read_csv(args.trf, comment='#', sep='\t', header=None).iloc[:, [0, 3, 4]]
    indel_df = pd.read_csv(args.indel, sep='\t', header=None)
    delins_df = pd.read_csv(args.delins, sep='\t', header=None)
    
    with open(args.output, 'w') as f:
        f.write("CHROM\tStart\tEnd\tID\tREF\tALT\tType\n")
    
    for _, row in vcf.iterrows():
        chr_name = row['CHROM']
        start_pos = int(row['POS'])
        
        if 'END=' not in row['INFO'] or 'SVTYPE=' not in row['INFO']:
            continue
        end_pos = int(row['INFO'].split('END=')[1].split(';')[0])
        sv_type = row['INFO'].split('SVTYPE=')[1].split(';')[0]
        
        var_frag = [chr_name, start_pos, end_pos]
        var_type = 'N/A'
        
        if sv_type == 'DEL':
            if check_overlap(var_frag, te_df)[0]:
                var_type = 'DEL_TEI'
            elif check_overlap(var_frag, trf_df)[0]:
                var_type = 'DEL_VNTR'
            elif check_overlap(var_frag, indel_df)[0]:
                ins_len = check_overlap(var_frag, indel_df)[1]
                var_type = 'DEL_FoSteS/MMBIR' if ins_len > 10 else 'DEL_NHEJ'
            else:
                homo = run_blast(chr_name, start_pos, end_pos)
                if homo <= 2:
                    var_type = 'DEL_NHEJ'
                elif 2 < homo <= 50:
                    var_type = 'DEL_alt-EJ'
                elif 50 < homo <= 200:
                    var_type = 'DEL_SSA'
                elif homo > 200:
                    var_type = 'DEL_NAHR'
        
        elif sv_type == 'INS':
            if check_overlap(var_frag, te_df)[0]:
                var_type = 'INS_TEI'
            elif check_overlap(var_frag, trf_df)[0]:
                var_type = 'INS_VNTR'
        
        elif sv_type == 'INV':
            indel_overlap, ins_len = check_overlap(var_frag, indel_df)
            if indel_overlap:
                if check_overlap(var_frag, delins_df)[0]:
                    var_type = 'INV_FoSteS/MMBIR' if ins_len > 10 else 'INV_NHEJ'
            else:
                homo = run_blast(chr_name, start_pos, end_pos)
                if homo <= 2:
                    var_type = 'INV_NHEJ'
                elif 2 < homo <= 200:
                    var_type = 'INV_alt-EJ'
                elif homo > 200:
                    var_type = 'INV_NAHR'
        
        elif sv_type == 'TRA':
            if check_overlap(var_frag, te_df)[0]:
                var_type = 'TRA_TEI'
            elif check_overlap(var_frag, indel_df)[0]:
                ins_len = check_overlap(var_frag, indel_df)[1]
                var_type = 'TRA_FoSteS/MMBIR' if ins_len > 10 else 'TRA_NHEJ'
            else:
                homo = run_blast(chr_name, start_pos, end_pos)
                if homo <= 2:
                    var_type = 'TRA_NHEJ'
                elif 2 < homo <= 200:
                    var_type = 'TRA_alt-EJ'
                elif homo > 200:
                    var_type = 'TRA_NAHR'
        
        with open(args.output, 'a') as f:
            f.write(f"{chr_name}\t{start_pos}\t{end_pos}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{var_type}\n")


if __name__ == '__main__':
    main()