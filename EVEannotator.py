# EVEannotator is designed for the identification of endogenous viral elements in animal genomes.

import os
import sys
import getopt

help_info = ''

def main(argv):
    try:
        opts,args = getopt.getopt(sys.argv[1:], '-h-v-i:-o:-d:-t:-s:-e:', ['help', 'version','genome =','outdir =','database =','threads =','sensitivity =','taxon-exclude'])
    except:
        print("Error")
    
    for opt, arg in opts:
        if opt in ['-h', '--help']:
            print(help_info)
        elif opt in ['-v', '--version']:
            print("EVEfinder: Version 1.0")
        elif opt in ['-i', '--genome']:
            genome = arg
        elif opt in ['-o', '--outdir']:
            result_path = arg
        elif opt in ['-d', '--database']:
            nr_db = arg
        elif opt in ['-t', '--threads']:
            threads = arg
        elif opt in ['-s', '--sensitivity']:
            sensitivity = arg
        elif opt in ['-e', '--taxon-exclude']:
            taxon_exclude = arg

    
    # call other function
    blastx_result = blastx(genome, nr_db, threads, sensitivity, result_path)
    blastx_seq = cluster_seq(blastx_result)
    blastp_tbl = blastp(blastx_seq, nr_db, threads, taxon_exclude)
    blastp_result_virus = getEVE(blastp_tbl)
    taxonkit(blastp_result_virus)
        
def blastx(genome_file, nr_db, threads, sensitivity, result_path):

    blastx = 'diamond blastx --db {nr} --query {genome} --out {blastx_result} --taxonlist 10239 --{sensitivity} --threads {thread} --evalue 0.1 -k0 --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue qseq_translated staxids sscinames sskingdoms skingdoms sphylums stitle'
    genome_file_name = genome_file.split('/')[-1]
    blastx_result = result_path + genome_file_name + '.blastx.tbl'
    blast_command = blastx.format(genome = genome_file, nr = nr_db, blastx_result = blastx_result, thread = threads, sensitivity = sensitivity)

    os.system(blast_command)
    return blastx_result

        
def cluster_seq(blastx_result):
    ori_fa = blastx_result + '.fa'
    ori_bed = blastx_result + '.bed'

    with open(blastx_result, 'r') as f1:
        lines = f1.readlines()
    
    f2 = open(ori_fa, 'w')
    f3 = open(ori_bed, 'w')

    for line in lines:
        if int(line.split()[4]) < int(line.split()[5]):
            seq_name = line.split()[0] + '__plus__' + line.split()[4] + '__' + line.split()[5] 
            bed_line = line.split()[0] + '\t' + line.split()[4] + '\t' + line.split()[5] + '\t' + seq_name + '\t' + '0' + '\t' + '+' + '\n'
        else:
            seq_name = line.split()[0] + '__minus__' + line.split()[5] + '__' + line.split()[4]
            bed_line = line.split()[0] + '\t' + line.split()[5] + '\t' + line.split()[4] + '\t' + seq_name + '\t' + '0' + '\t' + '-' + '\n'

        seq = '>' + seq_name +'\n' + line.split()[9] + '\n'
        f2.write(seq)
        f3.write(bed_line)


    # bedtools cluster
    new_bed = ori_bed[:-4] + '.grouped.bed'
    bedtools_command = 'bedtools cluster -s -i ' + ori_bed + ' > ' + new_bed
    os.system(bedtools_command)

    # extract the longest sequence in classified groups
    with open(new_bed) as f4:
        bed_lines = f4.readlines()
    
    bed_list = []
    n = 0

    for bed_line in bed_lines:
        if n == 0:
            bed_list.append(bed_line)
        else:
            new_groupID = bed_line.split('\t')[6]
            last_groupID = bed_list[-1].split('\t')[6]

            if new_groupID != last_groupID:
                bed_list.append(bed_line)
                
            else:
                last_length = int(bed_list[-1].split('\t')[2]) - int(bed_list[-1].split('\t')[1])
                new_length = int(bed_line.split('\t')[2]) - int(bed_line.split('\t')[1])

                if last_length <= new_length:
                    bed_list[-1] = bed_line

        n += 1
    
    # extract specified sequences from '.fa' according to names
    clustered_bed = ori_bed[:-4] + '.clustered.bed'

    with open('temp.bed', 'w') as f5:
        with open(clustered_bed, 'w') as f4:
            for i in bed_list:
                f4.write(i)

                f5_line = i.split('\t')[3] + '\n'
                f5.write(f5_line)

    clustered_fa = blastx_result + '.clustered.fa'
    seqkit_grep_rmdup_command = 'seqkit grep -n -f temp.bed ' + ori_fa + ' | seqkit rmdup >' + clustered_fa
    rm = 'rm temp.bed'
    os.system(seqkit_grep_rmdup_command)
    os.system(rm)

    return clustered_fa

def blastp(blastx_seq, nr_db, threads, taxon_exclude):
    blastp_result = blastx_seq[:-3] + '.blastp.tbl'
    query = blastx_seq
    nr = nr_db
    blastp = "diamond blastp --db {nr} --query {query} --out {blastp_result} --threads {threads} --taxon-exclude {exclude_taxon} --evalue 0.00001 -k 1 --outfmt 6 qseqid sseqid pident length evalue staxids sscinames sskingdoms stitle qseq"
    os.system(blastp.format(query = query, blastp_result = blastp_result, nr = nr, threads = threads, exclude_taxon = taxon_exclude))

    return blastp_result

def getEVE(blastp_tbl):
    blastp_result_virus = blastp_tbl + '.virus.tbl'

    grep = 'grep -i "Viruses" {blastp_result} > {blastp_result_virus}'
    os.system(grep.format(blastp_result = blastp_tbl, blastp_result_virus = blastp_result_virus))

    return(blastp_result_virus)

def taxonkit(blastp_result_virus):
    new_blastp_result_virus = blastp_result_virus[:-4] + '.classified.tbl'
    taxonkit_cmd = "awk '{print $6}' " + blastp_result_virus + " | awk -F ';' '{print $1}' | taxonkit lineage | taxonkit reformat -f '{f}' >temp.txt"
    merge = "csvtk join -t -T -f '6;1' --keep-unmatched " + blastp_result_virus + " temp.txt | uniq >" + new_blastp_result_virus

    os.system(taxonkit_cmd)
    os.system(merge)
    os.system('rm temp.txt')

if __name__ == "__main__":
    main(sys.argv[1:])

