import numpy as np
import os
import numpy as np
import random
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import lzma

random.seed(10)

class Simulation():


    def __init__(self):

        self.depth = 50
        self.reads_len = 150
        self.chosen_genome_num = 10

        # self.original_genome_list = "/mnt/d/breakpoints/assembly/simulation/ecoli/fna.list"
        # self.reference = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/ecoli_ref.fna"
        # self.sim_dir = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/"
        # self.res_dir = "/home/wangshuai/assembly_result/ecoli/"

        self.original_genome_list = "/mnt/d/breakpoints/assembly/simulation/assembly_test/Staphylococcus_lugdunensis/Staphylococcus_lugdunensis_fna.list"
        self.reference = "/mnt/d/breakpoints/assembly/simulation/assembly_test/Staphylococcus_lugdunensis/GCF_008728755.1/GCF_008728755.1_ASM872875v1_genomic.fna"
        self.sim_dir = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/sim_s_lugdunensis/"
        self.res_dir = "/home/wangshuai/assembly_result/s_lugdunensis/"

        self.result_dir = self.res_dir + "/our/"
        self.spades_dir = self.res_dir + "/spades/"
        self.fastas_dir = self.sim_dir + "/fasta/"
        self.fqdir = self.sim_dir + "/fastq/"
        self.truth_file = self.sim_dir + "/fastq/truth_file.txt"

        os.system(f"mkdir {self.sim_dir}")
        os.system(f"mkdir {self.fastas_dir}")
        os.system(f"mkdir {self.fqdir}")
        os.system(f"mkdir {self.res_dir}")
        os.system(f"mkdir {self.result_dir}")
        os.system(f"mkdir {self.spades_dir}")

       
        self.run_script = self.sim_dir + "/run.sh"
        self.spades_script = self.sim_dir + "/run_spades.sh"
        self.tool = "/mnt/d/breakpoints/assembly/scripts/workflow.sh"
        self.spades = "/mnt/d/breakpoints/assembly/scripts/spades.sh"
        self.ass_script = "/mnt/d/breakpoints/assembly/scripts/measure.py"


        

    def select_genomes(self):
        select_num = 0
        os.system(f"rm {self.fastas_dir}/*")
        scaffold_ID_dict = {}

        for line in open(self.original_genome_list):
            num = np.random.randint(10)
            if num < 2:
                origin_ref = line.strip()
                if not os.path.getsize(origin_ref):
                    continue
                fasta_sequences = SeqIO.parse(open(origin_ref),'fasta')       

                chrom_num = 0
                for record in fasta_sequences:
                    if chrom_num == 0:
                        scaffold_ID = record.id
                    chrom_num += 1

                if chrom_num < 3 and scaffold_ID not in scaffold_ID_dict: # make sure only one chromosome
                    # os.system(f"cp {origin_ref} {self.fastas_dir}")
                    os.system(f"head -n 1 {origin_ref} > /{self.fastas_dir}/{scaffold_ID}.fasta")
                    os.system(f"samtools faidx {origin_ref} {scaffold_ID}:1-100000 |grep -v \> >> /{self.fastas_dir}/{scaffold_ID}.fasta")
                    # os.system(f"samtools faidx {origin_ref} {scaffold_ID} |grep -v \> >> /{self.fastas_dir}/{scaffold_ID}.fasta")
                    select_num += 1
                    scaffold_ID_dict[scaffold_ID] = 1

            if select_num == self.chosen_genome_num:
                break
        print ("genome selection is done.")

    def generate_fastq(self, ID, genome):
        fq = f'wgsim -N 1000000 -e 0 -r 0 -R 0 -1 150 -2 150 {genome} {self.fqdir}/{ID}.1.fq {self.fqdir}/{ID}.2.fq'
        os.system(fq)

    def simulate_genomes(self):
        tru = open(self.truth_file, 'w')
        files = os.listdir(self.fastas_dir)
        for genome in files:
            genome = self.fastas_dir + genome
            print (genome)
            f = open(genome, 'r')
            # f = open(genome, encoding= 'unicode_escape')
            line = f.readline()
            mat = re.search("strain (.*?) chromosome", line)
            if mat:
                ID= mat.group(1).strip()
                ID = ID.replace("(", "_")
                ID = ID.replace(")", "_")
                ID = ID.replace("-", "_")
                ID = ID.replace(" ", "_")
                f.close()
                self.generate_fastq(ID, genome)
                print (ID, genome, file = tru)
            # print (ID)
            f.close()
        tru.close()

    def generate_running(self):
        f = open(self.run_script, 'w')
        for line in open(self.truth_file):
            array = line.strip().split()
            ID = array[0]
            fq1 = self.fqdir + ID + ".1.fq"
            fq2 = self.fqdir + ID + ".2.fq"
            genome = array[1]
            order = f"bash {self.tool} {self.reference} {fq1} {fq2} {ID} {self.result_dir} {genome}"
            # order = f"/usr/bin/time -v -o {self.result_dir}/{ID}.time bash {self.tool} {self.reference} {fq1} {fq2} {ID} {self.result_dir} {genome}"
            print (order, file = f)
        f.close()

    def generate_spades(self):
        f = open(self.spades_script, 'w')
        for line in open(self.truth_file):
            array = line.strip().split()
            ID = array[0]
            fq1 = self.fqdir + ID + ".1.fq"
            fq2 = self.fqdir + ID + ".2.fq"
            genome = array[1]
            order = f"bash {self.spades} {ID} {fq1} {fq2} {self.spades_dir} {genome}"
            # order = f"/usr/bin/time -v -o {self.spades_dir}/{ID}.time bash {self.spades} {ID} {fq1} {fq2} {self.spades_dir} {genome}"

            print (order, file = f)
        f.close()




if __name__ == "__main__":
    sim = Simulation()
    sim.select_genomes()


    sim.simulate_genomes()
    sim.generate_running()
    sim.generate_spades()




    