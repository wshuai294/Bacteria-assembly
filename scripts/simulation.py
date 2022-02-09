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

class Simulation():


    def __init__(self):

        self.depth = 50
        self.reads_len = 150
        self.chosen_genome_num = 10

        self.original_genome_list = "/mnt/d/breakpoints/assembly/simulation/ecoli/fna.list"
        self.fastas_dir = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/fasta/"
        self.fqdir = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/fastq/"
        self.truth_file = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/fastq/truth_file.txt"
        self.reference = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/ecoli_ref.fna"
        self.tool = "/mnt/d/breakpoints/assembly/scripts/workflow.sh"
        self.result_dir = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/result/"
        self.run_script = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/run.sh"
        self.spades_dir = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/spades/"
        self.spades_script = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/run_spades.sh"
        self.spades = "/mnt/d/breakpoints/assembly/scripts/spades.sh"
        self.ass_script = "/mnt/d/breakpoints/assembly/scripts/measure.py"


        

    def select_genomes(self):
        select_num = 0
        for line in open(self.original_genome_list):
            num = np.random.randint(10)
            if num < 1:
                origin_ref = line.strip()
                fasta_sequences = SeqIO.parse(open(origin_ref),'fasta')       

                chrom_num = 0
                for record in fasta_sequences:
                    chrom_num += 1
                if chrom_num == 1:
                    os.system(f"cp {origin_ref} {self.fastas_dir}")
                    select_num += 1
            if select_num == self.chosen_genome_num:
                break

    def generate_fastq(self, ID, genome):
        fq = f'wgsim -e 0 -N 500000 -e 0 -r 0 -R 0 -1 150 -2 150 {genome} {self.fqdir}/{ID}.1.fq {self.fqdir}/{ID}.2.fq'
        os.system(fq)

    def simulate_genomes(self):
        tru = open(self.truth_file, 'w')
        files = os.listdir(self.fastas_dir)
        for genome in files:
            genome = self.fastas_dir + genome
            f = open(genome)
            line = f.readline()
            mat = re.search("strain (.*?) chromosome", line)
            ID= mat.group(1).strip()
            f.close()
            self.generate_fastq(ID, genome)
            print (ID, genome, file = tru)
            # print (ID)
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

            print (order, file = f)
        f.close()





sim = Simulation()
# sim.select_genomes()
# sim.simulate_genomes()
# sim.generate_running()
sim.generate_spades()




    