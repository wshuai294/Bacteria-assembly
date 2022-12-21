"""
compute the CPU time (h) and max RAM (G)
on the log file generated by /usr/bin/time -v
"""

import re 
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from datetime import datetime
import os

class Com_Res():

    def __init__(self, time_file):
        self.time_file = time_file
    
    def extract_time(self): #log file obtained by /usr/bin/time -v
         #if no time available
        user_time, sys_time = None, None
        for line in open(self.time_file):
            time_re = re.search('User time \(seconds\):(.*?)$', line)
            if time_re:
                user_time =  time_re.group(1).strip()

            time_sys = re.search('System time \(seconds\):(.*?)$', line)
            if time_sys:
                sys_time = time_sys.group(1).strip()
        if user_time == None:
            print ("fail to extract time.")
        # print (user_time, sys_time)
        all_time = float(user_time) + float(sys_time)
        final_time = round(all_time, 1)
        return final_time

    def extract_mem(self):
        used_mem = 0 #if no time available
        for line in open(self.time_file):
            time_re = re.search('Maximum resident set size \(kbytes\):(.*?)$', line)
            if time_re:
                used_mem =  time_re.group(1).strip()
        final_mem = round(float(used_mem)/1000000, 2)
        return final_mem

class Sample():
    def __init__(self, ID, result_dir, true_fasta, method):
        print ("for sample", ID)
        self.ID = ID
        self.time_file = result_dir + "/" + ID + ".time"
        if not os.path.isfile(self.time_file):
            self.time_file = "/".join(result_dir.split("/")[:-1]) + "/" + ID + ".time"
        # self.N50_file = result_dir + "/" + ID + ".assessment"
        self.result_dir = result_dir
        self.N50 = None
        self.NGA50 = None
        self.base_error = None
        self.misassemblies = None
        self.cpu_time = None
        self.max_PAM = None
        self.method = method
        self.true_fasta = true_fasta
        self.result_fasta = None
        self.quast_dir = self.result_dir + "/quast_" + ID
        self.get_time()
        # self.get_N50()
        self.run_quast()

    def get_result(self):
        if self.method == "Our":
            self.result_fasta = self.result_dir + "/" + self.ID + ".contigs.consensus.fasta"
        else:
            self.result_fasta = self.result_dir + "/" + "/contigs.fasta"
    
    def run_quast(self):
        self.get_result()
        command = f"~/softwares/quast-5.2.0/quast.py {self.result_fasta} -r {self.true_fasta} -o {self.quast_dir}"
        # os.system(command)
        quast = Quast(self.quast_dir)
        self.NGA50 = quast.NGA50
        self.base_error = quast.base_error
        self.misassemblies = quast.misassemblies

        # print (quast.NGA50, quast.mismatch_rate, quast.indel_rate, quast.misassemblies)

    
    def get_time(self):
        com = Com_Res(self.time_file)
        self.cpu_time = com.extract_time()
        self.max_PAM = com.extract_mem()
    
    def get_N50(self):
        f = open(self.N50_file, 'r')
        line = f.readline()
        time_sys = re.search('Minimap_N50: (.*?),size:', line)
        if time_sys:
            N50 = time_sys.group(1).strip()
        else:
            print ("No pattern to extract N50!")
        self.N50 = int(N50)

class Benchmark():

    def __init__(self, run_script):
        self.sample_list = []
        self.our_dir = None
        self.true_fasta_dict = {}
        self.get_para(run_script)
        self.spades_dir = self.our_dir.replace("our", "spades")
        self.data = []
        

    def get_para(self, run_script):
        f = open(run_script, 'r')
        for line in f:
            array = line.strip().split()
            # print (array)
            ID = array[9]
            self.sample_list.append(ID)
            self.our_dir = array[10]
            self.true_fasta_dict[ID] = array[11]

        f.close()
    
    def main(self):
        for ID in self.sample_list:
            if ID == "VISLISI_33":
                continue
            print (ID)
            sample = Sample(ID, self.our_dir, self.true_fasta_dict[ID], "Our")
            self.data.append([sample.ID, sample.NGA50, sample.cpu_time, sample.max_PAM, "Our", sample.misassemblies, sample.base_error])
            sample = Sample(ID, self.spades_dir + "/" + ID, self.true_fasta_dict[ID], "Spades")
            self.data.append([sample.ID, sample.NGA50, sample.cpu_time, sample.max_PAM, "Spades", sample.misassemblies, sample.base_error])
            # break

        self.df=pd.DataFrame(self.data, columns=['ID', 'NGA50','CPU_Time', 'Peak_RAM', 'Methods', "Misassemblies", "Base_Error(/100k)"])
        self.make_figures()

    
    def make_figures(self):

        fig, axes = plt.subplots(5, 2, gridspec_kw={'width_ratios': [9, 1]}, figsize=(60,18))
        sns.barplot(ax = axes[0][0], x="ID",y='NGA50',hue= 'Methods',data=self.df, log=True)
        sns.boxplot(ax = axes[0][1], x="Methods",y='NGA50',data=self.df, showfliers =False, showmeans=True,meanprops={"marker":"o",
                       "markerfacecolor":"white", 
                       "markeredgecolor":"black",
                      "markersize":"10"})
        # sns.catplot(ax = axes[0][0],x="ID",y='NGA50',hue= 'Methods',data=self.df, log=True)
        sns.barplot(ax = axes[1][0], x="ID",y='CPU_Time',hue= 'Methods',data=self.df)
        sns.boxplot(ax = axes[1][1], x="Methods",y='CPU_Time',data=self.df, showmeans=True,meanprops={"marker":"o",
                       "markerfacecolor":"white", 
                       "markeredgecolor":"black",
                      "markersize":"10"})
        sns.barplot(ax = axes[2][0], x="ID",y='Peak_RAM',hue= 'Methods',data=self.df)
        sns.boxplot(ax = axes[2][1], x="Methods",y='Peak_RAM',data=self.df, showmeans=True,meanprops={"marker":"o",
                       "markerfacecolor":"white", 
                       "markeredgecolor":"black",
                      "markersize":"10"})
        sns.barplot(ax = axes[3][0], x="ID",y='Misassemblies',hue= 'Methods',data=self.df)
        sns.boxplot(ax = axes[3][1], x="Methods",y='Misassemblies',data=self.df, showmeans=True,meanprops={"marker":"o",
                       "markerfacecolor":"white", 
                       "markeredgecolor":"black",
                      "markersize":"10"})
        sns.barplot(ax = axes[4][0], x="ID",y='Base_Error(/100k)',hue= 'Methods',data=self.df)
        sns.boxplot(ax = axes[4][1], x="Methods",y='Base_Error(/100k)',data=self.df, showmeans=True,meanprops={"marker":"o",
                       "markerfacecolor":"white", 
                       "markeredgecolor":"black",
                      "markersize":"10"})

        give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
        # plt.xticks(rotation=30)
        plt.savefig('/home/wangshuai/assembly_result/figures/assembly_comparison_%s.pdf'%(give_time))

    def compar_ave(self):
        # self.df.groupby(['Methods'])['N50'].mean()
        mean_df = self.df.groupby(['Methods'])['NGA50'].mean().reset_index()
        our_mean_n50 = mean_df.at[0, 'NGA50']
        spades_mean_n50 = mean_df.at[1, 'NGA50']
        diff = our_mean_n50 - spades_mean_n50
        ratio = round(diff/spades_mean_n50, 2)
        print (mean_df)
        print (diff,ratio)

class Quast():

    def __init__(self, quast_result_dir):
        self.result_txt = quast_result_dir + "/report.txt" 
        self.NGA50 = None
        self.mismatch_rate = None
        self.indel_rate = None
        self.base_error = None
        self.misassemblies = None
        self.read_report()

        # 
    
    def read_report(self):
        f = open(self.result_txt, "r")
        for line in f:
            array = line.strip().split()
            if len(array) < 2:
                continue
            if array[0] == "#" and len(array) == 3 and array[1] == "misassemblies":
                self.misassemblies = int(array[2])
            elif array[0] == 'NGA50':
                # print (array)
                self.NGA50 = int(array[1])
            elif len(array) == 6 and array[1] == 'mismatches':
                self.mismatch_rate = float(array[5])
            elif len(array) == 6 and array[1] == 'indels':
                self.indel_rate = float(array[5])
            # print (array)
        f.close()
        self.base_error = self.mismatch_rate + self.indel_rate

if __name__ == "__main__":
    run_script = "/mnt/d/breakpoints/assembly/sim/run.sh"
    # run_script = "/mnt/d/breakpoints/assembly/sim/sim_s_lugdunensis/run.sh"
    ben = Benchmark(run_script)
    ben.main()
    ben.compar_ave()


    # quast = Quast("/home/wangshuai/assembly_result/ecoli/spades/104/quast_results/")
    

    # sa = Sample("2_0", "/home/wangshuai/assembly_result/ecoli/our/")
    # print (sa.N50, sa.cpu_time, sa.max_PAM)

    # time_file = sys.argv[1]
    # com = Com_Res(time_file)
    # cpu_time = com.extract_time()
    # max_PAM = com.extract_mem()
    # print ("CPU time: %s h; Maximum RAM: %s G."%(cpu_time, max_PAM))