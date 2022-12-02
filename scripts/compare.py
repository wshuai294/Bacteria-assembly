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
    def __init__(self,ID,result_dir):
        print ("for sample", ID)
        self.ID = ID
        self.time_file = result_dir + "/" + ID + ".time"
        if not os.path.isfile(self.time_file):
            self.time_file = "/".join(result_dir.split("/")[:-1]) + "/" + ID + ".time"
        self.N50_file = result_dir + "/" + ID + ".assessment"
        self.N50 = None
        self.cpu_time = None
        self.max_PAM = None
        self.get_time()
        self.get_N50()

    
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
        f.close()
    
    def main(self):
        for ID in self.sample_list:
            sample = Sample(ID, self.our_dir)
            print (ID)
            self.data.append([sample.ID, sample.N50, sample.cpu_time, sample.max_PAM, "Our"])
            sample = Sample(ID, self.spades_dir + "/" + ID)
            self.data.append([sample.ID, sample.N50, sample.cpu_time, sample.max_PAM, "Spades"])
            # break

        self.df=pd.DataFrame(self.data, columns=['ID', 'N50','CPU_Time', 'Peak_RAM', 'Methods'])
        self.make_figures()

    
    def make_figures(self):

        fig, axes = plt.subplots(2, 2, figsize=(40,10))
        sns.barplot(ax = axes[0][0], x="ID",y='N50',hue= 'Methods',data=self.df, log=True)
        # sns.catplot(ax = axes[0][0],x="ID",y='N50',hue= 'Methods',data=self.df, log=True)
        sns.barplot(ax = axes[0][1], x="ID",y='CPU_Time',hue= 'Methods',data=self.df)
        sns.barplot(ax = axes[1][0], x="ID",y='Peak_RAM',hue= 'Methods',data=self.df)
        # axes[0,1].set_ylim(0,0.05)
     
        #     plt.xticks(rotation=0)
        give_time = datetime.now().strftime("%Y_%m_%d_%H_%M")
        plt.savefig('/home/wangshuai/assembly_result/figures/assembly_comparison_%s.pdf'%(give_time))


if __name__ == "__main__":
    # run_script = "/mnt/d/breakpoints/assembly/sim/run.sh"
    run_script = "/mnt/d/breakpoints/assembly/sim/sim_s_lugdunensis/run.sh"
    ben = Benchmark(run_script)
    ben.main()

    # sa = Sample("2_0", "/home/wangshuai/assembly_result/ecoli/our/")
    # print (sa.N50, sa.cpu_time, sa.max_PAM)

    # time_file = sys.argv[1]
    # com = Com_Res(time_file)
    # cpu_time = com.extract_time()
    # max_PAM = com.extract_mem()
    # print ("CPU time: %s h; Maximum RAM: %s G."%(cpu_time, max_PAM))