import os
import shutil
import os, re



def download(species, raw_output):
    command = f"""
    ncbi-genome-download --formats fasta --genera "{species}" bacteria --assembly-levels complete -o {raw_output}
    gzip -d {raw_output}/refseq/bacteria/*/*fna.gz
    ls {raw_output}/refseq/bacteria/*/*fna >{raw_output}/fna.list
    """
    print (command)
    os.system(command)


def copy2database(raw_output, database_dir):
    with open(f"{raw_output}/fna.list", 'r') as file_list:
        num = 0
        for file_path in file_list:
            #print (file_path)
            if not os.path.isfile(file_path.strip()):
                print("Skipping file: {}".format(file_path))
                continue
            num += 1
            #if num > 100:
            #    continue
            with open(file_path.strip(), 'r') as f:
                line_count = 0
                chr_count = 0
                for line in f:
                    if line.startswith(">"):
                        chr_count += 1
                    if line_count == 0:
                        #print (line)
                        ID = line[1:].split()[0]
                    line_count += 1
                if chr_count <= 5:
                    shutil.copy2(file_path.strip(), f"{database_dir}/{ID}.fasta")
    #            else:
    #               print("Skipping file: {}".format(file_path))

if __name__ == "__main__":
    species = "Serratia marcescens"
    species_name = species.replace(" ", "_")
    print (species_name)
    raw_output = f"/mnt/d/breakpoints/downloads/{species_name}/"
    database_dir = f"/mnt/d/breakpoints/assembly/sim/database/{species_name}/"


    if not os.path.isdir(raw_output):
        os.system("mkdir %s"%(raw_output))
    if not os.path.isdir(database_dir):
        os.system("mkdir %s"%(database_dir))

    download(species, raw_output)
    copy2database(raw_output, database_dir)