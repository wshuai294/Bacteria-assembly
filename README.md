# Bacteria-assembly

## Install
```
git clone git@github.com:wshuai294/Bacteria-assembly.git --depth 1
conda create --name assembly -f environment.yml
conda activate assembly
make
```

## Dependencies
Please install Spades and add it to system path.


## Test
```
cd test/
sh test.sh
```

## Run
Given paired-end reads data and a list file with complete genomes of a specific species, run 
```
python scripts/main.py -h
```

