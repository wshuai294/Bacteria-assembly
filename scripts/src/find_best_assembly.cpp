#include <iostream>
#include <fstream>
#include <cmath> 
#include <cstdlib> 
#include <string>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <sys/stat.h> 
#include <pthread.h>
#include <thread>
#include <sstream>
#include <map>
#include <mutex>
#include <queue>


using namespace std;

const short c = 300;
const char coder_num = 3;
const unsigned char least_depth = 127;
const unsigned int window_size = 500;
const unsigned int window_min_hit = 450;

unsigned char low_depth;
int k;
char *kmer_count_table;
long array_size;
float record_match_rate[100000]; // at most 100,000 assemblies

int thread_num;
std::mutex mtx;  


class Encode{
    public:
    char coder [1000];
    char comple [256];
    int base [32];
    short choose_coder[100];
    int k;
    // short *choose_coder; 

    
    void generate_coder(void);
    void generate_complement(void);
    void generate_base(int k);
    void random_coder(int k);
    void constructer(int given_k);

};

void Encode::constructer(int given_k){
    k = given_k;
    this->generate_coder();
    this->generate_complement();
    this->generate_base(given_k);
    this->random_coder(given_k);
    cout <<"Encoder is established."<<endl;
}

void Encode::generate_coder(void){
    // A:65 97 T:116 84 C:99 67 G: 103 71
    // static char coder [1000];
    for (int j = 0; j < 1000; j++){
        coder[j] = 5;
    }
    coder[65] = 1;
    coder[97] = 1;

    coder[84] = 1;
    coder[116] = 1;

    coder[67] = 0;
    coder[99] = 0;

    coder[71] = 0;
    coder[103] = 0;

    coder[65+c] = 1;
    coder[97+c] = 1;

    coder[84+c] = 0;
    coder[116+c] = 0;

    coder[67+c] = 1;
    coder[99+c] = 1;

    coder[71+c] = 0;
    coder[103+c] = 0;

    coder[65+2*c] = 1;
    coder[97+2*c] = 1;

    coder[84+2*c] = 0;
    coder[116+2*c] = 0;

    coder[67+2*c] = 0;
    coder[99+2*c] = 0;

    coder[71+2*c] = 1;
    coder[103+2*c] = 1;

    // return coder;    
}

void Encode::generate_complement(void){
    // static char comple [256];
    for (int j = 0; j < 256; j++){
        comple[j] = 0;
    }
    comple[65] = 84;
    comple[97] = 84;
    comple[116] = 65;
    comple[84] = 65;
    comple[99] = 71;
    comple[67] = 71;
    comple[103] = 67;
    comple[71] = 67;
    // return comple;   
}

void Encode::generate_base(int k){
    
    for (int i = 0; i<k; i++){       
        base[i] = pow(2, k-i-1);
    }
    // return base;    
}

void Encode::random_coder(int k){
    short permu[18] = {0,1,2,0,2,1,1,2,0,1,0,2,2,0,1,2,1,0};
    
    int r;
    unsigned seed;
    seed = time(0);
    srand(seed);
    for (int j = 0; j < k; j++){
        r = rand() % 6;
        // cout << r << endl;
        for (int i = 0; i < coder_num; i++){
            choose_coder[j*3+i] = permu[r*3+i];
            // cout << choose_coder[j*3+i]<<"#"<<permu[r*3+i] << endl;
        }
    }
    // return choose_coder;    
}


string get_read_ID(string reads_seq){
    string delimiter = "/";
    string read_name_forward = reads_seq.substr(0, reads_seq.find(delimiter));
    delimiter = " ";
    read_name_forward = read_name_forward.substr(0, read_name_forward.find(delimiter));
    delimiter = "\t";
    read_name_forward = read_name_forward.substr(0, read_name_forward.find(delimiter));
    return read_name_forward;
}


class Fasta{
    public:
    string ref_seq;
    string line_seq, chr_name;
    string get_ref_seq(string fasta_file);
    string get_ref_name(string fasta_file);
};

string Fasta::get_ref_seq(string fasta_file){
    ifstream fa_file;
    fa_file.open(fasta_file, ios::in);
    // string line_seq;
    
    while (getline(fa_file,line_seq)){
        if (line_seq[0] == '>'){
            ref_seq = "\0";
        }
        else{
            ref_seq += line_seq;           
        }            
    }
    return ref_seq;
}

string Fasta::get_ref_name(string fasta_file){
    ifstream fa_file;
    fa_file.open(fasta_file, ios::in);
    
    while (getline(fa_file,line_seq)){
        if (line_seq[0] == '>'){
            chr_name = get_read_ID(line_seq).substr(1);   
            break; 
        }
    }
    return chr_name;
}

class Select_ref{
    public:
    
    float get_fitness(string fasta_file, Encode encoder);
    string check_each_genome(string genome_list_file, Encode encoder, string match_rate_file);
    void parallele_each_genome(string genome_list_file, Encode encoder, int start_g, int end_g);
    int assign_parallele(string genome_list_file);
    void output_match_rate(string genome_list_file, string match_rate_file);
};

float Select_ref::get_fitness(string fasta_file, Encode encoder){
    Fasta fasta;
    string ref_seq = fasta.get_ref_seq(fasta_file);
    double match_base_num = 0;

    int ref_len;
    char n;
    int m;
    int e;
    unsigned int kmer_index, comple_kmer_index, real_index;
    long extract_ref_len = 0;
    long slide_ref_len = 0;
    time_t t0 = time(0);
    int covert_num, comple_num;
    short convert_ref[300];
    short complemented_ref[300];

    // cout <<"check fitness..."<<endl;
    ref_len= ref_seq.length();
    // cout <<"genome len is "<<ref_len<<endl;
 

    int *ref_int = new int[ref_len];
    int *ref_comple_int = new int[ref_len];
    for (int j = 0; j < ref_len; j++){
        ref_int[j] = (int)ref_seq[j];
        ref_comple_int[j] = encoder.comple[ref_int[j]];
    }
    for (int j = 0; j < ref_len-encoder.k+1; j++){
        int max_dp = 0;
        for (int i = 0; i < 3; i++){
            kmer_index = 0;
            comple_kmer_index = 0;
            bool all_valid = true;
            for (int z = 0; z < encoder.k; z++){
                m = encoder.coder[c*encoder.choose_coder[z*3+i]+ref_int[j+z]];
                if (m == 5){
                    all_valid = false;
                    break;
                }
                kmer_index += m*encoder.base[z]; 
                comple_kmer_index += encoder.coder[c*encoder.choose_coder[(encoder.k-1-z)*3+i]+ref_comple_int[j+z]]*encoder.base[(k-1-z)];  
            }
            if (kmer_index > comple_kmer_index){ //use a smaller index
                real_index = comple_kmer_index;
            }   
            else{
                real_index = kmer_index;
            }
            if (all_valid == false){
                real_index = 0;
            }
            int dp = kmer_count_table[real_index];
            if (max_dp < dp){
                max_dp = dp;
            }
        }
        if (max_dp > low_depth){
            // cout <<"###" <<endl;
            match_base_num += 1;
        }
    }        
    delete [] ref_int;
    delete [] ref_comple_int;

    float match_rate = match_base_num/ref_len;
    // cout << "match rate is\t"<<match_rate<<endl;
    return match_rate;
}

string Select_ref::check_each_genome(string genome_list_file, Encode encoder, string match_rate_file){
    ifstream list_file;
    list_file.open(genome_list_file, ios::in);

    ofstream out_file;
    out_file.open(match_rate_file, ios::out);

    string each_genome, select_genome;
    float match_rate;
    float max_match_rate = 0;
    while (getline(list_file, each_genome)){
        match_rate = this-> get_fitness(each_genome, encoder);
        cout << each_genome<<" match rate is "<<match_rate<<endl;
        out_file << each_genome<<","<<match_rate<<endl;
        if (match_rate >= max_match_rate){
            max_match_rate = match_rate;
            select_genome = each_genome;
        }
    }
    cout << "select genome is "<<select_genome<<endl;
    list_file.close();
    out_file.close();

    return select_genome;
}

int Select_ref::assign_parallele(string genome_list_file){
    ifstream list_file;
    list_file.open(genome_list_file, ios::in);
    string each_genome;
    int genome_num = 0;
    while (getline(list_file, each_genome)){
        genome_num += 1;
    }
    list_file.close();

    int each_thread_g_num = genome_num/thread_num + 1;
    cout<<"There is "<<genome_num<<" genomes."<<endl;
    cout<<"Each thread process "<<each_thread_g_num<<" genomes"<<endl;
    return each_thread_g_num;

}

void Select_ref::parallele_each_genome(string genome_list_file, Encode encoder, int start_g, int end_g){
    ifstream list_file;
    list_file.open(genome_list_file, ios::in);

    string each_genome;
    float match_rate;
    int genome_index = 0;
    while (getline(list_file, each_genome)){
        if (genome_index >= start_g & genome_index < end_g){
            match_rate = this-> get_fitness(each_genome, encoder);
            record_match_rate[genome_index] = match_rate;
        }
        genome_index += 1;
    }

    list_file.close();
}

void Select_ref::output_match_rate(string genome_list_file, string match_rate_file){
    ifstream list_file;
    list_file.open(genome_list_file, ios::in);
    ofstream out_file;
    out_file.open(match_rate_file, ios::out);
    int genome_index = 0;
    string each_genome;

    while (getline(list_file, each_genome)){
        float match_rate = record_match_rate[genome_index];
        out_file << each_genome<<","<<match_rate<<endl;
        genome_index += 1;
    }

    list_file.close();
    out_file.close();

}

void read_fastq(string fastq_file, Encode encoder, int down_sam_ratio, long start, long end)
{
    time_t t0 = time(0);
    ifstream fq_file; 
    fq_file.open(fastq_file);

    long pos = 0;
    for (long i = start; i>0; i--){
        fq_file.seekg(i, ios::beg);
        char j;
        fq_file.get(j);
        if (j == '@'){ //only read name has this symbol.
            pos = i;
            break;
        }       
    }
    fq_file.seekg(pos, ios::beg);
    long add_size = start;

    string reads_seq;
    int reads_int [300];
    int reads_comple_int [300];

    unsigned int lines = 0;
    int converted_reads [900];
    int complemented_reads [900];
    int m;
    int n;
    unsigned int kmer_index, comple_kmer_index, real_index, b;   
    int r ;
    short read_len = 0;

    while (getline(fq_file,reads_seq))
    {
        if (add_size>=end){
            break;
        }
        add_size += reads_seq.length();

        if (lines % 4 == 1){
            time_t t1 = time(0);
            if (lines % 10000000 == 10000000-1){
                cout <<lines<<"reads\t" << t1-t0 <<endl;
            }
            read_len = reads_seq.length();//cal read length
            r = rand() % 100 ;
            // cout <<r << "r"<<lines<<endl;
            
            if (r < down_sam_ratio){
                for (int j = 0; j < read_len; j++){
                    reads_int[j] = (int)reads_seq[j];
                    reads_comple_int[j] = encoder.comple[reads_int[j]];
                }
                
                for (int j = 0; j < read_len-encoder.k+1; j++){
                    
                    for (int i = 0; i < 3; i++){
                        // cout <<j<<"hi"<<i<<endl;
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        
                        bool all_valid = true;
                        // cout <<j<<"\t"<<i<<"\t"<<endl;
                        for (int z = 0; z<encoder.k; z++){
                            m = encoder.coder[c*encoder.choose_coder[z*3+i]+reads_int[j+z]]; // choose_coder[z*3+i] indicate which coder
                            n = encoder.coder[c*encoder.choose_coder[(k-1-z)*3+i]+reads_comple_int[j+z]];
                            if (m == 5){
                                all_valid = false;
                                // cout << "N !" << endl;
                                break;
                            }
                            kmer_index += m*encoder.base[z]; 
                            comple_kmer_index += n*encoder.base[(encoder.k-1-z)];
                              
                        }
                        
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }   
                        else{
                            real_index = kmer_index;
                        }
                        if ((int)kmer_count_table[real_index] < least_depth & all_valid == true ){
                            kmer_count_table[real_index] += 1;
                            // cout << (int)kmer_count_table[real_index] << "\t" << endl;
                        }  
                    }
                }
            }         
        }
        lines++;
    }
    fq_file.close();
    // return kmer_count_table;

}


short * saved_random_coder(string index_name){
    ifstream index_file;
    static short read_choose_coder[100];
    long i = 0;
    unsigned int real_index;
    index_file.open(index_name, ios::in | ios::binary);   
    while (!index_file.eof()){
        index_file.read(reinterpret_cast<char*>(&real_index), sizeof(unsigned int));
        if (i < 100){
            read_choose_coder[i] = real_index;
        }
        else{
            break;
        }
        i += 1;  
    }
    index_file.close(); 
    return read_choose_coder;
}

int cal_sam_ratio(string fq1, long down_sampling_size){
    int down_sam_ratio = 0;
    int read_len = 0;
    long i = 0;
    long sample_size = 0;
    string reads_seq;

    ifstream fq_file; 
    fq_file.open(fq1);
    while (getline(fq_file, reads_seq)){
        if (i == 1){
            read_len = reads_seq.length();
            cout <<"The first read length is "<<read_len<<endl;
        }
        if (i % 4 == 1){
            sample_size += reads_seq.length();
        }
        i += 1;
    }
    fq_file.close();
    sample_size = sample_size * 2; //pair
    down_sam_ratio = 100*down_sampling_size/sample_size;
    cout<<i<<"\t"<<sample_size<<" Down-sampling ratio is "<<down_sam_ratio<< "%." <<endl;
    return down_sam_ratio;
}

long file_size(string filename)  
{  
    struct stat statbuf;  
    stat(filename.c_str(),&statbuf);  
    long size=statbuf.st_size;   
    return size;  
}  

int main( int argc, char *argv[])
{
    unsigned seed;
    // seed = time(0);
    seed = 1;
    srand(seed);
    time_t now1 = time(0);

    string fq1 = argv[1];
    string fq2 = argv[2];  
    string fasta_list = argv[3]; // a fasta file list, where each fasta is a single chromosome
    k = stod(argv[4]); // kmer length
    low_depth = stod(argv[5]); // lower than this value will not be regarded as k-mer hit
    thread_num = stod(argv[6]);
    string match_rate_file = argv[7]; // output file, the match rate of each fasta
    array_size = pow(2, k);
    kmer_count_table = new char[array_size];
    int down_sam_ratio = stod(argv[8]); // percentage (0-100), randomly select reads with this probability, 30

    
    Encode encoder;
    encoder.constructer(k);

    string index_name = "no use"; //fasta_file + ".k" + to_string(k) + ".index.dat";
    
    
    long start = 0;
    long end = 0;
    long size = file_size(fq1);
    long each_size = size/thread_num;
    
    cout <<"size\t"<<size<<endl;
    
    
    std::vector<std::thread>threads;

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(read_fastq, fq1, encoder, down_sam_ratio, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();
  
    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(read_fastq, fq2, encoder, down_sam_ratio, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();
    time_t now2 = time(0);
    cout << "reads are loaded. "<< now2 - now1<<endl;



    Select_ref select;
    // string select_genome = select.check_each_genome(fasta_list, encoder, match_rate_file); //"/mnt/d/breakpoints/assembly/sim/sim_s_lugdunensis/fasta/s_lugdunensis.database"
    int each_thread_g_num = select.assign_parallele(fasta_list);
    int start_g = 0;
    int end_g = 0;
    for (int i=0; i<thread_num; i++){
        start_g = i*each_thread_g_num;
        end_g = (i+1)*each_thread_g_num;
        threads.push_back(thread(&Select_ref::parallele_each_genome, select, fasta_list, encoder, start_g, end_g));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();
    select.output_match_rate(fasta_list, match_rate_file);

    return 0;
}
