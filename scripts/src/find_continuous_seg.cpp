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
#include <list>
#include <algorithm>


using namespace std;

const short c = 300;
const char coder_num = 3;
int k = 26;
int *kmer_position_table;
char *kmer_pos_num_table;


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
    // unsigned seed;
    // seed = time(0);
    // srand(seed);
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
    string line_seq, ref_name;
    void get_ref_seq(string fasta_file);
    // string get_ref_name(string fasta_file);
};

void Fasta::get_ref_seq(string fasta_file){
    ifstream fa_file;
    fa_file.open(fasta_file, ios::in);
    // string line_seq;

    ref_name = "undefined";
    
    while (getline(fa_file,line_seq)){
        if (line_seq[0] == '>'){
            if (ref_name == "undefined"){
                ref_seq = "\0";
                ref_name = get_read_ID(line_seq).substr(1); 
            }
            else{
                break;  // only consider one contig
            }
            
        }
        else{
            ref_seq += line_seq;           
        }            
    }
    fa_file.close();
    // return ref_seq;
}

void read_ref(string ref_seq, Encode encoder)
{
    int ref_len;
    char n;
    int m;
    int e;
    int dp;
    unsigned int kmer_index, comple_kmer_index, real_index;
    long extract_ref_len = 0;
    long slide_ref_len = 0;
    // string chr_name ;
    time_t t0 = time(0);
    int covert_num, comple_num;


    ref_len= ref_seq.length();

    int *ref_int = new int[ref_len];
    int *ref_comple_int = new int[ref_len];
    for (int j = 0; j < ref_len; j++){
        ref_int[j] = (int)ref_seq[j];
        ref_comple_int[j] = encoder.comple[ref_int[j]];
    }
    for (int j = 0; j < ref_len-encoder.k+1; j++){
        int max_dp = 0;
        // int min_dp = 0;
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
                comple_kmer_index += encoder.coder[c*encoder.choose_coder[(encoder.k-1-z)*3+i]+ref_comple_int[j+z]]*encoder.base[(encoder.k-1-z)];  
            }
            if (kmer_index > comple_kmer_index){ //use a smaller index
                real_index = comple_kmer_index;
            }   
            else{
                real_index = kmer_index;
            }
            if (all_valid == false){
                real_index = 0;
                dp = 0;
            }
            else{
                // if (kmer_position_table[real_index] == 0){
                //     kmer_position_table[real_index] = j;
                // }
                // else{
                //     int r = rand() % 100 ;
                //     if (r < 50){
                //         kmer_position_table[real_index] = j;
                //     }
                // }
                kmer_position_table[real_index] = j;
                if (kmer_pos_num_table[real_index] < 10){
                    kmer_pos_num_table[real_index] += 1;
                }
            }
        }
    }        
    delete [] ref_int;
    delete [] ref_comple_int;
}

long get_fq_start(ifstream& fq_file, long start){ // find the read name line of a fastq file
    long pos = 0;
    bool flag = false;
    bool done = false;
    int x = 0;
    for (long i = start; i>0; i--){

        for (long j = i; j < i + 1000; j ++){
            fq_file.seekg(j, ios::beg);
            char chr1, chr2;
            fq_file.get(chr1);
            fq_file.get(chr2);

            if (chr1 == '\n' and chr2 == '+'){ //third field
                flag = true;
                x = 0;
            }  
            if (flag) {
                if (chr1 == '\n'){
                    x += 1;
                }
                if (chr1 == '\n' and chr2 == '@' and x == 3){
                    pos = j + 1;
                    done = true;
                    break;
                }
            }  
            else{
                if (chr1 == '\n'){
                    x += 1;
                }
            }
            if (x == 3){
                break;
            } 
        }
        // cout <<  i << "\t" << pos <<endl; 
        if (done){
            break;
        }
    }
    return pos;
}


void read_fastq(string fastq_file, Encode encoder, int down_sam_ratio, string out_file, long start, long end)
{
    time_t t0 = time(0);
    ifstream fq_file; 
    fq_file.open(fastq_file);

    ofstream out; 
    out.open(out_file);

    string reads_seq;
    int reads_int [500];
    int reads_comple_int [500];

    unsigned int lines = 0;
    int m;
    int n;
    unsigned int kmer_index, comple_kmer_index, real_index, b;   
    int r ;
    short read_len = 0;

    long pos = get_fq_start(fq_file, start);
    fq_file.seekg(pos, ios::beg);
    long add_size = start;

    while (getline(fq_file,reads_seq)){


        if (add_size>=end){
            break;
        }
        add_size += reads_seq.length() + 1;
        if (lines == 0){
            cout << reads_seq << " read name " <<endl;
        }

        if (lines % 4 == 1){
            time_t t1 = time(0);
            if (lines % 50000 == 1){
                cout <<lines<<" reads " << t1-t0 <<endl;
            }
            read_len = reads_seq.length();//cal read length
            r = rand() % 100 ;
            // cout <<r << "r"<<lines<<endl;
            if (read_len > 499){
                cout <<"read is too long!"<<endl;
            }
            
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
                        if (all_valid == true ){
                            
                            // cout << lines<< " "<< j<< " " <<i << " "  <<(int)kmer_position_table[real_index] << "\t" << (int)kmer_pos_num_table[real_index] << endl;
                            out << lines<< " "<< j<< " " <<i << " "  <<(int)kmer_position_table[real_index] << "\t" << (int)kmer_pos_num_table[real_index] << endl;
                        }  
                    }
                }
            }         
        }
        lines++;
    }
    fq_file.close();
}

long file_size(string filename){  
    struct stat statbuf;  
    stat(filename.c_str(),&statbuf);  
    long size=statbuf.st_size;   
    return size;  
}  

int main( int argc, char *argv[])
{
    unsigned seed;
    seed = 1;
    srand(seed);



    string reference = argv[1];
    string fastq_file = argv[2];
    string fastq_file_2 = argv[3];
    string out_file = argv[4];
    int down_sam_ratio = stod(argv[5]);
    k = stod(argv[6]);
    int thread_num = stod(argv[7]);

    // int down_sam_ratio = 50;
    Encode encoder;
    encoder.constructer(k);



    Fasta fasta;
    fasta.get_ref_seq(reference);
    int array_size = pow(2, k);
    kmer_position_table = new int[array_size];
    kmer_pos_num_table = new char[array_size]; 
    memset(kmer_position_table, 0, sizeof(kmer_position_table));
    memset(kmer_pos_num_table, 0, sizeof(kmer_pos_num_table));

    read_ref(fasta.ref_seq, encoder);
    cout <<"ref is loaded."<<endl;
    
    // read_fastq(fastq_file, encoder, down_sam_ratio, out_file);


    long size = file_size(fastq_file);
    long each_size = size/thread_num;  
    long start = 0;
    long end = 0;
    cout <<"first read file size:\t"<<size<<endl;
    std::vector<std::thread>threads;
    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        string part_outfile = out_file + "_" + to_string(start)+ "_" + to_string(end) + ".fq1.part.txt";
        threads.push_back(thread(read_fastq, fastq_file, encoder, down_sam_ratio, part_outfile, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();

    // size = file_size(fastq_file_2);
    // each_size = size/thread_num;  
    // start = 0;
    // end = 0;
    // cout <<"second read file size:\t"<<size<<endl;

    // for (int i=0; i<thread_num; i++){
    //     start = i*each_size;
    //     end = (i+1)*each_size;
    //     if (i == thread_num-1){
    //         end = size;
    //     }
    //     string part_outfile = out_file + "_" + to_string(start)+ "_" + to_string(end) + ".fq2.part.txt";
    //     threads.push_back(thread(read_fastq, fastq_file_2, encoder, down_sam_ratio, part_outfile, start, end));
    // }
	// for (auto&th : threads)
	// 	th.join();
    // threads.clear();




    delete [] kmer_position_table;
    delete [] kmer_pos_num_table;
}