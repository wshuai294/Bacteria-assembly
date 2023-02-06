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

// const unsigned char low_depth = 10;
// const int k = 26;
// long array_size = pow(2, k);
// char *kmer_count_table = new char[array_size];

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

class Unmap{
    public:
        int save_map_region[1000]; // save segment intervals [start end start end ...]
        int save_map_index = 0;
        char *map_kmer_table = new char[array_size];
        string ref_seq;
        void get_unmap(char* kmer_hit_array, int ref_len);
        void get_map_kmer(Encode encoder, string ref_seq);
        void get_unmap_reads(string fastq_file, Encode encoder);
        void output_map_segments(string fasta_file, string ref_seq, string outdir, string ID);

};

void Unmap::get_unmap(char* kmer_hit_array, int ref_len){
    int start = 0;
    int end = 0;
    int good_base_num = 0; // in a windows
    bool good_flag = false;
    int each_base = 0;
    // short record_each_base[window_size];
    queue<int>record_each_base;
    for (int j = 0; j < ref_len; j++){
        if (kmer_hit_array[j] >= low_depth){
            each_base = 1;
        }
        else{
            each_base = 0;
        }

        if (j < window_size){
            record_each_base.push(each_base);
            good_base_num += each_base;
        }
        else{
            good_base_num = good_base_num - record_each_base.front() + each_base;
            record_each_base.pop();
            record_each_base.push(each_base);
        }
        end += 1;
        // cout << j <<"\t" <<good_base_num<<endl; 
        if (good_base_num >= window_min_hit){
            if (good_flag == false){
                start = end + 1 - window_size;
                if (start < 0){
                    start = 1;
                }
                good_flag = true;
            } 
        }
        else{
            if (good_flag == true){
                // cout << start <<"\t" <<end<<endl;
                if (save_map_index == 0){
                    save_map_region[save_map_index*2] = start;
                    save_map_region[save_map_index*2+1] = end; 
                    save_map_index += 1;
                }
                else{
                    if (start < save_map_region[save_map_index*2-1]){
                        save_map_region[save_map_index*2-1] = end;
                    }
                    else{
                        save_map_region[save_map_index*2] = start;
                        save_map_region[save_map_index*2+1] = end; 
                        save_map_index += 1;
                    }  
                }
            }
            good_flag = false;
        }
    }
    if (good_flag == true){
        if (start < save_map_region[save_map_index*2-1]){
            save_map_region[save_map_index*2-1] = end;
        }
        else{
            save_map_region[save_map_index*2] = start;
            save_map_region[save_map_index*2+1] = end; 
            save_map_index += 1;
            
        }
    }

}

void Unmap::get_map_kmer(Encode encoder, string ref_seq){
    unsigned int kmer_index, comple_kmer_index, real_index;
    int m;
    for (int z = 0; z < save_map_index; z++){
        int start = save_map_region[z*2];
        int end = save_map_region[z*2+1];
        string segment = ref_seq.substr(start, end);
        int seg_len = segment.length();

        // cout <<z <<"map seg\t" <<start <<"\t" <<end<<endl;

        int *ref_int = new int[seg_len];
        int *ref_comple_int = new int[seg_len];
        for (int j = 0; j < seg_len; j++){
            ref_int[j] = (int)ref_seq[j];
            ref_comple_int[j] = encoder.comple[ref_int[j]];
        }

        for (int j = 0; j < seg_len-k+1; j++){
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
                    comple_kmer_index += encoder.coder[c*encoder.choose_coder[(k-1-z)*3+i]+ref_comple_int[j+z]]*encoder.base[(k-1-z)];  
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
                map_kmer_table[real_index] = 1;

            } 

        }
        delete [] ref_int;
        delete [] ref_comple_int;
    }

}

void Unmap::get_unmap_reads(string fastq_file, Encode encoder){
    ifstream fq_file; 
    fq_file.open(fastq_file);
    ofstream out_fq_file;
    out_fq_file.open(fastq_file+".unmap.fq", ios::out);
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
    string read_fq_info;
    bool unmap_flag;
    while (getline(fq_file, reads_seq))
    {
        if (lines % 4 == 1){
            read_len = reads_seq.length();//cal read length
            for (int j = 0; j < read_len; j++){
                reads_int[j] = (int)reads_seq[j];
                reads_comple_int[j] = encoder.comple[reads_int[j]];
            }
            int hit_kmer_num = 0;
            for (int j = 0; j < read_len-encoder.k+1; j++){
                bool hit_flag = false;
                for (int i = 0; i < 3; i++){
                    kmer_index = 0;
                    comple_kmer_index = 0;
                    bool all_valid = true;
                    // cout <<j<<"\t"<<i<<"\t"<<endl;
                    for (int z = 0; z < encoder.k; z++){
                        m = encoder.coder[c*encoder.choose_coder[z*3+i]+reads_int[j+z]]; // choose_coder[z*3+i] indicate which coder
                        n = encoder.coder[c*encoder.choose_coder[(encoder.k-1-z)*3+i]+reads_comple_int[j+z]];
                        if (m == 5){
                            all_valid = false;
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
                    if (all_valid == true & map_kmer_table[real_index] > 0){
                        hit_flag = true;
                    }  
                }
                if (hit_flag){
                    hit_kmer_num += 1;
                }
            } 
            if (hit_kmer_num < 50){ // determine if a read is unmap 
                unmap_flag = true;
            }
            // cout << lines << "\t" << hit_kmer_num <<endl;    
        }
        else{
            if (lines % 4 == 0){
                read_fq_info =  "\0";
                unmap_flag = false;
            }
        }        
        read_fq_info += reads_seq + '\n';
        if (lines % 4 == 3 & unmap_flag){
            // cout << read_fq_info << endl;
            out_fq_file << read_fq_info; 
        }
        lines++;
    }
    fq_file.close();
    out_fq_file.close();

}

void Unmap::output_map_segments(string fasta_file, string ref_seq, string outdir, string ID){
    int start = 0;

    int end = 0;
    ofstream out_fa_file;
    out_fa_file.open(outdir+ "/"+ ID + ".map.fasta", ios::out);
    string map_segment_seq, map_segment_name;
    for (int j = 0; j < save_map_index; j++){
        start = save_map_region[j*2];
        end = save_map_region[j*2+1];
        // cout <<j <<"final\t" <<start <<"\t" <<end<<endl;
        map_segment_seq = ref_seq.substr(start, end-start);
        map_segment_name = ">chr:"+ to_string(start) + "-" + to_string(end);

        out_fa_file << map_segment_name<<endl; 
        out_fa_file << map_segment_seq<<endl;
    }
    out_fa_file.close();
}

char* read_ref(string ref_seq, Encode encoder, string index_name)
{
    // ifstream fa_file;
    // ofstream index_file;
    // ofstream len_file;
    // index_file.open(index_name, ios::out | ios::binary);
    // len_file.open(fasta_file+".genome.len.txt", ios::out);
    // fa_file.open(fasta_file, ios::in);
    // string ref_seq, line_seq;
    // ref_seq = "\0";
    int ref_len;
    char n;
    int m;
    int e;
    unsigned int kmer_index, comple_kmer_index, real_index;
    long extract_ref_len = 0;
    long slide_ref_len = 0;
    // string chr_name ;
    time_t t0 = time(0);
    int covert_num, comple_num;
    short convert_ref[300];
    short complemented_ref[300];

    // save random coder
    // for (int j = 0; j < 100; j++){
    //     index_file.write((char *)(&choose_coder[j]), sizeof(unsigned int));
    // }
    cout <<"Start index ref..."<<endl;

    // while (getline(fa_file,line_seq)){
    //     if (line_seq[0] == '>'){
    //         chr_name = get_read_ID(line_seq).substr(1);
    //         ref_seq = "\0";
    //     }
    //     else{
    //         ref_seq += line_seq;           
    //     }            
    // }
    ref_len= ref_seq.length();
    char *kmer_hit_array = new char[ref_len];
 
    // for (int i = 0; i < 3; i++){
    //     kmer_index = 0; // start kmer
    //     comple_kmer_index = 0;
    // }
    int *ref_int = new int[ref_len];
    int *ref_comple_int = new int[ref_len];
    for (int j = 0; j < ref_len; j++){
        ref_int[j] = (int)ref_seq[j];
        ref_comple_int[j] = encoder.comple[ref_int[j]];
    }
    // index_file.write((char *)(&ref_len), sizeof(unsigned int));
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
            }
            int dp = kmer_count_table[real_index];
            if (max_dp < dp){
                max_dp = dp;
            }
            // cout <<j<<"\t"<< i<<"\t"<<dp<<"\t" <<endl;
            // index_file.write((char *)(&real_index), sizeof(real_index));  
        }
        kmer_hit_array[j] = max_dp; 
    }        
    delete [] ref_int;
    delete [] ref_comple_int;
    delete [] kmer_count_table;
    time_t t1 = time(0);
    // cout << "chr:\t" << chr_name<<"\t" <<ref_len <<" bp\t" <<t1-t0<< endl;


    // cout << "Index is done."<< "\t" << extract_ref_len;

    // fa_file.close();
    // index_file.close();
    // len_file.close();
    return kmer_hit_array;
}

class Select_ref{
    public:
    float get_fitness(string fasta_file, Encode encoder);
    string check_each_genome(string genome_list_file, Encode encoder, string match_rate_file);

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

    cout <<"check fitness..."<<endl;
    ref_len= ref_seq.length();
    cout <<"genome len is "<<ref_len<<endl;
 

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

char * generate_coder(bool coder_num)
{
    // A:65 97 T:116 84 C:99 67 G: 103 71
    static char coder [1000];
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

    return coder;
}

int * generate_base(int k)
{
    static int base [32];
    for (int i = 0; i<k; i++){       
        base[i] = pow(2, k-i-1);
    }
    return base;
}

char * generate_complement(void)
{
    static char comple [256];
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
    return comple;
}

short * random_coder(int k)
{
    static short permu[18] = {0,1,2,0,2,1,1,2,0,1,0,2,2,0,1,2,1,0};
    static short choose_coder[100];
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
    // for (int j = 0; j < 100; j++){
    //         cout << j <<"\t" << choose_coder[j] << endl;
    // }
    return choose_coder;
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

long * split_ref(string index_name, string fasta_file, int thread_num){
    ifstream len_file;
    len_file.open(fasta_file+".genome.len.txt", ios::in);

    static long split_ref_cutoffs[300];
    int cut_index = 0;
    long index_size = file_size(index_name);
    long each_index_size = index_size/thread_num + 1;

    int ref_len, ref_index;
    long slide_ref_len;
    string chr_name;

    long add;
    long pos = 100 * 4;
    long start_byte, end_byte;
    start_byte = pos;
    long start_ref_index = 1;
    long count_ref_index = 0;

    while(len_file >> chr_name >> ref_index >> ref_len >> slide_ref_len ){
        count_ref_index += 1;

        add = 4*((ref_len-k+1)*coder_num+1); //the size of the genome.
        if (pos -start_byte > each_index_size){
            end_byte = pos + add;
            cout << start_byte <<"--"<<end_byte<<"\t"<<ref_len<<"index"<< start_ref_index<<endl;
            split_ref_cutoffs[3*cut_index] = start_byte;
            split_ref_cutoffs[3*cut_index+1] = end_byte;
            split_ref_cutoffs[3*cut_index+2] = start_ref_index;
            cut_index += 1;
            start_byte = end_byte;     
            start_ref_index = count_ref_index + 1;   
        }
        pos += add;
        if (pos >= index_size){
            break;
        }
        
    }
    if (start_byte != index_size){
        end_byte = index_size;
        split_ref_cutoffs[3*cut_index] = start_byte;
        split_ref_cutoffs[3*cut_index+1] = end_byte;
        split_ref_cutoffs[3*cut_index+2] = start_ref_index;
        split_ref_cutoffs[3*cut_index+3] = 0; // indicate the end
        cout << start_byte <<"\t--\t"<<end_byte<<"\t"<<ref_len<<"index"<< start_ref_index<<endl;
    }
    len_file.close();
    return split_ref_cutoffs;
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
    int down_sam_ratio = stod(argv[8]); // percentage 30

    // string fq1 = argv[1];
    // string fq2 = argv[2];  
    // string fasta_file = argv[3];
    // string outdir = argv[4];
    // string ID = argv[5];
    // thread_num = 1; 
    // k = 26;
    // low_depth = 10;
    // array_size = pow(2, k);
    // kmer_count_table = new char[array_size];
    

    Encode encoder;
    encoder.constructer(k);

    


    
    // string fq1 = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/sim_s_lugdunensis//fastq/FDAARGOS_222.1.fq";
    // string fq2 = "/mnt/d/breakpoints/assembly/simulation/assembly_test/sim/sim_s_lugdunensis//fastq/FDAARGOS_222.2.fq";
    // string fasta_file = "/mnt/d/breakpoints/assembly/simulation/assembly_test/Staphylococcus_lugdunensis/GCF_008728755.1/GCF_008728755.1_ASM872875v1_genomic.fna";

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
    string select_genome = select.check_each_genome(fasta_list, encoder, match_rate_file); //"/mnt/d/breakpoints/assembly/sim/sim_s_lugdunensis/fasta/s_lugdunensis.database"

    return 0;
}

// Fasta fasta;
// string ref_seq = fasta.get_ref_seq(fasta_file);
// int ref_len = ref_seq.length();
// char* kmer_hit_array = read_ref(ref_seq, encoder, index_name);
// Unmap unmap;
// unmap.get_unmap(kmer_hit_array, ref_len);
// unmap.get_map_kmer(encoder, ref_seq);
// unmap.get_unmap_reads(fq1, encoder);
// time_t now3 = time(0);
// cout << "read ref done. "<< now3 - now2 << endl;
// unmap.get_unmap_reads(fq2, encoder);
// unmap.output_map_segments(fasta_file, ref_seq, outdir, ID);
// time_t now4 = time(0);
// cout << "Get unmap reads. "<< now4 - now3 << endl;