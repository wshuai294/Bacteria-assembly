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
const unsigned char least_depth = 255; // 127
const unsigned int window_size = 1000; // 1000
const unsigned int window_min_hit = 1000; // 999

double low_depth;
int k;
unsigned char *kmer_count_table;
unsigned char *map_kmer_table; // all kmers of mapped region
long array_size;
float least_read_kmer_ratio;
// vector<string> record_read_mapping_list;

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

double findMedian(unsigned char *arr) {
    int n = sizeof(arr) / sizeof(arr[0]);
    // Sort the array in ascending order
    sort(arr, arr + n);

    // If the array has odd length, return the middle element
    if (n % 2 == 1) {
        return arr[n / 2];
    }
    // If the array has even length, return the average of the two middle elements
    else {
        return (arr[n / 2 - 1] + arr[n / 2]) / 2.0;
    }
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
        cout <<  i << "\t" << pos <<endl; 
        if (done){
            break;
        }
    }
    return pos;
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

class Unmap{
    public:
        int save_map_region[10000]; // save segment intervals [start end start end ...]
        int save_map_index = 0;
        int final_map_region[10000];
        int final_map_index = 0;
        // char *map_kmer_table = new char[array_size];
        string ref_seq;
        void get_unmap(unsigned char* kmer_hit_array, int ref_len, string prefix);
        void get_map_kmer(Encode encoder, string ref_seq, std::vector<std::pair<int, int>> uniq_interval);
        void output_map_segments(string fasta_file, string ref_seq, string outdir, string ID, string ref_name);

};

void Unmap::get_unmap(unsigned char* kmer_hit_array, int ref_len, string prefix){
    cout << "ref len is: \t"<<ref_len<<endl; 
    double median_count = findMedian(kmer_hit_array);
    double depth_cutoff = 1.9 * median_count;
    low_depth = 0.3 * median_count;
    cout <<"\n<<<<<Median kmer count in ref is:\t" << median_count <<endl; 
    ofstream out_count_file;
    out_count_file.open(prefix + ".ref.kmer.count.tab", ios::out);
    
    int start = 0;
    int end = 0;
    int good_base_num = 0; // in a windows
    int window_depth_sum = 0; 
    bool good_flag = false;
    bool depth_flag = false;
    int each_base = 0;
    int small_window = 200;
    int overlap_cutoff = 0;
    // short record_each_base[window_size];
    queue<int>record_each_base;
    queue<int>record_each_depth;
    for (int j = 0; j < ref_len; j++){
        if (kmer_hit_array[j] >= low_depth){
            each_base = 1;
        }
        else{
            each_base = 0;
        }
        out_count_file << j << "\t" << int(kmer_hit_array[j]) <<endl; 

        if (j < window_size){
            record_each_base.push(each_base);
            good_base_num += each_base;
        }
        else{
            good_base_num = good_base_num - record_each_base.front() + each_base;
            record_each_base.pop();
            record_each_base.push(each_base);
        }

        if (j < small_window){
            record_each_depth.push(kmer_hit_array[j]);
            window_depth_sum += kmer_hit_array[j];
        }
        else{
            window_depth_sum = window_depth_sum - record_each_depth.front() + kmer_hit_array[j];
            record_each_depth.pop();
            record_each_depth.push(kmer_hit_array[j]);
        }

        double window_depth = window_depth_sum / small_window;

        // if (window_depth > 100 ){
        //     cout << j << "\t" << window_depth <<endl;
        // }

        depth_flag = true;
        if (window_depth < depth_cutoff){
            depth_flag = true;
        }
        else{
            depth_flag = false;
        }
        
        // cout << window_depth <<"\t" <<depth_flag<<"\t"<<depth_cutoff<<endl; 
            
        end += 1;
        // cout << j <<"\t" <<good_base_num<<"\t"<<ref_len<<endl; 
        if (good_base_num >= window_min_hit && depth_flag){
            if (good_flag == false){
                start = end + 1 - window_size;
                if (start < 0){
                    start = 1;
                }  
            } 
            good_flag = true;
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
                    if (start - save_map_region[save_map_index*2-1] < overlap_cutoff){
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
        if (start - save_map_region[save_map_index*2-1] < overlap_cutoff){
            save_map_region[save_map_index*2-1] = end;
        }
        else{
            save_map_region[save_map_index*2] = start;
            save_map_region[save_map_index*2+1] = end; 
            save_map_index += 1;  
        }
    }
    cout <<"\n mapped segment num is \t" <<save_map_index << "\twindow_min_hit:\t"<<window_min_hit <<endl;
    delete []kmer_hit_array;
    out_count_file.close();
}

void Unmap::get_map_kmer(Encode encoder, string ref_seq, std::vector<std::pair<int, int>> uniq_interval){
    unsigned int kmer_index, comple_kmer_index, real_index;
    int m;
    int map_len = 0;
    // check whether this region locate in the unique seq of ref
    for (int z = 0; z < save_map_index; z++){
        int start = save_map_region[z*2];
        int end = save_map_region[z*2+1];
        int clean_start = 0;
        int clean_end = 0;
        map_len += (end - start);
        cout <<z <<" map seg\t" <<start <<"\t" <<end<<endl;
        for (const auto& interval : uniq_interval) {
            bool overlap = false;

            if (start >= interval.first &&  start <= interval.second){
                if (end <= interval.second){
                    clean_start = start;
                    clean_end = end;
                }
                else{
                    clean_start = start;
                    clean_end = interval.second;                   
                }
                overlap = true;
            }
            else{
                if (end >= interval.first &&  end <= interval.second){
                    clean_start = interval.first;
                    clean_end = end;
                    overlap = true;
                } 
                else{
                    if (interval.first >= start && interval.first <= end ){
                        overlap = true;
                        clean_start = interval.first;
                        if (interval.second <= end){
                            clean_end = interval.second;
                        }
                        else{
                            clean_end = end;
                        }
                    }
                }
            }

            if (overlap == true){
                cout <<final_map_index <<"\t" << clean_start <<"\t" << clean_end <<"\t" << interval.first << "\t" << interval.second<< "\t" << start << "\t" << end <<endl;
                if (clean_end - clean_start > 200){
                    final_map_region[2*final_map_index] = clean_start;
                    final_map_region[2*final_map_index+1] = clean_end;
                    final_map_index += 1;
                }
            } 
        }
    }
    cout <<" map len\t" << map_len << endl;
    cout << "<<<<<<<<<<<<<<<< final segment number "<< final_map_index<<endl;
    for (int z = 0; z < final_map_index; z++){
        int start = final_map_region[z*2];
        int end = final_map_region[z*2+1];      

        string segment = ref_seq.substr(start, end);
        int seg_len = segment.length();

        cout <<z <<" map and uniq seg\t" <<start <<"\t" <<end<<endl;

        int *ref_int = new int[seg_len];
        int *ref_comple_int = new int[seg_len];
        for (int j = 0; j < seg_len; j++){
            ref_int[j] = (int)segment[j];
            ref_comple_int[j] = encoder.comple[ref_int[j]];
        }
        for (int j = 0; j < seg_len-encoder.k+1; j++){
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
                map_kmer_table[real_index] = 1;
            } 
        }
        delete [] ref_int;
        delete [] ref_comple_int;
    }
}

void Unmap::output_map_segments(string fasta_file, string ref_seq, string outdir, string ID, string ref_name){
    int start = 0;

    int end = 0;
    ofstream out_fa_file;
    out_fa_file.open(outdir+ "/"+ ID + ".map.fasta", ios::out);
    string map_segment_seq, map_segment_name;
    for (int j = 0; j < final_map_index; j++){
        start = final_map_region[j*2];
        end = final_map_region[j*2+1];
        // cout <<j <<"final\t" <<start <<"\t" <<end<<endl;
        map_segment_seq = ref_seq.substr(start, end-start);
        map_segment_name = ">" + ref_name + ":" + to_string(start) + "-" + to_string(end);

        out_fa_file << map_segment_name<<endl; 
        out_fa_file << map_segment_seq<<endl;
    }
    out_fa_file.close();
}

std::vector<std::pair<int, int>> get_ref_repeat(string ref_seq, Encode encoder, string prefix)
{
    unsigned char *kmer_count_table;
    kmer_count_table = new unsigned char[array_size];
    memset(kmer_count_table, 0, sizeof(kmer_count_table));


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
    short convert_ref[300];
    short complemented_ref[300];

    ofstream out_count_file;
    out_count_file.open(prefix + ".ref.origin.kmer.count.tab", ios::out);

    ref_len= ref_seq.length();
    unsigned char *ref_kmer_count = new unsigned char[ref_len];

    int *ref_int = new int[ref_len];
    int *ref_comple_int = new int[ref_len];
    for (int j = 0; j < ref_len; j++){
        ref_int[j] = (int)ref_seq[j];
        ref_comple_int[j] = encoder.comple[ref_int[j]];
    }

    for (int j = 0; j < ref_len-encoder.k+1; j++){
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
                kmer_count_table[real_index] += 1;
            }  
        }
    }    


    for (int j = 0; j < ref_len-encoder.k+1; j++){
        int max_count = 0;
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
                if (kmer_count_table[real_index] > max_count){
                    max_count = kmer_count_table[real_index];
                }
            }  
        }
        out_count_file << j << "\t" << max_count << endl;
        ref_kmer_count[j] = max_count;
    }   



    delete [] ref_int;
    delete [] ref_comple_int;
    delete [] kmer_count_table;
    out_count_file.close();

    // return ref_kmer_count;


    int window = 200;
    int repeat_num = 0;
    std::vector<int> count_list;

    std::cout << "# get unique \n"<< std::endl;
    int length = sizeof(ref_kmer_count) / sizeof(ref_kmer_count[0]);
    std::cout << "ref length " << length << std::endl;
    std::cout << "ref length " << sizeof(ref_kmer_count) << std::endl;

    for (int j = 0; j < ref_len; ++j ){
        int count = int(ref_kmer_count[j]);
        if (count > 1) {
            count_list.push_back(1);
        }
        else {
            count_list.push_back(0);
        }
    }

    int start = 0, end = 0;
    bool continuous = false;

    std::vector<std::pair<int, int>> uniq_interval;
    for (int i = 0; i < ref_len - window + 1; ++i) {
        repeat_num = 0;
        for (int j = i; j < i + window; ++j) {
            repeat_num += count_list[j];
        }
        double repeat_ratio = repeat_num / static_cast<double>(window);
        // std::cout << i << "\t" << repeat_ratio << std::endl;

        if (repeat_ratio < 0.9) {
            if (!continuous) {
                start = i;
                end = i + window;
            }
            else {
                end = i + window;
            }
            continuous = true;
        }
        else {
            if (continuous) {
                if (!uniq_interval.empty() && start <= uniq_interval.back().second) {
                    uniq_interval.back().second = end;
                }
                else {
                    uniq_interval.push_back(std::make_pair(start, end));
                }
            }
            continuous = false;
        }
    }

    if (continuous) {
        if (!uniq_interval.empty() && start <= uniq_interval.back().second) {
            uniq_interval.back().second = end;
        }
        else {
            uniq_interval.push_back(std::make_pair(start, end));
        }
    }

    // int uniq_len = 0;
    // for (const auto& interval : uniq_interval) {
    //     int interval_len = interval.second - interval.first;
    //     if (interval_len < 10000) {
    //         continue;
    //     }
    //     std::cout << "# "<< interval.first << ", " << interval.second << std::endl;
    //     uniq_len += interval_len;
    // }

    // std::cout << "Unique ref length is: " << uniq_len << std::endl;
    // Printing count_list[2510352:2512217] is not provided in the original code

    // return 0;
    return uniq_interval;

}

unsigned char* read_ref(string ref_seq, Encode encoder)
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
    short convert_ref[300];
    short complemented_ref[300];

    ref_len= ref_seq.length();
    unsigned char *kmer_hit_array = new unsigned char[ref_len];

    int *ref_int = new int[ref_len];
    int *ref_comple_int = new int[ref_len];
    for (int j = 0; j < ref_len; j++){
        ref_int[j] = (int)ref_seq[j];
        ref_comple_int[j] = encoder.comple[ref_int[j]];
    }
    // index_file.write((char *)(&ref_len), sizeof(unsigned int));
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
                dp = kmer_count_table[real_index];
            }
            
            // if (max_dp < dp){  // get maximum kmer count across hash functions
            //     max_dp = dp;
            // }

            if (i == 0){  // get minimum kmer count across hash functions
                max_dp = dp;
            }
            else{
                if (dp < max_dp){
                    max_dp = dp;
                }
            }
            // cout <<j<<"\t"<< i<<"\t"<<dp<<"\t" <<endl;
            // index_file.write((char *)(&real_index), sizeof(real_index));  
        }
        // kmer_hit_array[j] = min_dp; //max_dp; 
        kmer_hit_array[j] = max_dp;
    }        
    delete [] ref_int;
    delete [] ref_comple_int;

    time_t t1 = time(0);
    return kmer_hit_array;
}

class Select_ref{
    public:
    float get_fitness(string fasta_file, Encode encoder);
    string check_each_genome(string genome_list_file, Encode encoder, string match_rate_file);
};

float Select_ref::get_fitness(string fasta_file, Encode encoder){
    Fasta fasta;
    fasta.get_ref_seq(fasta_file);
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
    ref_len= fasta.ref_seq.length();
    // cout <<"genome len is "<<ref_len<<endl;
 

    int *ref_int = new int[ref_len];
    int *ref_comple_int = new int[ref_len];
    for (int j = 0; j < ref_len; j++){
        ref_int[j] = (int)fasta.ref_seq[j];
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
        if (match_rate > 0.99){ // stop searching if the match rate is over 0.99
            break;
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

    // long pos = 0;
    // for (long i = start; i>0; i--){
    //     fq_file.seekg(i, ios::beg);
    //     char j;
    //     fq_file.get(j);
    //     if (j == '@'){ //only read name has this symbol.
    //         pos = i;
    //         break;
    //     }       
    // }
    long pos = get_fq_start(fq_file, start);
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

void get_unmap_fastq(string fastq_file, Encode encoder, int down_sam_ratio, long start, long end, string prefix, string fq_index)
{
    time_t t0 = time(0);
    ifstream fq_file; 
    fq_file.open(fastq_file);

    ofstream unmap_fq_file;
    string unmap_file = prefix + "_" + fq_index + "_part_" + to_string(start)+".txt";
    unmap_fq_file.open(unmap_file, ios::out);
    cout << unmap_file << endl;

    // long pos = 0;
    // for (long i = start; i>0; i--){
    //     fq_file.seekg(i, ios::beg);
    //     char chr1, chr2, chr3;
    //     fq_file.get(chr1);
    //     fq_file.get(chr2);
    //     fq_file.get(chr3);

    //     if (chr1 != '+' and chr2 == '\n' and chr3 == '@'){ //only read name has this symbol at the front.
    //         pos = i + 2;
    //         break;
    //     }       
    // }
    long pos = get_fq_start(fq_file, start);

    fq_file.seekg(pos, ios::beg);
    long add_size = start;

    string reads_seq;
    int reads_int [300];
    int reads_comple_int [300];

    unsigned int lines = 0;
    int converted_reads [900];
    int complemented_reads [900];
    int m,n,r;
    unsigned int kmer_index, comple_kmer_index, real_index, b;   
    short read_len = 0;
    bool all_valid;
    int read_index = 0;
    bool unmap_flag;
    string read_name;

    while (getline(fq_file,reads_seq))
    {
        if (add_size>=end){
            break;
        }
        add_size += reads_seq.length();

        if (lines % 4 == 1){
            read_len = reads_seq.length();//cal read length
            r = rand() % 100 ;
            if (r < down_sam_ratio){
                for (int j = 0; j < read_len; j++){
                    reads_int[j] = (int)reads_seq[j];
                    reads_comple_int[j] = encoder.comple[reads_int[j]];
                }
                int hit_kmer_num = 0;
                int N_num = 0;
                for (int j = 0; j < read_len-encoder.k+1; j++){
                    bool hit_flag = false;
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        all_valid = true;
                        for (int z = 0; z<encoder.k; z++){
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
                        if (all_valid == false){
                            real_index = 0;
                        }
                        if ((int)map_kmer_table[real_index] > 0 and all_valid == true ){
                            hit_flag = true;
                        }  
                    }
                    if (hit_flag){
                        hit_kmer_num += 1;
                    }
                    if (! all_valid){
                        N_num += 1; // count num of N
                    }
                }
                int cutoff = least_read_kmer_ratio * (read_len - k + 1);
                int max_N = 0.1 * (read_len - k + 1);
                // cout << max_N << " "<< N_num << endl;  
                if ( hit_kmer_num < cutoff && N_num < max_N){ // determine if a read is unmap , should not contain too many Ns
                    unmap_flag = true;
                }
                else {
                    unmap_flag = false;
                }
                // cout << lines << " " << hit_kmer_num << " " << unmap_flag << " " << cutoff << endl;  
            }         
        }

        if (lines % 4 == 0){
            unmap_flag = false;
            read_index += 1;
            read_name = get_read_ID(reads_seq).substr(1);
        }
      
        if (lines % 4 == 3){
            // cout << unmap_flag<< "\tread_name\t" << read_name << "\t" << endl;
            if (unmap_flag){
                unmap_fq_file << read_name << endl;
            }
        }
        lines++;
    }
    fq_file.close();
    unmap_fq_file.close();
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

long file_size(string filename){  
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
    string outdir = argv[7]; // 
    array_size = pow(2, k);
    kmer_count_table = new unsigned char[array_size];
    int down_sam_ratio = stod(argv[8]); // percentage (0-100), randomly select reads with this probability, 30
    string ID = argv[9]; // sample ID
    least_read_kmer_ratio = stod(argv[10]); // kmer num less than this will be regarded as unmapped
    string match_rate_file = outdir + "/" + ID + ".match_rate.csv"; //output file, the match rate of each fasta
    string prefix = outdir + "/" + ID ;

    
    Encode encoder;
    encoder.constructer(k);
    
    std::vector<std::thread>threads;
    long start = 0;
    long end = 0;
    long size = file_size(fq1);
    long each_size = size/thread_num;  
    cout <<"reads size\t"<<size<<endl;
    
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
    time_t now3 = time(0);
    cout << "ref is selected. "<< now3 - now2<<endl;

    memset(kmer_count_table, 0, sizeof(kmer_count_table));
    down_sam_ratio = 100;
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

    time_t now4 = time(0);
    cout << "reads are reloaded. "<< now4 - now3<<endl;

    Fasta fasta;
    fasta.get_ref_seq(select_genome);
    // string ref_name = fasta.get_ref_name(select_genome);
    int ref_len = fasta.ref_seq.length();
    cout << "ref\t" << fasta.ref_name << "\tlen:\t" << ref_len << "\t" << select_genome <<endl;
    unsigned char* kmer_hit_array = read_ref(fasta.ref_seq, encoder);
    delete [] kmer_count_table;
    std::vector<std::pair<int, int>> uniq_interval = get_ref_repeat(fasta.ref_seq, encoder, prefix);
    // get_uniq_intervals(ref_kmer_count);
    int uniq_len = 0;
    for (const auto& interval : uniq_interval) {
        int interval_len = interval.second - interval.first;
        std::cout << "# "<< interval.first << ", " << interval.second << std::endl;
        uniq_len += interval_len;
    }

    std::cout << "Unique ref length is: " << uniq_len << std::endl;
    
    
    Unmap unmap;
    unmap.get_unmap(kmer_hit_array, ref_len, prefix);
    map_kmer_table = new unsigned char[array_size];
    unmap.get_map_kmer(encoder, fasta.ref_seq, uniq_interval);
    unmap.output_map_segments(select_genome, fasta.ref_seq, outdir, ID, fasta.ref_name);
    time_t now5 = time(0);
    cout << "Getting map interval is done. "<< now5 - now4 << endl;

    down_sam_ratio = 100;
    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(get_unmap_fastq, fq1, encoder, down_sam_ratio, start, end, prefix, "fq1"));
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
        threads.push_back(thread(get_unmap_fastq, fq2, encoder, down_sam_ratio, start, end, prefix, "fq2"));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();

    time_t now6 = time(0);
    cout << "Get unmap reads. "<< now6 - now5 << endl;
    cout << "Whole process takes  "<< now6 - now1 << endl;

    delete [] map_kmer_table;

    return 0;
}


