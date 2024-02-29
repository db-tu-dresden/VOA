#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>

#include <immintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

#include <chrono>

#include <algorithm>

#include <omp.h> // unused rn

#define VECTOR_ELEMENT_COUNT 8

using time_stamp = std::chrono::high_resolution_clock::time_point;
time_stamp time_now(){
    return std::chrono::high_resolution_clock::now();
}

uint64_t duration_time (time_stamp begin, time_stamp end){
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
}

uint64_t duration_time_milliseconds (std::chrono::high_resolution_clock::time_point begin, std::chrono::high_resolution_clock::time_point end){
    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
}

uint64_t duration_time_seconds (time_stamp begin, time_stamp end){
    return std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
}


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


size_t noise(size_t position, size_t seed){
    size_t BIT_NOISE1 = 0x68E31DA4;
    size_t BIT_NOISE2 = 0xB5297A4D;
    size_t BIT_NOISE3 = 0x1B56C4E9;

    uint64_t mangled = position;
    mangled *= BIT_NOISE1;
    mangled += seed;
    mangled *= BIT_NOISE2;
    mangled ^= (mangled << 13);
    mangled += BIT_NOISE1;
    mangled ^= (mangled >> 7);
    mangled *= BIT_NOISE3;
    mangled ^= (mangled << 17);
    return mangled;
}

void generate_permutation(uint64_t* permutations, size_t permutation_size, size_t permutation_count, size_t total_elements, size_t seed){
    size_t run_nr = 0;

    for(size_t p_nr = 0; p_nr < permutation_count; p_nr++){
start_generation:
        uint64_t * help = &permutations[p_nr * permutation_size];
        help[0] = (p_nr % total_elements) + 1;
        for(size_t p_element = 1; p_element < permutation_size; p_element++){
            size_t element;
            bool fine;
            do{
                fine = true;
                element = (noise(run_nr++ + help[p_element - 1],  seed + p_nr) % total_elements) + 1;
                for(size_t check = 0; check < p_element; check++){
                    if(element == help[check]){
                        fine = false;
                        break;
                    }
                }
            }while(!fine);
            permutations[p_nr * permutation_size + p_element] = element;
        }
        for(size_t check_1 = 0; check_1 < p_nr; check_1++){
            bool same = true;
            for(size_t check_2 = 0; check_2 < permutation_size; check_2 ++){
                if(permutations[check_1 * permutation_size + check_2] != help[check_2]){
                    same = false;
                    break;
                }
            }
            if(same){ 
                goto start_generation;
            }
        }
    }
}

void generate_build_data(uint64_t* values, size_t size){
    for(size_t i = 0; i < size; ++i){
        values[i] = i + 1;
    }
}

size_t generate_probe_data(uint64_t*& values, size_t size, size_t key_amount, size_t seed){
    size_t permutation_count = key_amount + 1;
    size_t permutation_length = VECTOR_ELEMENT_COUNT; 
    uint64_t* permutations = (uint64_t*) malloc(permutation_count * permutation_length * sizeof(uint64_t));
    generate_permutation(permutations, permutation_length, permutation_count, key_amount, seed);
    
    size_t i = 0, pos = 0;
    for(i = 0, pos = 0; i < (size / permutation_length) &&  pos < size; i++){
        size_t p_sel = noise(i, seed + 1) % permutation_count;
        uint64_t * help = &permutations[p_sel * permutation_length];
        for(size_t e = 0; e < permutation_length; e++, pos++){
            values[pos] = help[e];
        }
    }

    free(permutations);
    return pos; 
}

void generate_table(uint64_t *table, size_t t_size, uint64_t* build, size_t b_size, size_t chunk_size_elements){
    size_t shift_val = __builtin_ctz(chunk_size_elements);
    for(size_t i = 0; i < t_size; i++){
        table[i] = 0;
    }

    for(size_t i = 0; i < b_size; i++){
        size_t pos = (build[i] - 1) << shift_val;
        table[pos] = build[i];
    }

}

size_t probe_scalar(uint64_t *table, size_t t_size, uint64_t *lookup, size_t l_size, size_t chunk_size_elements){
    size_t shift_val = __builtin_ctz(chunk_size_elements);
    size_t hit_count = 0;
    for(size_t i = 0; i < l_size; i++){
        uint64_t current_val = lookup[i];
        size_t pos = (current_val - 1) << shift_val;
        bool res_val = table[pos] == current_val;
        hit_count += res_val;
    }
    return hit_count;
}

size_t probe_simd_horizontal(uint64_t *table, size_t t_size, uint64_t *lookup, size_t l_size, size_t chunk_size_elements){
    size_t shift_val = __builtin_ctz(chunk_size_elements);
    size_t hit_count = 0;
    for(size_t i = 0; i < l_size; i++){
        uint64_t current_val = lookup[i];
        __m512i current_vec = _mm512_set1_epi64(current_val);
        size_t pos = (current_val - 1) << shift_val;
        __m512i table_val = _mm512_loadu_epi64(&table[pos]);
        __mmask8 res_val = _mm512_cmp_epi64_mask(current_vec, table_val, _MM_CMPINT_EQ);
        
        hit_count += __builtin_popcount(res_val);
    }
    return hit_count;
}

size_t probe_simd_voa(uint64_t *table, size_t t_size, uint64_t * lookup, size_t l_size, size_t chunk_size_elements){
    size_t shift_val = __builtin_ctz(chunk_size_elements);
    __m512i zero_512i = _mm512_setzero_si512();
    size_t hit_count = 0;
    for(size_t i = 0; i < l_size; i+= VECTOR_ELEMENT_COUNT){
        __m512i current_vec = _mm512_loadu_epi64(&lookup[i]);
        __m512i position_vec = _mm512_slli_epi64(_mm512_sub_epi64(current_vec, _mm512_set1_epi64(1)), shift_val);
        __m512i table_val = _mm512_i64gather_epi64( position_vec, table, 8);
        __mmask8 res_val = _mm512_cmp_epi64_mask(current_vec, table_val, _MM_CMPINT_EQ);
        
        hit_count += __builtin_popcount(res_val);
    }
    return hit_count;
}

// TODO
// size_t probe_simd_horizontal_multi_threaded(uint64_t *table, size_t t_size, uint64_t *lookup, size_t l_size, size_t chunk_size_elements){
//     size_t shift_val = __builtin_ctz(chunk_size_elements);
//     size_t hit_count = 0;
//     for(size_t i = 0; i < l_size; i++){
//         uint64_t current_val = lookup[i];
//         __m512i current_vec = _mm512_set1_epi64(current_val);
//         size_t pos = (current_val - 1) << shift_val;
//         __m512i table_val = _mm512_loadu_epi64(&table[pos]);
//         __mmask8 res_val = _mm512_cmp_epi64_mask(current_vec, table_val, _MM_CMPINT_EQ);
        
//         hit_count += __builtin_popcount(res_val);
//     }
//     return hit_count;
// }

// size_t probe_simd_voa_multi_threaded(uint64_t *table, size_t t_size, uint64_t * lookup, size_t l_size, size_t chunk_size_elements){
//     size_t shift_val = __builtin_ctz(chunk_size_elements);
//     __m512i zero_512i = _mm512_setzero_si512();
//     size_t hit_count = 0;
//     for(size_t i = 0; i < l_size; i+= VECTOR_ELEMENT_COUNT){
//         __m512i current_vec = _mm512_loadu_epi64(&lookup[i]);
//         __m512i position_vec = _mm512_slli_epi64(_mm512_sub_epi64(current_vec, _mm512_set1_epi64(1)), shift_val);
//         __m512i table_val = _mm512_i64gather_epi64( position_vec, table, 8);
//         __mmask8 res_val = _mm512_cmp_epi64_mask(current_vec, table_val, _MM_CMPINT_EQ);
        
//         hit_count += __builtin_popcount(res_val);
//     }
//     return hit_count;
// }

void sort(uint64_t* vals, size_t size){
    for(size_t i = 0; i < size; i ++){
        for(size_t e = i + 1; e < size; e++){
            if(vals[i] > vals[e]){
                uint64_t help = vals[i];
                vals[i] = vals[e];
                vals[e] = help;
            } 
        }
    }
}

int main(int argc, char** argv){
    if(argc >= 2 && cmdOptionExists(argv, argv+argc, "-h")){
        std::cout << "This is a simple baseline project for VOA.\n" 
            << "Currently, we measure  random access times of uload and gather in the context of Hash Tables\n"
            << "the options you have to run this project are:\n"
            << "\t-k <amount>\tdefault: 32\t\tthe amount is how many different integers are in the table\n\t\t\t\t\t\t\t(influences the table size)\n"
            << "\t-kb <amount>\t\t\t\tsimilar to -k but sets the table size in MiByte\n"
            << "\t-s <amount>\tdefault: 9 max: 9\tthe shift amount to generate a stride\n"
            << "\t-p <amount>\tdefault: 1\t\tthe amount of data that should be probed for in GiByte\n"
            << "\t-r <amount>\tdefault: 7\t\tthe number of repeats each measurement should be repeated\n"
            << "\t-t <amount>\tdefault: 1\t\tthe number of threats used to run the benchmark \\unused rn\n"
            << "\t-c <amount>\tdefault: 0\t\tthe number of equal elements in a permutation\n"
            << "\t-h\t\t\t\t\tprints this help\n";
        exit(0);
    }   
    char * input_k = getCmdOption(argv, argv+argc,"-k");
    char * input_kb = getCmdOption(argv, argv+argc,"-kb");
    char * input_s = getCmdOption(argv, argv+argc,"-s");
    char * input_p = getCmdOption(argv, argv+argc,"-p");
    char * input_r = getCmdOption(argv, argv+argc,"-r");
    char * input_t = getCmdOption(argv, argv+argc,"-t");
    char * input_t = getCmdOption(argv, argv+argc,"-c");

    size_t key_amount = 32;
    size_t shift_size = 9;
    size_t probe_amount = 1024 * 1024 * 128; // 1 GiByte
    size_t repeats = 7;

    if(input_p){
        double a = atof(input_p);
        probe_amount *= a;
    }

    if(input_r){
        size_t a = atoi(input_r);
        repeats = a;
        if(repeats <= 0){
            repeats = 1;
        }
    }

    if(input_s){
        shift_size = atoi(input_s);
        if(shift_size < 0){
            shift_size = 0;
        }else if(shift_size > 9){
            shift_size = 9;
        }
    }

    if(input_k){
        key_amount =  atol(input_k);
    }else if(input_kb){
        double nkey = atof(input_kb);
        size_t nnkey = nkey * 1024 * 1024;
        nnkey >>= shift_size;
        nnkey /= sizeof(uint64_t);
        if(nnkey < VECTOR_ELEMENT_COUNT){
            nnkey = VECTOR_ELEMENT_COUNT;
        }
        key_amount = nnkey;
    }

    size_t create_amount = key_amount;
    size_t chunk_size = 1 << shift_size;

    size_t table_size = (key_amount << shift_size); 
    size_t table_size_alloc = table_size + VECTOR_ELEMENT_COUNT; // we allocate more memory than we want to use to avoid segmentation faults

    uint64_t* table = (uint64_t*)malloc(table_size_alloc * sizeof(uint64_t));
    uint64_t* table_data = (uint64_t*)malloc(create_amount * sizeof(uint64_t));
    uint64_t* probe_data = (uint64_t*)malloc(probe_amount * sizeof(uint64_t));
    
    std::cout << "---------------------------------------------------------------------------\n";
//---Generating-Hash-Table-and-Printing-Info---------
    size_t table_size_byte = table_size * sizeof(uint64_t);
    double table_size_kibyte = table_size_byte / 1024;
    double table_size_mibyte = table_size_kibyte / 1024;
    double table_size_gibyte = table_size_mibyte / 1024;


    std::cout << "Hash Table Information:\t\t"; 
    if(table_size_byte < 1024){
        std::cout << table_size_byte << " Byte";
    }else if(table_size_byte >= 1024 && table_size_kibyte < 1024){
        std::cout << table_size_kibyte << " kiB";
    }else if(table_size_kibyte >= 1024 && table_size_mibyte < 1024){
        std::cout << table_size_mibyte << " MiB"; 
    }else{
        std::cout << table_size_gibyte << " GiB";
    }
    std::cout << "\t\t" << table_size << " Buckets\n";
    std::cout << "\t" << key_amount << " different keys\t" << chunk_size << " stride size\n";

    std::cout << "---------------------------------------------------------------------------\n";

    generate_build_data(table_data, create_amount);
    generate_table(table, table_size, table_data, create_amount, chunk_size);
//---Generating-Probe-Data-and-Printing-Info---------
    probe_amount = generate_probe_data(probe_data, probe_amount, key_amount, 0xadd230b);

    size_t probe_amount_byte = probe_amount * sizeof(uint64_t);
    double probe_amount_kibyte = probe_amount_byte / 1024;
    double probe_amount_mibyte = probe_amount_kibyte / 1024;
    double probe_amount_gibyte = probe_amount_mibyte / 1024;

    std::cout << "Probe Data Information:\t\t"; 
    if(probe_amount_byte < 1024){
        std::cout << probe_amount_byte << " Byte";
    }else if(probe_amount_byte >= 1024 && probe_amount_kibyte < 1024){
        std::cout << probe_amount_kibyte << " kiB";
    }else if(probe_amount_kibyte >= 1024 && probe_amount_mibyte < 1024){
        std::cout << probe_amount_mibyte << " MiB"; 
    }else{
        std::cout << probe_amount_gibyte << " GiB";
    }
    std::cout << "\t" << probe_amount << " Values" << std::endl;
    std::cout << "---------------------------------------------------------------------------\n";


//---Sanity-Checking-that-both-algs-have-same-result-
    size_t a,b,c;
    std::cout << "sanity check: " << std::flush;
    a = probe_scalar(table, table_size, probe_data, probe_amount, chunk_size);
    std::cout << a << "\t" << std::flush;
    b = probe_simd_horizontal(table, table_size, probe_data, probe_amount, chunk_size);
    std::cout << b << "\t" << std::flush;
    c = probe_simd_voa(table, table_size, probe_data, probe_amount, chunk_size); 
    std::cout << c << "\t" << std::flush; 
    std::cout << a - b << " " << a - c << " " << b - c << std::endl;    

//---Timeing-----------------------------------------
    uint64_t median, mean, min, max;
    size_t ignore_best_and_worst_x = (repeats/6);
    size_t time_ms[repeats];
    size_t x = 0; // used for evading compiler optimization
    std::cout << "---------------------------------------------------------------------------\n";
    std::cout << "  " << repeats <<" repeats ± " << ignore_best_and_worst_x <<"\tmedian\tmean\tmax\tmin\n";
//---Scalar-Testing------------------------------
    std::cout << "Scalar" << std::flush;
    for(size_t i = 0; i < repeats; i++){
        time_stamp b = time_now();
        x += probe_scalar(table, table_size, probe_data, probe_amount, chunk_size);
        time_stamp e = time_now();
        time_ms[i] = duration_time_milliseconds(b, e);
    }

    sort(time_ms, repeats);
    median = time_ms[repeats/2];
    min = time_ms[ignore_best_and_worst_x];
    max = time_ms[repeats - ignore_best_and_worst_x - 1];
    mean = 0;
    for(size_t i = ignore_best_and_worst_x; i < repeats -ignore_best_and_worst_x; i ++){
        mean += time_ms[i];
    }
    mean /= (repeats - ignore_best_and_worst_x * 2);
    std::cout <<" ms\t\t" << median << "\t" << mean << "\t" << max << "\t" << min << std::endl;
    std::cout << "   Million Probes/s\t" << (probe_amount * 1000.) / (median * 1000 * 1000) << "\t\t\t\t" << (x & 0xF) << std::endl;
//---Horizontal-Testing------------------------------
    std::cout << "Horizontal" << std::flush;
    for(size_t i = 0; i < repeats; i++){
        time_stamp b = time_now();
        x += probe_simd_horizontal(table, table_size, probe_data, probe_amount, chunk_size);
        time_stamp e = time_now();
        time_ms[i] = duration_time_milliseconds(b, e);
    }

    sort(time_ms, repeats);
    median = time_ms[repeats/2];
    min = time_ms[ignore_best_and_worst_x];
    max = time_ms[repeats - ignore_best_and_worst_x - 1];
    mean = 0;
    for(size_t i = ignore_best_and_worst_x; i < repeats -ignore_best_and_worst_x; i ++){
        mean += time_ms[i];
    }
    mean /= (repeats - ignore_best_and_worst_x * 2);
    std::cout <<" ms\t\t" << median << "\t" << mean << "\t" << max << "\t" << min << std::endl;
    std::cout << "   Million Probes/s\t" << (probe_amount * 1000.) / (median * 1000 * 1000) << "\t\t\t\t" << (x & 0xF) << std::endl;

//---VOA-Testing-------------------------------------
    std::cout << "VoA" << std::flush;
    for(size_t i = 0; i < repeats; i++){
        time_stamp b = time_now();
        x += probe_simd_voa(table, table_size, probe_data, probe_amount, chunk_size);
        time_stamp e = time_now();
        time_ms[i] = duration_time_milliseconds(b, e);
    }

    sort(time_ms, repeats);
    median = time_ms[repeats/2];
    min = time_ms[ignore_best_and_worst_x];
    max = time_ms[repeats - ignore_best_and_worst_x - 1];
    mean = 0;
    for(size_t i = ignore_best_and_worst_x; i < repeats -ignore_best_and_worst_x; i ++){
        mean += time_ms[i];
    }
    mean /= (repeats - ignore_best_and_worst_x * 2);
    std::cout <<" ms\t\t\t" << median << "\t" << mean << "\t" << max << "\t" << min << std::endl;
    std::cout << "   Million Probes/s\t" << (probe_amount * 1000.) / (median * 1000 * 1000) << "\t\t\t\t" << (x & 0xF) << std::endl;
    std::cout << "---------------------------------------------------------------------------\n\n";
}