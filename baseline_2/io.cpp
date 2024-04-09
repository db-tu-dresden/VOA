#include <filesystem>

#include "io.hpp"

bool file_exists (std::string name){
    std::ifstream f(name.c_str());
    return f.good();
}

void write_to_file(
    std::string filename,
    std::string content,
    bool override
){
    std::ofstream myfile;
    if(override)
        myfile.open(filename, std::ios::out);
    else
        myfile.open(filename, std::ios::app);
    if(myfile.is_open()){
        myfile << content << "\n";
        myfile.close();
    }else{
        throw std::runtime_error("Unable to open the file to write the results!\n");
    }
}

void write_to_file(
    std::string filename,
    std::string content,
    uint64_t* time_in_ms,
    size_t count,
    bool override
){
    std::ofstream myfile;
    if(override)
        myfile.open(filename, std::ios::out);
    else
        myfile.open(filename, std::ios::app);
    if(myfile.is_open()){
        for(size_t i = 0; i < count; i++){
            myfile << content << time_in_ms[i] << "\n";
        }
        myfile.close();
    }else{
        throw std::runtime_error("Unable to open the file to write the results!\n");
    }
}

void create_raw_file(std::string filename){
    std::stringstream c_filename_ss;
    c_filename_ss << filename << "_raw.csv";
    std::string fname = c_filename_ss.str();

    if(!file_exists(fname)){
        std::stringstream header;
        header << "table_size_byte,probe_size_byte,key_count,collisions,thread_count,step_size,algorithm,time_ms";
        write_to_file(fname, header.str(), true);
    } 
}

void create_summary_file(std::string filename){
    std::stringstream c_filename_ss;
    c_filename_ss << filename << "_summary.csv";
    std::string fname = c_filename_ss.str();

    if(!file_exists(fname)){
        std::stringstream header;
        header << "table_size_byte,probe_size_byte,collisions,thread_count,step_size,algorithm,time_ms,mpps";
        write_to_file(fname, header.str(), true);
    } 
}

void write_summary(
    std::string filename, 
    size_t table_size_byte, 
    size_t probe_size_byte, 
    size_t collisions,
    size_t thread_count, 
    size_t step_size, 
    std::string algorithm, 
    uint64_t time_in_ms, 
    double mpps
){
    std::stringstream c_filename_ss;
    c_filename_ss << filename << "_summary.csv";
    std::string fname = c_filename_ss.str();

    std::stringstream content;
    content << table_size_byte << "," << probe_size_byte << "," << collisions << "," << thread_count << "," << step_size << "," << algorithm << "," << time_in_ms << "," << mpps;
    write_to_file(fname, content.str());
}

void write_raw(
    std::string filename, 
    size_t table_size_byte, 
    size_t probe_size_byte, 
    size_t key_count, 
    size_t collisions, 
    size_t thread_count, 
    size_t step_size, 
    std::string algorithm,
    uint64_t* time_in_ms,
    size_t number_results
){
    std::stringstream c_filename_ss;
    c_filename_ss << filename << "_raw.csv";
    std::string fname = c_filename_ss.str();

    std::stringstream content;
    content << table_size_byte << "," << probe_size_byte << "," << key_count << "," << collisions << "," << thread_count << "," << step_size << "," << algorithm << ",";
    write_to_file(fname, content.str(), time_in_ms, number_results);
}