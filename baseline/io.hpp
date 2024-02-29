#ifndef TUD_VOA_IO_HPP
#define TUD_VOA_IO_HPP

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>

bool file_exists (std::string name);

void write_to_file(std::string filename, std::string content, bool override = false);
void write_to_file(std::string filename, std::string content, uint64_t* result, size_t count, bool override = false);

void create_summary_file(std::string filename);
void create_raw_file(std::string filename);

void write_summary(
    std::string filename, 
    size_t table_size_byte,
    size_t probe_size_byte, 
    size_t collisions, 
    size_t thread_count, 
    size_t step_size, 
    std::string algorithmus, 
    uint64_t time_in_ms, 
    double mpps);

void write_raw(
    std::string filename, 
    size_t table_size_byte, 
    size_t probe_size_byte, 
    size_t key_count, 
    size_t collisions, 
    size_t thread_count, 
    size_t step_size, 
    std::string algorithmus, 
    uint64_t* time_in_ms,
    size_t number_results
);


#endif //TUD_VOA_IO_HPP