#ifndef TINY_VCF_HPP
#define TINY_VCF_HPP
#include "pliib.hpp"
#include <sstream>
#include <vector>
#include <string>
#include <cstdint>
#include <assert.h>


namespace TVCF{


    // TODO: break into multi-allelic and mono-allelic Variant types,
    // where multiallelic just wraps the mono-allelic case.
    struct variant{
        char* chrom;
        std::uint64_t pos;
        char* id;
        char* ref;
        std::vector<char*> alt;
        std::uint8_t qual;
        std::vector<std::string> filters;
        std::map<std::string, std::string> infos;

        inline std::string to_string(){

        };

        inline string make_id(){
            std::stringstream st;
            st << chrom << '_' << pos << ref << '_';
            for (auto a : alt){
                st << a;
            }
        };

    };
    

    inline void parse_line(char* line){

    };

    inline void parse_header_line(char* line){

    };

    inline void parse_variant_line(char* line, variant*& var){
        
        
        char** splits;
        int num_splits;
        int* split_sizes;

        pliib::split(line, '\t', splits, num_splits, split_sizes);
        assert(num_splits >= 7);
        
        pliib::strcopy(splits[0], var->chrom);
        var->pos = stoull(string(splits[1]));
        pliib::strcopy(splits[2], var->id);

    };
}

#endif
