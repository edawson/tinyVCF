#ifndef TINY_VCF_HPP
#define TINY_VCF_HPP
#include "pliib.hpp"
#include "sha1.h"
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
        std::map<std::string, std::string> samples;

        variant(){

        };

	~variant(){
	  //delete [] chrom;
	  //delete [] id;
	  //delete [] ref;
	  //for (auto& i : alt){
	  //  delete [] i;
	  //}

	}

        variant(const variant& h){
            pliib::strcopy(h.chrom, this->chrom);
            this->pos = h.pos;
            pliib::strcopy(h.id, this->id);
            pliib::strcopy(h.ref, this->ref);
            std::vector<char*> a (h.alt.begin(), h.alt.end());
            this->alt = a;

            this->qual = h.qual;

            std::vector<std::string> f (h.filters.begin(), h.filters.end());
            this->filters = f;
            std::map<std::string, std::string> i(h.infos.begin(), h.infos.end());
            this->infos = i;
            std::map<std::string, std::string> s(h.samples.begin(), h.samples.end());
            this->samples = s;
        };

        uint64_t zero_based_position(){
            return pos - 1;
        };

        inline std::string to_vcf(){
            std::stringstream st;
            st << chrom << '\t' << pos << '\t' <<
                id << '\t' << ref << '\t';
            std::vector<std::string> altstrs;
            altstrs.resize(alt.size());
            for (size_t i = 0; i < alt.size(); ++i){
                 
            }
            return st.str();
        };

        inline std::string make_id(){
            Chocobo1::SHA1 s;
            std::stringstream st;
            st << this->chrom << '\n';
            st << this->pos << '\n';
            std::string refcop(this->ref);
            st << pliib::to_upper(refcop) << '\n';
            for (auto& a : this->alt){
                std::string altcop(a);
                st << pliib::to_upper(altcop) << '\n';
            }
            std::string fin = st.str();
            s.addData(fin.c_str(), fin.length());

            return s.finalize().toString();

        };

        inline std::string get_info(std::string infotag){
            if (infos.find(infotag) != infos.end()){
                return infos.at(infotag);
            }
            std::cerr << "INFO TAG NOT FOUND: " << infotag << "." << std::endl;
            return "";
        };

        inline std::uint64_t get_sv_end(){
            // NB: "END" tag is one-based
            // This is amusingly confusing....
            std::string s = get_info("END");
            if (s != ""){
                return std::stoull(s); 
            }
            return 0;
        };
        
        inline std::string get_sv_type(){
            std::string s = get_info("SVTYPE");
            if (s != ""){
                return s;
            }
            return "";
        };

        inline std::uint64_t get_sv_span(int altnum){
            std::string s = get_info("SPAN");
            if (s != ""){
                return std::stoull(s);
            }
            return 0;
        };

    
        // wraps substitutions, indels, and SVS
        // returns the 1-based end position of a variant regardless of type.
        inline std::uint64_t get_reference_end(int altnum){
            std::uint64_t val = this->get_sv_end();
            std::string svtype = get_sv_type();
            // Enforce zero-length SV types for INS variants
            if (val != 0 && svtype == "INS"){
                return 0;
            }
            // IF we don't have an SV end, take the diff between ref and alt.
            // TODO: implement that logic....
            else{
                if (pliib::canonical(this->ref, std::strlen(this->ref)) && 
                        pliib::canonical(this->alt[altnum], std::strlen(this->alt[altnum]))){
                    
                }
            }
            return 0;
        }
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
        var->pos = std::stoull(std::string(splits[1]));
        pliib::strcopy(splits[2], var->id);
        pliib::strcopy(splits[3], var->ref);
        
        char** alt_splits;
        int alt_split_num;
        int* alt_split_sizes;
        pliib::split(splits[4], '\t', alt_splits, alt_split_num, alt_split_sizes);
        var->alt.resize(alt_split_num);
        for (size_t i = 0; i < alt_split_num; ++i){
            pliib::strcopy(alt_splits[i], var->alt[i]);
        }

        pliib::destroy_splits(alt_splits, alt_split_num, alt_split_sizes);

        var->qual = (std::uint8_t) std::atoi(splits[5]);

        char** filter_splits;
        int filter_split_num;
        int* filter_split_sizes;
        
        pliib::split(splits[6], ';', filter_splits, filter_split_num, filter_split_sizes);
        var->filters.resize(filter_split_num);
        for (size_t i = 0; i < filter_split_num; ++i){
            var->filters[i].assign(filter_splits[i]);
        }
        pliib::destroy_splits(filter_splits, filter_split_num, filter_split_sizes);
        
        char** info_splits;
        int info_split_num;
        int* info_split_sizes;
        pliib::split(splits[7], ';', info_splits, info_split_num, info_split_sizes);

        for (size_t i = 0; i < info_split_num; ++i){
            char** i_splits;
            int i_split_num;
            int* i_split_sizes;
            pliib::split(info_splits[i], '=', i_splits, i_split_num, i_split_sizes);
            if (i_split_num == 1){
                var->infos[std::string(i_splits[0])] = std::string(i_splits[0]);
            }
            else if (i_split_num == 2){
                var->infos[std::string(i_splits[0])] = std::string(i_splits[1]);
            }
            else{
                std::cerr << "WARNING: INVALID INFO FIELD LENGTH FOR VARIANT AT " << var->chrom << " " << var->pos << ". IGNORING FIELD." << std::endl;
            }
            pliib::destroy_splits(i_splits, i_split_num, i_split_sizes);
        }

        pliib::destroy_splits(info_splits, info_split_num, info_split_sizes);

        pliib::destroy_splits(splits, num_splits, split_sizes);

    };
}

#endif
