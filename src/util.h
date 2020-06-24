#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

namespace rdxon
{

  #ifndef RDXON_MAX_HASH
  #define RDXON_MAX_HASH 4294967296
  #endif

  #ifndef RDXON_CHUNK_SIZE
  #define RDXON_CHUNK_SIZE 10000000
  #endif

  #ifndef RDXON_KMER_MAXFREQ
  #define RDXON_KMER_MAXFREQ 5000
  #endif


  inline int32_t     // -1: failure, 0: bam, 1: gzipped fastq, 2: fastq file
  inputType(std::string const& path) {
    std::ifstream ifile(path.c_str(), std::ios::binary | std::ios::in);
    if (ifile.is_open()) {
      char fcode[4];
      ifile.seekg(0);
      ifile.read(fcode, 4);
      ifile.close();
      bool inputBam = false;
      samFile* samfile = sam_open(path.c_str(), "r");
      if (samfile != NULL) {
	bam_hdr_t* hdr = sam_hdr_read(samfile);
	if (hdr != NULL) {
	  inputBam = true;
	  bam_hdr_destroy(hdr);
	}
	sam_close(samfile);
      }
      if (inputBam) return 0; // Bam
      else if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) return 1; // Gzipped fastq
      else if (fcode[0] == '@') return 2; // Fastq file
      else if (fcode[0] == '>') return 3; // Fasta file
    }
    return -1;
  }

  
  inline void
  reverseComplement(std::string& sequence) {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }

  inline uint16_t
  avgQual(std::string const& s) {
    uint32_t aq = 0;
    for(uint32_t i = 0; i < s.size(); ++i) aq += (int32_t) s[i];
    return (uint16_t) ((aq / s.size()) - 33);
  }

  inline uint16_t
  avgQual(std::vector<uint8_t> const& s) {
    uint32_t aq = 0;
    for(uint32_t i = 0; i < s.size(); ++i) aq += (int32_t) s[i];
    return (uint16_t) (aq / s.size());
  }
  
  inline bool
  nContent(std::string const& s) {
    for(uint32_t i = 0; i < s.size(); ++i) {
      if ((s[i] == 'N') || (s[i] == 'n')) return true;
    }
    return false;
  }

  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }

}

#endif
