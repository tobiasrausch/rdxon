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

  static const unsigned char cpl[256] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
    48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
    96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
  };

  inline bool
  is_gz(std::string const& f) {
    std::ifstream bfile(f.c_str(), std::ios_base::binary | std::ios::ate);
    bfile.seekg(0, std::ios::beg);
    char byte1;
    bfile.read(&byte1, 1);
    char byte2;
    bfile.read(&byte2, 1);
    bfile.close();
    if ((byte1 == '\x1F') && (byte2 == '\x8B')) return true;
    else return false;
  }

  inline int32_t     // -1: failure, 0: bam, 1: gzipped fastq, 2: fastq file
  fafqtype(std::string const& path) {
    std::ifstream file(path.c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    int32_t retVal = -1;
    if(std::getline(instream, gline)) {
      if ((gline.size()) && (gline[0] == '@')) retVal = 1; // Gzipped FASTQ
      else if ((gline.size()) && (gline[0] == '>')) retVal = 2; // Gzipped FASTA
    }
    dataIn.pop();
    dataIn.pop();
    file.close();	
    return retVal;
  }

  inline int32_t     // -1: failure, 0: bam, 1: gzipped fastq, 2: fastq file
  inputType(std::string const& path) {
    if (is_gz(path)) {
      return fafqtype(path);
    } else {
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
	else if (fcode[0] == '@') return 3; // Fastq file
	else if (fcode[0] == '>') return 4; // Fasta file
      }
    }
    return -1;
  }

  
  inline uint16_t
  avgQual(std::string const& s) {
    uint32_t aq = 0;
    for(uint32_t i = 0; i < s.size(); ++i) aq += (int32_t) s[i];
    return (uint16_t) ((aq / s.size()) - 33);
  }
  
}

#endif
