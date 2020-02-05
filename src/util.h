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
  #define RDXON_CHUNK_SIZE 1000000
  #endif

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


  template<typename TConfigStruct, typename TBitSet>
  inline bool
  _loadKmerDB(TConfigStruct const& c, TBitSet& bitH1, TBitSet& bitH2) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] K-mer DB parsing." << std::endl;;
    
    // Parse k-mer db
    std::ifstream file(c.kmerdb.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    uint64_t lcount = 0;
    while(std::getline(instream, gline)) {
      if ((gline.size()) && (gline[0] != '#')) {
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(gline, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter != tokens.end()) {
	  uint32_t kcount = boost::lexical_cast<uint32_t>(*tokIter++);
	  uint32_t h1 = boost::lexical_cast<uint32_t>(*tokIter++);
	  uint32_t h2 = boost::lexical_cast<uint32_t>(*tokIter++);
	  if (kcount >= c.minFreq) {
	    bitH1[h1] = true;
	    bitH2[h2] = true;
	  }
	}
	++lcount;
	if (lcount % RDXON_CHUNK_SIZE == 0) {
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << lcount << " records." << std::endl;
	}
	if (lcount > RDXON_CHUNK_SIZE) break;
      }
    }
    dataIn.pop();
    dataIn.pop();
    file.close();
    
    uint64_t bh1c = 0;
    uint64_t bh2c = 0;
    for(uint64_t i = 0; i < RDXON_MAX_HASH; ++i) {
      if (bitH1[i]) ++bh1c;
      if (bitH2[i]) ++bh2c;
    }
    std::cout << "Total flagged k-mers: " << bh1c << ',' << bh2c << std::endl;
    return true;
  }


  template<typename TConfigStruct, typename TBitSet>
  inline bool
  _flagSingletons(TConfigStruct const& c, TBitSet& bitH1, TBitSet& bitH2) {
    // Singleton masks
    TBitSet multiH1(RDXON_MAX_HASH, false);
    TBitSet multiH2(RDXON_MAX_HASH, false);
    //std::bitset<RDXON_MAX_HASH> multiH1;
    //std::bitset<RDXON_MAX_HASH> multiH2;
    
    // Parse FASTQ
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Flag singletons." << std::endl;
    std::ifstream file(c.infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::string seq;
    uint64_t lcount = 0;
    while(std::getline(instream, gline)) {
      if (lcount % 4 == 1) seq = gline;
      else if (lcount % 4 == 3) {
	std::string rcseq(seq);
	reverseComplement(rcseq);
	uint32_t seqlen = seq.size();
	for (uint32_t pos = 0; pos + c.kmerLength <= seqlen; ++pos) {
	  if (nContent(seq.substr(pos, c.kmerLength))) continue;
	  unsigned h1 = hash_string(seq.substr(pos, c.kmerLength).c_str());
	  unsigned h2 = hash_string(rcseq.substr(seqlen - c.kmerLength - pos, c.kmerLength).c_str());
	  if (h1 > h2) {
	    unsigned tmp = h1;
	    h1 = h2;
	    h2 = tmp;
	  }
	  if ((!bitH1[h1]) || (!bitH2[h2])) {
	    bitH1[h1] = true;
	    bitH2[h2] = true;
	  } else {
	    multiH1[h1] = true;
	    multiH2[h2] = true;
	  }
	}
      }
      ++lcount;
      if (lcount % (RDXON_CHUNK_SIZE * 4) == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
      }
      if (lcount > RDXON_CHUNK_SIZE) break;
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;

    // Reset multi k-mers
    uint64_t bh1c = 0;
    uint64_t bh2c = 0;
    for(uint64_t i = 0; i < RDXON_MAX_HASH; ++i) {
      if (multiH1[i]) {
	++bh1c;
	bitH1[i] = false;
      }
      if (multiH2[i]) {
	++bh2c;
	bitH2[i] = false;
      }
    }    
    std::cout << "Multi k-mers: " << bh1c << ',' << bh2c << std::endl;
    bh1c = 0;
    bh2c = 0;
    for(uint64_t i = 0; i < RDXON_MAX_HASH; ++i) {
      if (bitH1[i]) ++bh1c;
      if (bitH2[i]) ++bh2c;
    }
    std::cout << "Singleton k-mers: " << bh1c << ',' << bh2c << std::endl;
    dataIn.pop();
    dataIn.pop();
    file.close();
    return true;
  }
  

  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmer(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, TMissingKmers& hp) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ parsing." << std::endl;;  
    std::ifstream file(c.infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::string seq;
    uint64_t lcount = 0;
    uint32_t filterCount = 0;
    uint32_t passCount = 0;
    while(std::getline(instream, gline)) {
      if (lcount % 4 == 1) seq = gline;
      else if (lcount % 4 == 3) {
	std::string rcseq(seq);
	reverseComplement(rcseq);
	uint32_t seqlen = seq.size();
	bool filterSeq = true;
	for (uint32_t pos = 0; pos + c.kmerLength <= seqlen; ++pos) {
	  if (nContent(seq.substr(pos, c.kmerLength))) continue;
	  unsigned h1 = hash_string(seq.substr(pos, c.kmerLength).c_str());
	  unsigned h2 = hash_string(rcseq.substr(seqlen - c.kmerLength - pos, c.kmerLength).c_str());
	  if (h1 > h2) {
	    unsigned tmp = h1;
	    h1 = h2;
	    h2 = tmp;
	  }
	  if ((!bitH1[h1]) || (!bitH2[h2])) {
	    filterSeq = false;
	    typename TMissingKmers::iterator it = hp.find(std::make_pair(h1, h2));
	    if (it == hp.end()) hp.insert(std::make_pair(std::make_pair(h1, h2), 1));
	    else ++it->second;
	  }
	}
	if (filterSeq) ++filterCount;
	else ++passCount;
      }
      ++lcount;
      if (lcount % (RDXON_CHUNK_SIZE * 4) == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
      }
      if (lcount > RDXON_CHUNK_SIZE) break;
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
    std::cout << "Filtered reads: " << filterCount << std::endl;
    std::cout << "Passed reads: " << passCount << std::endl;

    dataIn.pop();
    dataIn.pop();
    file.close();
    return true;
  }

}

#endif
