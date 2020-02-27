#ifndef KMER_H
#define KMER_H

#include <boost/filesystem.hpp>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

namespace rdxon
{

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
	  if (c.minFreq) {
	    uint32_t kcount = boost::lexical_cast<uint32_t>(*tokIter++);
	    uint32_t h1 = boost::lexical_cast<uint32_t>(*tokIter++);
	    uint32_t h2 = boost::lexical_cast<uint32_t>(*tokIter++);
	    if (kcount >= c.minFreq) {
	      bitH1[h1] = true;
	      bitH2[h2] = true;
	    }
	  } else {
	    uint32_t h1 = boost::lexical_cast<uint32_t>(*tokIter++);
	    uint32_t h2 = boost::lexical_cast<uint32_t>(*tokIter++);
	    bitH1[h1] = true;
	    bitH2[h2] = true;
	  }
	}
	++lcount;
	if (lcount % RDXON_CHUNK_SIZE == 0) {
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << lcount << " records." << std::endl;
	}
	//if (lcount > RDXON_CHUNK_SIZE) break;
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


  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmer(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, TMissingKmers& hp) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ parsing." << std::endl;;  
    std::ifstream file(c.infile.string().c_str(), std::ios_base::in | std::ios_base::binary);

    // Singleton masks
    TBitSet singleH1(RDXON_MAX_HASH, false);
    TBitSet singleH2(RDXON_MAX_HASH, false);
    //std::bitset<RDXON_MAX_HASH>& singleH1 = *(new std::bitset<RDXON_MAX_HASH>());
    //std::bitset<RDXON_MAX_HASH>& singleH2 = *(new std::bitset<RDXON_MAX_HASH>());
    
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
      else if ((lcount % 4 == 3) && (avgQual(gline) >= c.minQual)) {
	std::string rcseq(seq);
	reverseComplement(rcseq);
	uint32_t seqlen = seq.size();
	bool filterSeq = true;
	for (uint32_t pos = 0; pos + c.kmerLength <= seqlen; ++pos) {
	  std::string kmerStr = seq.substr(pos, c.kmerLength);
	  if ((nContent(kmerStr)) || (avgQual(kmerStr) < c.minQual)) continue;
	  unsigned h1 = hash_string(kmerStr.c_str());
	  unsigned h2 = hash_string(rcseq.substr(seqlen - c.kmerLength - pos, c.kmerLength).c_str());
	  if (h1 > h2) {
	    unsigned tmp = h1;
	    h1 = h2;
	    h2 = tmp;
	  }
	  if ((!bitH1[h1]) || (!bitH2[h2])) {
	    // K-mer not in DB
	    if ((!singleH1[h1]) || (!singleH2[h2])) {
	      // Potential singleton k-mer due to seq. error
	      singleH1[h1] = true;
	      singleH2[h2] = true;
	    } else {
	      // K-mer not a singleton and not in DB
	      filterSeq = false;
	      typename TMissingKmers::iterator it = hp.find(std::make_pair(h1, h2));
	      if (it == hp.end()) hp.insert(std::make_pair(std::make_pair(h1, h2), 2));
	      else ++it->second;
	    }
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
      //if (lcount > RDXON_CHUNK_SIZE) break;
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
    std::cout << "Filtered reads: " << filterCount << " (" << ((double) filterCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    std::cout << "Passed reads: " << passCount << " (" << ((double) passCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;

    dataIn.pop();
    dataIn.pop();
    file.close();
    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _filterForTheRare(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ filtering." << std::endl;;

    // Dump file
    boost::iostreams::filtering_ostream dumpOut;
    if (c.hasDumpFile) {
      dumpOut.push(boost::iostreams::gzip_compressor());
      dumpOut.push(boost::iostreams::file_sink(c.dumpfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    }
    
    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Data in
    std::ifstream file(c.infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::string header;
    std::string seq;
    uint64_t lcount = 0;
    uint32_t filterCount = 0;
    uint32_t passCount = 0;
    while(std::getline(instream, gline)) {
      if (lcount % 4 == 0) header = gline;
      else if (lcount % 4 == 1) seq = gline;
      else if (lcount % 4 == 2) {} // Skip the spacing line
      else if ((lcount % 4 == 3) && (avgQual(gline) >= c.minQual)) {
	std::string rcseq(seq);
	reverseComplement(rcseq);
	uint32_t seqlen = seq.size();
	bool filterSeq = true;
	for (uint32_t pos = 0; pos + c.kmerLength <= seqlen; ++pos) {
	  std::string kmerStr = seq.substr(pos, c.kmerLength);
	  if ((nContent(kmerStr)) || (avgQual(kmerStr) < c.minQual)) continue;
	  unsigned h1Raw = hash_string(kmerStr.c_str());
	  unsigned h2Raw = hash_string(rcseq.substr(seqlen - c.kmerLength - pos, c.kmerLength).c_str());
	  unsigned h1 = h1Raw;
	  unsigned h2 = h2Raw;
	  if (h1 > h2) {
	    h1 = h2Raw;
	    h2 = h1Raw;
	  }
	  if ((bitH1[h1]) && (bitH2[h2]) && (hs.find(std::make_pair(h1, h2)) != hs.end())) {
	    filterSeq = false;
	    if (c.hasDumpFile) {
	      if (h1Raw < h2Raw) dumpOut << kmerStr.c_str() << std::endl;
	      else dumpOut << rcseq.substr(seqlen - c.kmerLength - pos, c.kmerLength) << std::endl;
	    }
	  }
	}
	if (filterSeq) ++filterCount;
	else {
	  ++passCount;
	  // Output FASTQ
	  dataOut << header << std::endl;
	  dataOut << seq << std::endl;
	  dataOut << "+" << std::endl;
	  dataOut << gline << std::endl;
	}
      }
      ++lcount;
      if (lcount % (RDXON_CHUNK_SIZE * 4) == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
      }
      //if (lcount > RDXON_CHUNK_SIZE) break;
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
    std::cout << "Filtered reads: " << filterCount << " (" << ((double) filterCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    std::cout << "Passed reads: " << passCount << " (" << ((double) passCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    
    // Close input
    dataIn.pop();
    dataIn.pop();
    file.close();

    // Close output
    dataOut.pop();
    dataOut.pop();

    // Close dump file
    if (c.hasDumpFile) {
      dumpOut.pop();
      dumpOut.pop();
    }
    return true;
  }
    
}

#endif
