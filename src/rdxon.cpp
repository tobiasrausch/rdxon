#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>

#define BOOST_DISABLE_ASSERTS

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/dynamic_bitset.hpp>

#include "version.h"
#include "util.h"
#include "kmer.h"


using namespace rdxon;

// Config arguments
struct Config {
  uint16_t minFreq;
  uint16_t minQual;
  uint16_t minOccur;
  uint16_t maxOccur;
  uint16_t kmerLength;
  boost::filesystem::path outfile;
  boost::filesystem::path kmerdb;
  boost::filesystem::path infile;
};


template<typename TConfigStruct>
inline int32_t
rdxonRun(TConfigStruct const& c) {
#ifdef PROFILE
  ProfilerStart("delly.prof");
#endif

  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    
  // K-mer set
  typedef std::pair<uint32_t, uint32_t> THashPair;
  typedef std::set<THashPair> THashSet;
  THashSet hs;

  // Bitmask for filtering
  typedef boost::dynamic_bitset<> TBitSet;
  TBitSet bitH1(RDXON_MAX_HASH, false);
  TBitSet bitH2(RDXON_MAX_HASH, false);

  // Fill hash set
  if (hs.empty()) {
    // K-mer map
    typedef std::map<THashPair, uint32_t> THashMap;
    THashMap hp;  

    // Flag singletons
    if (!_flagSingletons(c, bitH1, bitH2)) {
      std::cerr << "Couldn't parse FASTQ file!" << std::endl;
      return 1;
    }    
    
    // DB parsing
    if (!_loadKmerDB(c, bitH1, bitH2)) {
      std::cerr << "Couldn't parse k-mer DB!" << std::endl;
      return 1;
    }
    
    // Fastq counting step
    if (!_countMissingKmer(c, bitH1, bitH2, hp)) {
      std::cerr << "Couldn't parse FASTQ file!" << std::endl;
      return 1;
    }

    // Clean bit arrays
    for(uint64_t i = 0; i < RDXON_MAX_HASH; ++i) bitH1[i] = false;
    for(uint64_t i = 0; i < RDXON_MAX_HASH; ++i) bitH2[i] = false;
    
    // Process hash table
    std::vector<uint32_t> kmerFreqDist(RDXON_KMER_MAXFREQ, 0);
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Process hash map." << std::endl;;
    uint32_t filterKmerMin = 0;
    uint32_t filterKmerMax = 0;
    uint32_t passKmer = 0;
    for(typename THashMap::iterator it = hp.begin(); it != hp.end(); ++it) {
      if (it->second < RDXON_KMER_MAXFREQ) ++kmerFreqDist[it->second];
      else ++kmerFreqDist[RDXON_KMER_MAXFREQ - 1];
      if (it->second < c.minOccur) ++filterKmerMin;
      else if (it->second > c.maxOccur) ++filterKmerMax;
      else {
	hs.insert(it->first);
	bitH1[it->first.first] = true;
	bitH2[it->first.second] = true;
	++passKmer;
      }
    }
    std::cout << "Filtered k-mers (<min): " << filterKmerMin << std::endl;
    std::cout << "Filtered k-mers (>max): " << filterKmerMax << std::endl;
    std::cout << "Passed hashed k-mers: " << passKmer << std::endl;
    std::cout << "Rare k-mer distribution (^RKD)" << std::endl;
    for(uint32_t i = 0; i < RDXON_KMER_MAXFREQ; ++i) {
      std::cout << "RKD\t" << i << "\t" << kmerFreqDist[i] << std::endl;
    }
  }

  // Filter for the rare
  if (!_filterForTheRare(c, bitH1, bitH2, hs)) {
    std::cerr << "Couldn't parse FASTQ file!" << std::endl;
    return 1;
  }  
  
#ifdef PROFILE
  ProfilerStop();
#endif
  
  // End
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
  return 0;
}


int main(int argc, char **argv) {
  Config c;

  // Define generic options
  std::string svtype;
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("kmer,k", boost::program_options::value<uint16_t>(&c.kmerLength)->default_value(61), "k-mer length")
    ("frequency,f", boost::program_options::value<uint16_t>(&c.minFreq)->default_value(1), "min. k-mer frequency in DB")
    ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(30), "min. avg. base quality of k-mer")
    ("recurrence,r", boost::program_options::value<uint16_t>(&c.minOccur)->default_value(3), "min. k-mer recurrence in FASTQ")
    ("maxrecur,s", boost::program_options::value<uint16_t>(&c.maxOccur)->default_value(500), "max. k-mer recurrence in FASTQ")
    ("database,d", boost::program_options::value<boost::filesystem::path>(&c.kmerdb), "k-mer database")
    ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fq.gz"), "output file")
    ;

  // Define hidden options
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  // Set the visibility
  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);


  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("database"))) { 
    std::cout << std::endl;
    std::cout << "Usage: " <<  argv[0] << " [OPTIONS] -d <kmer.db> <input.fq.gz>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return rdxonRun(c);
}

