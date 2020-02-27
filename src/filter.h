#ifndef FILTER_H
#define FILTER_H

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
#include <boost/dynamic_bitset/serialization.hpp>
#include <boost/serialization/bitset.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


#include "util.h"
#include "kmer.h"


namespace rdxon {

  // Config arguments
  struct FilterConfig {
    uint16_t minFreq;
    uint16_t minQual;
    uint16_t minOccur;
    uint16_t maxOccur;
    uint16_t kmerLength;
    bool hasKmerTable;
    bool hasDumpFile;
    boost::filesystem::path outfile;
    boost::filesystem::path kmerdb;
    boost::filesystem::path kmerX;
    boost::filesystem::path kmerY;
    boost::filesystem::path infile;
    boost::filesystem::path dumpfile;
  };


  template<typename TConfigStruct>
  inline int32_t
    rdxonRun(TConfigStruct const& c) {
#ifdef PROFILE
    ProfilerStart("rdxon.prof");
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
    //std::bitset<RDXON_MAX_HASH>& bitH1 = *(new std::bitset<RDXON_MAX_HASH>());
    //std::bitset<RDXON_MAX_HASH>& bitH2 = *(new std::bitset<RDXON_MAX_HASH>());
    
    // Fill hash set
    if (hs.empty()) {
      // K-mer map
      typedef std::map<THashPair, uint32_t> THashMap;
      THashMap hp;  
      
      // DB parsing
      if (c.hasKmerTable) {
	if (!_loadKmerDB(c, bitH1, bitH2)) {
	  std::cerr << "Couldn't parse k-mer DB!" << std::endl;
	  return 1;
	}
	// Save DB
	//std::ofstream outf1("kmer.x.map", std::ios::binary);
	//boost::archive::binary_oarchive outar1(outf1);
	//outar1 << bitH1;
	//outf1.close();
	//std::ofstream outf2("kmer.y.map", std::ios::binary);
	//boost::archive::binary_oarchive outar2(outf2);
	//outar2 << bitH2;
	//outf2.close();
      } else {
	// Load k-mer DB
	std::ifstream inf1(c.kmerX.string().c_str(), std::ios::binary);
	boost::archive::binary_iarchive inar1(inf1);
	inar1 >> bitH1;
	inf1.close();
	std::ifstream inf2(c.kmerY.string().c_str(), std::ios::binary);
	boost::archive::binary_iarchive inar2(inf2);
	inar2 >> bitH2;
	inf2.close();
      }
      
      // Count k-mers not in DB
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
      
      // Debug
      //std::cerr << "Rare k-mer distribution (^RKD)" << std::endl;
      //for(uint32_t i = 0; i < RDXON_KMER_MAXFREQ; ++i) {
      //std::cerr << "RKD\t" << i << "\t" << kmerFreqDist[i] << std::endl;
      //}
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


  int filter(int argc, char **argv) {
    FilterConfig c;

    // Define generic options
    std::string svtype;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(30), "min. avg. base quality of k-mer")
      ("recurrence,r", boost::program_options::value<uint16_t>(&c.minOccur)->default_value(3), "min. k-mer recurrence in FASTQ")
      ("maxrecur,s", boost::program_options::value<uint16_t>(&c.maxOccur)->default_value(500), "max. k-mer recurrence in FASTQ")
      ("kmerX,x", boost::program_options::value<boost::filesystem::path>(&c.kmerX), "k-mer.x map file")
      ("kmerY,y", boost::program_options::value<boost::filesystem::path>(&c.kmerY), "k-mer.y map file")
      ("dump,u", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for rare k-mers (optional)")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fq.gz"), "output file")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input file")
      ("frequency,f", boost::program_options::value<uint16_t>(&c.minFreq)->default_value(1), "min. k-mer frequency in DB [0: two-column input]")
      ("database,d", boost::program_options::value<boost::filesystem::path>(&c.kmerdb), "k-mer database")
      ("kmer,k", boost::program_options::value<uint16_t>(&c.kmerLength)->default_value(61), "k-mer length")
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
    
    // Tabular k-mer file
    if (vm.count("database")) c.hasKmerTable = true;
    else c.hasKmerTable = false;
    
    // Check command line arguments
    bool showHelp = false;
    if (c.hasKmerTable) {
      if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("database"))) showHelp = true;
    } else {
      if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("kmerX")) || (!vm.count("kmerY"))) showHelp = true;
    }
    if (showHelp) {
      std::cout << std::endl;
      std::cout << "Usage: " <<  argv[0] << " [OPTIONS] -x <kmer.x.map> -y <kmer.y.map> <input.fq.gz>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }
    
    // Dump file
    if (vm.count("dump")) c.hasDumpFile = true;
    else c.hasDumpFile = false;
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return rdxonRun(c);
  }

}

#endif
