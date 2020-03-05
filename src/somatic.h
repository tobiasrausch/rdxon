#ifndef SOMATIC_H
#define SOMATIC_H

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

#include "config.h"
#include "util.h"
#include "kmer.h"


namespace rdxon {

  template<typename TConfigStruct>
  inline int32_t
    somaticRun(TConfigStruct const& c) {
#ifdef PROFILE
    ProfilerStart("rdxon.prof");
#endif

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    
    // Bitmask for filtering
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet bitH1(RDXON_MAX_HASH, false);
    TBitSet bitH2(RDXON_MAX_HASH, false);

    // Germline k-mer map
    typedef std::pair<uint32_t, uint32_t> THashPair;
    typedef std::map<THashPair, uint32_t> THashMap;
    THashMap germHp;
    FilterConfig germC(c, false);
    _fillHashMap(c, bitH1, bitH2, germHp); 
      
    // Clean bit arrays
    for(uint64_t i = 0; i < RDXON_MAX_HASH; ++i) bitH1[i] = false;
    for(uint64_t i = 0; i < RDXON_MAX_HASH; ++i) bitH2[i] = false;


    /*
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
    bool filterRet = false;
    if (c.files.size() == 1) filterRet = _filterForTheRare(c, bitH1, bitH2, hs);
    else filterRet = _filterForTheRarePE(c, bitH1, bitH2, hs);
    if (!filterRet) {
      std::cerr << "Couldn't parse FASTQ files!" << std::endl;
      return 1;
    }
    */

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }


  int somatic(int argc, char **argv) {
    SomaticConfig c;

    // Define generic options
    std::string svtype;
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(30), "min. avg. base quality of k-mer")
      ("recurrence,r", boost::program_options::value<uint16_t>(&c.minOccur)->default_value(3), "min. k-mer recurrence in tumor FASTQ")
      ("maxrecur,s", boost::program_options::value<uint16_t>(&c.maxOccur)->default_value(500), "max. k-mer recurrence in tumor FASTQ")
      ("ctrecurrence,t", boost::program_options::value<uint16_t>(&c.minControlOccur)->default_value(2), "min. k-mer recurrence in control FASTQ")
      ("kmerX,x", boost::program_options::value<boost::filesystem::path>(&c.kmerX), "k-mer.x map file")
      ("kmerY,y", boost::program_options::value<boost::filesystem::path>(&c.kmerY), "k-mer.y map file")
      ("dump,u", boost::program_options::value<boost::filesystem::path>(&c.dumpfile), "gzipped output file for rare k-mers (optional)")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fq.gz"), "output file")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.files), "input file")
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
      std::cout << "Usage:" << std::endl;
      std::cout << " Single-end mode: rdxon " <<  argv[0] << " [OPTIONS] -x <kmer.x.map> -y <kmer.y.map> <tumor.fq.gz> <control.fq.gz>" << std::endl;
      std::cout << " Paired-end mode: rdxon " <<  argv[0] << " [OPTIONS] -x <kmer.x.map> -y <kmer.y.map> -o <outprefix> <tumor.1.fq.gz> <tumor.2.fq.gz> <control.1.fq.gz> <control.2.fq.gz>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }
    
    // Dump file
    if (vm.count("dump")) c.hasDumpFile = true;
    else c.hasDumpFile = false;

    // Check input files
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
	std::cerr << "FASTQ file is missing: " << c.files[file_c].string() << std::endl;
	return 1;
      }
    }
    if ((c.files.size() != 2) && (c.files.size()==4)) {
      std::cerr << "Please specify only 2 FASTQ files (single-end mode, tumor-normal) or 4 FASTQ files (paired-end mode, tumor-normal)!" << std::endl;
      return 1;
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "rdxon ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return somaticRun(c);
  }

}

#endif
