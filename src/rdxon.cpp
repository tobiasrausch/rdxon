#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>

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


using namespace rdxon;

// Config arguments
struct Config {
  uint16_t minFreq;
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

  // DB
  typedef boost::dynamic_bitset<> TBitSet;
  TBitSet bitH1(4294967296, false);
  TBitSet bitH2(4294967296, false);
  if (!_loadKmerDB(c, bitH1, bitH2)) {
    std::cerr << "Couldn't parse k-mer DB!" << std::endl;
    return 1;
  }

  // Fastq
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ parsing." << std::endl;;  
  
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
    ("frequency,f", boost::program_options::value<uint16_t>(&c.minFreq)->default_value(1), "min. k-mer frequency in DB")
    ("kmerdb,k", boost::program_options::value<boost::filesystem::path>(&c.kmerdb), "k-mer database")
    ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "output file")
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
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("kmerdb"))) { 
    std::cout << std::endl;
    std::cout << "Usage: " <<  argv[0] << " [OPTIONS] -k <kmer.db> <input.fq.gz>" << std::endl;
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

