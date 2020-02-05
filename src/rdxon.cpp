#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

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
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>


#include "version.h"


using namespace rdxon;

// Config arguments
struct Config {
  uint16_t minMapQual;
  boost::filesystem::path outfile;
  boost::filesystem::path kmerdb;
  boost::filesystem::path infile;
};


template<typename TConfigStruct>
inline int32_t
rdxonRun(TConfigStruct& c) {
#ifdef PROFILE
  ProfilerStart("delly.prof");
#endif
  
  
#ifdef PROFILE
  ProfilerStop();
#endif
  
  // End
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
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
    ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. paired-end (PE) mapping quality")
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

