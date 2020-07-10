#ifndef FETCH_H
#define FETCH_H

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

  template<typename TStream, typename TConfigStruct, typename TBitSet, typename THashSet>
  inline void
  _extractHashedKmer(TStream& dataOut, TConfigStruct const& c, std::string const& seqName, std::string const& seq, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs, TBitSet& obshash) {
    int32_t seqlen = seq.size();
    uint32_t nsum = 0;
    for (int32_t pos = 0; ((pos < c.kmerLength) && (pos < seqlen)); ++pos) {
      if (seq[pos] == 'N') ++nsum;
    }
    for (int32_t pos = 0; pos + c.kmerLength <= seqlen; ++pos) {
      if (pos) {
	if (seq[pos - 1] == 'N') --nsum;
	if (seq[pos + c.kmerLength - 1] == 'N') ++nsum;
      }
      if (nsum) continue;
      unsigned h1Raw = 37;
      for(int32_t i = pos; (i < (pos+c.kmerLength)); ++i) h1Raw = (h1Raw * 54059) ^ (seq[i] * 76963);
      unsigned h2Raw = 37;
      for(int32_t i = pos+c.kmerLength-1; i>=(int32_t)pos; --i) h2Raw = (h2Raw * 54059) ^ (cpl[(uint8_t) seq[i]] * 76963);
      unsigned h1 = h1Raw;
      unsigned h2 = h2Raw;
      if (h1 > h2) {
	h1 = h2Raw;
	h2 = h1Raw;
      }
      if ((bitH1[h1]) && (bitH2[h2])) {
	typename THashSet::const_iterator it = std::lower_bound(hs.begin(), hs.end(), std::make_pair(h1, h2));
	if ((it->first == h1) && (it->second == h2)) {
	  if (c.firstHitOnly) {
	    uint32_t idx = it - hs.begin();
	    if (obshash[idx]) continue;
	    else obshash[idx] = true;
	  }
	  //dataOut << h1 << '\t' << h2 << '\t' << seqName << ':' << (pos + 1) << '-' << (pos + c.kmerLength);
	  dataOut << h1 << '\t' << h2 << '\t' << seqName << ':' << (pos + 1);
	  if (c.includeSeq) {
	    dataOut << '\t';
	    if (h1Raw < h2Raw) {
	      for(int32_t i = pos; (i < (pos+c.kmerLength)); ++i) dataOut << seq[i];
	    } else {
	      for(int32_t i = pos+c.kmerLength-1; i>=(int32_t)pos; --i) dataOut << cpl[(uint8_t) seq[i]];
	    }
	  }
	  dataOut << std::endl;
	}
      }
    }
  }

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _extractHashedKmerBAM(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs, TBitSet& obshash) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] BAM: Extract matching k-mers" << std::endl;

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    // Open BAM
    samFile* samfile = sam_open(c.infile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)) continue;
      
      // Get the read sequence
      std::string seq(rec->core.l_qseq, 'N');
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      _extractHashedKmer(dataOut, c, bam_get_qname(rec), seq, bitH1, bitH2, hs, obshash);
    }
    // Clean-up
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    // Close output
    dataOut.pop();
    dataOut.pop();
    
    return true;
  }

  
  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _extractHashedKmerFasta(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs, TBitSet& obshash) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTA: Extract matching k-mers" << std::endl;
    
    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    // Open file
    std::ifstream fafile(c.infile.string().c_str());
    if (fafile.good()) {
      std::string line;
      std::string faname = "";
      std::string tmpfasta = "";
      while(std::getline(fafile, line)) {
	if (!line.empty()) {
	  if (line[0] == '>') {
	    if (!faname.empty()) {
	      _extractHashedKmer(dataOut, c, faname, tmpfasta, bitH1, bitH2, hs, obshash);
	      // Reset
	      tmpfasta = "";
	      faname = "";
	    }
	    faname = line.substr(1);
	    std::size_t found = faname.find_first_of(" \n\r\t");
	    if (found != std::string::npos) faname.erase(found);
	  } else {
	    if (line.at(line.length() - 1) == '\r' ){
	      tmpfasta += boost::to_upper_copy(line.substr(0, line.length() - 1));
	    } else {
	      tmpfasta += boost::to_upper_copy(line);
	    }
	  }
	}
      }
      _extractHashedKmer(dataOut, c, faname, tmpfasta, bitH1, bitH2, hs, obshash);
      fafile.close();
    }

    // Close output
    dataOut.pop();
    dataOut.pop();
    
    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _extractHashedKmerFastaGZ(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs, TBitSet& obshash) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTA.gz: Extract matching k-mers" << std::endl;

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    // Open file
    std::ifstream file(c.infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string line;
    std::string faname = "";
    std::string tmpfasta = "";
    while(std::getline(instream, line)) {
      if (!line.empty()) {
	if (line[0] == '>') {
	  if (!faname.empty()) {
	    _extractHashedKmer(dataOut, c, faname, tmpfasta, bitH1, bitH2, hs, obshash);
	    // Reset
	    tmpfasta = "";
	    faname = "";
	  }
	  faname = line.substr(1);
	  std::size_t found = faname.find_first_of(" \n\r\t");
	  if (found != std::string::npos) faname.erase(found);
	} else {
	  if (line.at(line.length() - 1) == '\r' ){
	    tmpfasta += boost::to_upper_copy(line.substr(0, line.length() - 1));
	  } else {
	    tmpfasta += boost::to_upper_copy(line);
	  }
	}
      }
    }
    _extractHashedKmer(dataOut, c, faname, tmpfasta, bitH1, bitH2, hs, obshash);
    dataIn.pop();
    dataIn.pop();

    // Close output
    dataOut.pop();
    dataOut.pop();
    
    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _extractHashedKmerFastq(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs, TBitSet& obshash) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ: Extract matching k-mers" << std::endl;
    
    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    // Open file
    std::ifstream fqfile(c.infile.string().c_str());
    if (fqfile.good()) {
      std::string line;
      std::string header;
      uint64_t lcount = 0;
      while(std::getline(fqfile, line)) {
	if (lcount % 4 == 0) {
	  header = line.substr(1);
	  std::size_t found = header.find_first_of(" \n\r\t");
	  if (found != std::string::npos) header.erase(found);
	} else if (lcount % 4 == 1) {
	  _extractHashedKmer(dataOut, c, header, line, bitH1, bitH2, hs, obshash);
	}
	++lcount;
      }
      fqfile.close();
    }

    // Close output
    dataOut.pop();
    dataOut.pop();
    
    return true;
  }


  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _extractHashedKmerFastqGZ(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs, TBitSet& obshash) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ.gz: Extract matching k-mers" << std::endl;
    
    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    // Open file
    std::ifstream file(c.infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string line;
    std::string header;
    uint64_t lcount = 0;
    while(std::getline(instream, line)) {
      if (lcount % 4 == 0) {
	header = line.substr(1);
	std::size_t found = header.find_first_of(" \n\r\t");
	if (found != std::string::npos) header.erase(found);
      } else if (lcount % 4 == 1) {
	_extractHashedKmer(dataOut, c, header, line, bitH1, bitH2, hs, obshash);
      }
      ++lcount;
    }
    dataIn.pop();
    dataIn.pop();

    // Close output
    dataOut.pop();
    dataOut.pop();
    
    return true;
  }

  template<typename TConfigStruct>
  inline int32_t
    fetchRun(TConfigStruct const& c) {
#ifdef PROFILE
    ProfilerStart("rdxon.prof");
#endif

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Read gzipped, sorted hash table." << std::endl;;
    
    // K-mer set
    typedef std::pair<uint32_t, uint32_t> THashPair;
    typedef std::vector<THashPair> THashSet;
    THashSet hs;

    // Bitmask for filtering
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet bitH1(RDXON_MAX_HASH, false);
    TBitSet bitH2(RDXON_MAX_HASH, false);

    // Parse hash table
    if (c.hasHashTable) {
      std::ifstream file(c.htable.string().c_str(), std::ios_base::in | std::ios_base::binary);
      boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
      dataIn.push(boost::iostreams::gzip_decompressor());
      dataIn.push(file);
      std::istream instream(&dataIn);
      std::string gline;
      while(std::getline(instream, gline)) {
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(gline, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter != tokens.end()) {
	  uint32_t h1 = boost::lexical_cast<uint32_t>(*tokIter++);
	  uint32_t h2 = boost::lexical_cast<uint32_t>(*tokIter++);
	  hs.push_back(std::make_pair(h1, h2));
	  bitH1[h1] = true;
	  bitH2[h2] = true;
	}
      }
      dataIn.pop();
      dataIn.pop();
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Loaded " << hs.size() << " hashes." << std::endl;

    // Extract hashed k-mers
    if (!hs.empty()) {
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet obshash;
      if (c.firstHitOnly) obshash.resize(hs.size(), false);
      bool filterRet = false;
      if (!c.intype) filterRet = _extractHashedKmerBAM(c, bitH1, bitH2, hs, obshash);
      else if (c.intype == 1) filterRet = _extractHashedKmerFastqGZ(c, bitH1, bitH2, hs, obshash);
      else if (c.intype == 2) filterRet = _extractHashedKmerFastaGZ(c, bitH1, bitH2, hs, obshash);
      else if (c.intype == 3) filterRet = _extractHashedKmerFastq(c, bitH1, bitH2, hs, obshash);
      else if (c.intype == 4) filterRet = _extractHashedKmerFasta(c, bitH1, bitH2, hs, obshash);
      else {
	std::cerr << "Unsupported file format!" << std::endl;
      }
      if (!filterRet) {
	std::cerr << "Couldn't parse input files!" << std::endl;
	return 1;
      }
    }
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;;
    return 0;
  }

  int fetch(int argc, char **argv) {
    FetchConfig c;

    // Define generic options
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("table,t", boost::program_options::value<boost::filesystem::path>(&c.htable), "gzipped, sorted hash table")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file (only required for BAM input)")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.tsv.gz"), "output file")
      ("sequence,s", "include k-mer sequence in output")
      ("unique,u", "output only first hit for each unique k-mer")
      ;
    
    // Define hidden options
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input file")
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

    // Hash table present?
    if (vm.count("table")) c.hasHashTable = true;
    else c.hasHashTable = false;

    // K-mer sequence?
    if (vm.count("sequence")) c.includeSeq = true;
    else c.includeSeq = false;

    // Unique k-mers only?
    if (vm.count("unique")) c.firstHitOnly = true;
    else c.firstHitOnly = false;
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << std::endl;
      std::cout << "Usage:" << std::endl;
      std::cout << "   rdxon " <<  argv[0] << " [OPTIONS] -t hash.table.gz <input.fa.gz>" << std::endl;
      std::cout << "   rdxon " <<  argv[0] << " [OPTIONS] -t hash.table.gz <input.fq.gz>" << std::endl;
      std::cout << "   rdxon " <<  argv[0] << " [OPTIONS] -g <genome.fa> -t hash.table.gz <input.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 0;
    }
    
    // Check input files
    if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
      std::cerr << "Input file is missing: " << c.infile.string() << std::endl;
      return 1;
    } else {
      c.intype = inputType(c.infile.string());
      if (c.intype == -1) {
	std::cerr << "Unrecognized input file format: " << c.infile.string() << std::endl;
	return 1;
      }
      if (!c.intype) {
	// BAM input
	if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
	  std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
	  return 1;
	}
      }
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "rdxon ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return fetchRun(c);
  }

}

#endif
