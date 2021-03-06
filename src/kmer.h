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
  inline void
  _chopSeqAndQual(TConfigStruct const& c, std::string const& seq, std::string const& qual, TBitSet const& bitH1, TBitSet const& bitH2, TBitSet& singleH1, TBitSet& singleH2, TMissingKmers& hp) {
    int32_t seqlen = seq.size();
    uint32_t aq = 0;
    uint32_t nsum = 0;
    for (int32_t pos = 0; ((pos < c.kmerLength) && (pos < seqlen)); ++pos) {
      if (seq[pos] == 'N') ++nsum;
      aq += (uint8_t) qual[pos];
    }
    for (int32_t pos = 0; pos + c.kmerLength <= seqlen; ++pos) {
      if (pos) {
	if (seq[pos - 1] == 'N') --nsum;
	if (seq[pos + c.kmerLength - 1] == 'N') ++nsum;
	aq -= (uint8_t) qual[pos - 1];
	aq += (uint8_t) qual[pos + c.kmerLength - 1];
      }
      if ((uint16_t) (aq /c.kmerLength) < c.minQual + 33) continue;
      if (nsum) continue;
      unsigned h1 = 37;
      for(int32_t i = pos; (i < (pos+c.kmerLength)); ++i) h1 = (h1 * 54059) ^ (seq[i] * 76963);
      unsigned h2 = 37;
      for(int32_t i = pos+c.kmerLength-1; i>=(int32_t)pos; --i) h2 = (h2 * 54059) ^ (cpl[(uint8_t) seq[i]] * 76963);
      if (h1 > h2) {
	unsigned tmp = h1;
	h1 = h2;
	h2 = tmp;
      }
      if ((!bitH1[h1]) || (!bitH2[h2])) {
	// K-mer not in DB
	if (((!singleH1[h1]) || (!singleH2[h2])) && (c.minOccur >= 2)) {
	  // Potential singleton k-mer due to seq. error
	  singleH1[h1] = true;
	  singleH2[h2] = true;
	} else {
	  // K-mer not a singleton and not in DB
	  typename TMissingKmers::iterator it = hp.find(std::make_pair(h1, h2));
	  if (it == hp.end()) {
	    if (c.minOccur<2) hp.insert(std::make_pair(std::make_pair(h1, h2), 1));
	    else hp.insert(std::make_pair(std::make_pair(h1, h2), 2));
	  } else ++it->second;
	}
      }
    }
  }
  
  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmerFastaGZ(TConfigStruct const& c, boost::filesystem::path const& infile, TBitSet const& bitH1, TBitSet const& bitH2, TBitSet& singleH1, TBitSet& singleH2, TMissingKmers& hp) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

    // Open file
    std::ifstream file(infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    uint64_t lcount = 0;
    while(std::getline(instream, gline)) {
      if (lcount % 2 == 1) {
	std::string qual(gline.size(), 'E');
	_chopSeqAndQual(c, gline, qual, bitH1, bitH2, singleH1, singleH2, hp);
      }
      ++lcount;
      if (lcount % (RDXON_CHUNK_SIZE * 2) == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
      }
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
    dataIn.pop();
    dataIn.pop();
    file.close();
    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmerFasta(TConfigStruct const& c, boost::filesystem::path const& infile, TBitSet const& bitH1, TBitSet const& bitH2, TBitSet& singleH1, TBitSet& singleH2, TMissingKmers& hp) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

    // Open file
    std::ifstream fafile(infile.string().c_str());
    std::string seq;
    uint64_t lcount = 0;
    if (fafile.good()) {
      std::string gline;
      while(std::getline(fafile, gline)) {
	if (lcount % 2 == 1) {
	  std::string qual(gline.size(), 'E');
	  _chopSeqAndQual(c, gline, qual, bitH1, bitH2, singleH1, singleH2, hp);
	}
	++lcount;
	if (lcount % (RDXON_CHUNK_SIZE * 2) == 0) {
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
	}
      }
      fafile.close();
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
    return true;
  }
  
  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmerFastqGZ(TConfigStruct const& c, boost::filesystem::path const& infile, TBitSet const& bitH1, TBitSet const& bitH2, TBitSet& singleH1, TBitSet& singleH2, TMissingKmers& hp) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

    // Open file
    std::ifstream file(infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::string seq;
    uint64_t lcount = 0;
    while(std::getline(instream, gline)) {
      if (lcount % 4 == 1) seq = gline;
      else if ((lcount % 4 == 3) && (avgQual(gline) >= c.minQual)) {
	_chopSeqAndQual(c, seq, gline, bitH1, bitH2, singleH1, singleH2, hp);
      }
      ++lcount;
      if (lcount % (RDXON_CHUNK_SIZE * 4) == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
      }
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
    dataIn.pop();
    dataIn.pop();
    file.close();
    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmerFastq(TConfigStruct const& c, boost::filesystem::path const& infile, TBitSet const& bitH1, TBitSet const& bitH2, TBitSet& singleH1, TBitSet& singleH2, TMissingKmers& hp) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

    // Open file
    std::ifstream fqfile(infile.string().c_str());
    std::string seq;
    uint64_t lcount = 0;
    if (fqfile.good()) {
      std::string gline;
      while(std::getline(fqfile, gline)) {
	if (lcount % 4 == 1) seq = gline;
	else if ((lcount % 4 == 3) && (avgQual(gline) >= c.minQual)) {
	  _chopSeqAndQual(c, seq, gline, bitH1, bitH2, singleH1, singleH2, hp);
	}
	++lcount;
	if (lcount % (RDXON_CHUNK_SIZE * 4) == 0) {
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
	}
      }
      fqfile.close();
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
    return true;
  }
  
  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmerBAM(TConfigStruct const& c, boost::filesystem::path const& infile, TBitSet const& bitH1, TBitSet const& bitH2, TBitSet& singleH1, TBitSet& singleH2, TMissingKmers& hp) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();

    // Open BAM
    samFile* samfile = sam_open(infile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    uint64_t lcount = 0;

    // Parse BAM
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)) continue;
      
      // Get the read sequence
      uint8_t* qualptr = bam_get_qual(rec);
      uint32_t aq = 0;
      uint32_t nsum = 0;
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) aq += (uint32_t)((uint8_t) qualptr[i]);
      if ((uint16_t) (aq / rec->core.l_qseq) >= c.minQual) {
	std::string seq(rec->core.l_qseq, 'N');
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	nsum = 0;
	aq = 0;
	for (int32_t pos = 0; ((pos < c.kmerLength) && (pos < rec->core.l_qseq)); ++pos) {
	  if (seq[pos] == 'N') ++nsum;
	  aq += (uint8_t) qualptr[pos];
	}
	for (int32_t pos = 0; pos + c.kmerLength <= rec->core.l_qseq; ++pos) {
	  if (pos) {
	    if (seq[pos - 1] == 'N') --nsum;
	    if (seq[pos + c.kmerLength - 1] == 'N') ++nsum;
	    aq -= (uint8_t) qualptr[pos - 1];
	    aq += (uint8_t) qualptr[pos + c.kmerLength - 1];
	  }
	  if ((uint16_t) (aq /c.kmerLength) < c.minQual) continue;
	  if (nsum) continue;
	  unsigned h1 = 37;
	  for(int32_t i = pos; (i < (pos+c.kmerLength)); ++i) h1 = (h1 * 54059) ^ (seq[i] * 76963);
	  unsigned h2 = 37;
	  for(int32_t i = pos+c.kmerLength-1; i>=(int32_t)pos; --i) h2 = (h2 * 54059) ^ (cpl[(uint8_t) seq[i]] * 76963);
	  if (h1 > h2) {
	    unsigned tmp = h1;
	    h1 = h2;
	    h2 = tmp;
	  }
	  if ((!bitH1[h1]) || (!bitH2[h2])) {
	    // K-mer not in DB
	    if (((!singleH1[h1]) || (!singleH2[h2])) && (c.minOccur >= 2)) {
	      // Potential singleton k-mer due to seq. error
	      singleH1[h1] = true;
	      singleH2[h2] = true;
	    } else {
	      // K-mer not a singleton and not in DB
	      typename TMissingKmers::iterator it = hp.find(std::make_pair(h1, h2));
	      if (it == hp.end())
		if (c.minOccur < 2) hp.insert(std::make_pair(std::make_pair(h1, h2), 1));
		else hp.insert(std::make_pair(std::make_pair(h1, h2), 2));
	      else ++it->second;
	    }
	  }
	}
      }
      ++lcount;
      if (lcount % RDXON_CHUNK_SIZE == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << lcount << " reads." << std::endl;
      }
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << lcount << " reads." << std::endl;

    // Clean-up
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    return true;
  }
  
  template<typename TConfigStruct, typename TBitSet, typename TMissingKmers>
  inline bool
  _countMissingKmer(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, TMissingKmers& hp) {
    // Singleton masks
    TBitSet singleH1(RDXON_MAX_HASH, false);
    TBitSet singleH2(RDXON_MAX_HASH, false);
    //std::bitset<RDXON_MAX_HASH>& singleH1 = *(new std::bitset<RDXON_MAX_HASH>());
    //std::bitset<RDXON_MAX_HASH>& singleH2 = *(new std::bitset<RDXON_MAX_HASH>());

    // Parse FASTQs
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      if (!c.intype) {
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] BAM parsing: " << c.files[file_c].string() << std::endl;
	if (!_countMissingKmerBAM(c, c.files[file_c], bitH1, bitH2, singleH1, singleH2, hp)) return false;
      } else if (c.intype == 1) {
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ.gz parsing: " << c.files[file_c].string() << std::endl;
	if (!_countMissingKmerFastqGZ(c, c.files[file_c], bitH1, bitH2, singleH1, singleH2, hp)) return false;
      } else if (c.intype == 2) {
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTA.gz parsing: " << c.files[file_c].string() << std::endl;
	if (!_countMissingKmerFastaGZ(c, c.files[file_c], bitH1, bitH2, singleH1, singleH2, hp)) return false;
      } else if (c.intype == 3) {
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ parsing: " << c.files[file_c].string() << std::endl;
	if (!_countMissingKmerFastq(c, c.files[file_c], bitH1, bitH2, singleH1, singleH2, hp)) return false;
      } else if (c.intype == 4) {
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTA parsing: " << c.files[file_c].string() << std::endl;
	if (!_countMissingKmerFasta(c, c.files[file_c], bitH1, bitH2, singleH1, singleH2, hp)) return false;
      } else {
	std::cerr << "Unsupported file format!" << std::endl;
	return false;
      }
    }
    return true;
  }


  template<typename TConfigStruct, typename TBitSet, typename THashMap>
  inline bool
  _fillHashMap(TConfigStruct const& c, TBitSet& bitH1, TBitSet& bitH2, THashMap& hp) {
    // DB parsing
    if (c.hasKmerTable) {
      if (!_loadKmerDB(c, bitH1, bitH2)) {
	std::cerr << "Couldn't parse k-mer DB!" << std::endl;
	return false;
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
      std::cerr << "Couldn't parse input files!" << std::endl;
      return false;
    }
    return true;
  }



  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _spotifyTheRare(TConfigStruct const& c, std::string const& seq, std::string const& qual, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    int32_t seqlen = seq.size();
    uint32_t aq = 0;
    uint32_t nsum = 0;
    for (int32_t pos = 0; ((pos < c.kmerLength) && (pos < seqlen)); ++pos) {
      if (seq[pos] == 'N') ++nsum;
      aq += (uint8_t) qual[pos];
    }
    for (int32_t pos = 0; pos + c.kmerLength <= seqlen; ++pos) {
      if (pos) {
	if (seq[pos - 1] == 'N') --nsum;
	if (seq[pos + c.kmerLength - 1] == 'N') ++nsum;
	aq -= (uint8_t) qual[pos - 1];
	aq += (uint8_t) qual[pos + c.kmerLength - 1];
      }
      if ((uint16_t) (aq /c.kmerLength) < c.minQual + 33) continue;
      if (nsum) continue;
      unsigned h1 = 37;
      for(int32_t i = pos; (i < (pos+c.kmerLength)); ++i) h1 = (h1 * 54059) ^ (seq[i] * 76963);
      unsigned h2 = 37;
      for(int32_t i = pos+c.kmerLength-1; i>=(int32_t)pos; --i) h2 = (h2 * 54059) ^ (cpl[(uint8_t) seq[i]] * 76963);
      if (h1 > h2) {
	unsigned tmp = h1;
	h1 = h2;
	h2 = tmp;
      }
      if ((bitH1[h1]) && (bitH2[h2]) && (hs.find(std::make_pair(h1, h2)) != hs.end())) return false;
    }
    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _filterForTheRareFasta(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTA filtering: " << c.files[0].string() << std::endl;

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Data in
    std::ifstream fafile(c.files[0].string().c_str());
    std::string header;
    std::string seq;
    uint64_t lcount = 0;
    uint32_t filterCount = 0;
    uint32_t passCount = 0;
    if (fafile.good()) {
      std::string gline;
      while(std::getline(fafile, gline)) {
	if (lcount % 2 == 0) header = gline;
	else if (lcount % 2 == 1) {
	  std::string qual(gline.size(), 'E');
	  bool filterSeq = _spotifyTheRare(c, gline, qual, bitH1, bitH2, hs);
	  if (filterSeq) ++filterCount;
	  else {
	    ++passCount;
	    // Output FASTA
	    dataOut << header << std::endl;
	    dataOut << gline << std::endl;
	  }
	}
	++lcount;
	if (lcount % (RDXON_CHUNK_SIZE * 2) == 0) {
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
	}
      }
      fafile.close();
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
    std::cout << "Filtered reads: " << filterCount << " (" << ((double) filterCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    std::cout << "Passed reads: " << passCount << " (" << ((double) passCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    
    // Close output
    dataOut.pop();
    dataOut.pop();

    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _filterForTheRareFastaGZ(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTA.gz filtering: " << c.files[0].string() << std::endl;

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Data in
    std::ifstream file(c.files[0].string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::string header;
    uint64_t lcount = 0;
    uint32_t filterCount = 0;
    uint32_t passCount = 0;
    while(std::getline(instream, gline)) {
      if (lcount % 2 == 0) header = gline;
      else if (lcount % 2 == 1) {
	std::string qual(gline.size(), 'E');
	bool filterSeq = _spotifyTheRare(c, gline, qual, bitH1, bitH2, hs);
	if (filterSeq) ++filterCount;
	else {
	  ++passCount;
	  // Output FASTA
	  dataOut << header << std::endl;
	  dataOut << gline << std::endl;
	}
      }
      ++lcount;
      if (lcount % (RDXON_CHUNK_SIZE * 2) == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
      }
      //if (lcount > RDXON_CHUNK_SIZE) break;
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 2) << " reads." << std::endl;
    std::cout << "Filtered reads: " << filterCount << " (" << ((double) filterCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    std::cout << "Passed reads: " << passCount << " (" << ((double) passCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    
    // Close input
    dataIn.pop();
    dataIn.pop();
    file.close();

    // Close output
    dataOut.pop();
    dataOut.pop();

    return true;
  }
  

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _filterForTheRareFastq(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ filtering: " << c.files[0].string() << std::endl;

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Data in
    std::ifstream fqfile(c.files[0].string().c_str());
    std::string header;
    std::string seq;
    uint64_t lcount = 0;
    uint32_t filterCount = 0;
    uint32_t passCount = 0;
    if (fqfile.good()) {
      std::string gline;
      while(std::getline(fqfile, gline)) {
	if (lcount % 4 == 0) header = gline;
	else if (lcount % 4 == 1) seq = gline;
	else if (lcount % 4 == 3) {
	  if (avgQual(gline) >= c.minQual) {
	    bool filterSeq = _spotifyTheRare(c, seq, gline, bitH1, bitH2, hs);
	    if (filterSeq) ++filterCount;
	    else {
	      ++passCount;
	      // Output FASTQ
	      dataOut << header << std::endl;
	      dataOut << seq << std::endl;
	      dataOut << "+" << std::endl;
	      dataOut << gline << std::endl;
	    }
	  } else ++filterCount;
	}
	++lcount;
	if (lcount % (RDXON_CHUNK_SIZE * 4) == 0) {
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
	}
      }
      fqfile.close();
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " reads." << std::endl;
    std::cout << "Filtered reads: " << filterCount << " (" << ((double) filterCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    std::cout << "Passed reads: " << passCount << " (" << ((double) passCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    
    // Close output
    dataOut.pop();
    dataOut.pop();

    return true;
  }
  
  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _filterForTheRareFastqGZ(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] FASTQ.gz filtering: " << c.files[0].string() << std::endl;

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Data in
    std::ifstream file(c.files[0].string().c_str(), std::ios_base::in | std::ios_base::binary);
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
      else if (lcount % 4 == 3) {
	if (avgQual(gline) >= c.minQual) {
	  bool filterSeq = _spotifyTheRare(c, seq, gline, bitH1, bitH2, hs);
	  if (filterSeq) ++filterCount;
	  else {
	    ++passCount;
	    // Output FASTQ
	    dataOut << header << std::endl;
	    dataOut << seq << std::endl;
	    dataOut << "+" << std::endl;
	    dataOut << gline << std::endl;
	  }
	} else ++filterCount;
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

    return true;
  }

  template<typename TConfigStruct, typename TBitSet, typename THashSet>
  inline bool
  _filterForTheRareBAM(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] BAM filtering: " << c.files[0].string() << std::endl;

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Data in
    samFile* samfile = sam_open(c.files[0].string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    uint64_t lcount = 0;
    uint32_t filterCount = 0;
    uint32_t passCount = 0;

        // Parse BAM
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)) continue;

      // Get the read sequence
      uint8_t* qualptr = bam_get_qual(rec);
      uint32_t aq = 0;
      uint32_t nsum = 0;
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) aq += (uint32_t)((uint8_t) qualptr[i]);
      if ((uint16_t) (aq / rec->core.l_qseq) >= c.minQual) {
	std::string seq(rec->core.l_qseq, 'N');
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	nsum = 0;
	aq = 0;
	for (int32_t pos = 0; ((pos < c.kmerLength) && (pos < rec->core.l_qseq)); ++pos) {
	  if (seq[pos] == 'N') ++nsum;
	  aq += (uint8_t) qualptr[pos];
	}
	bool filterSeq = true;
	for (int32_t pos = 0; pos + c.kmerLength <= rec->core.l_qseq; ++pos) {
	  if (pos) {
	    if (seq[pos - 1] == 'N') --nsum;
	    if (seq[pos + c.kmerLength - 1] == 'N') ++nsum;
	    aq -= (uint8_t) qualptr[pos - 1];
	    aq += (uint8_t) qualptr[pos + c.kmerLength - 1];
	  }
	  if ((uint16_t) (aq /c.kmerLength) < c.minQual) continue;
	  if (nsum) continue;
	  unsigned h1 = 37;
	  for(int32_t i = pos; (i < (pos+c.kmerLength)); ++i) h1 = (h1 * 54059) ^ (seq[i] * 76963);
	  unsigned h2 = 37;
	  for(int32_t i = pos+c.kmerLength-1; i>=(int32_t)pos; --i) h2 = (h2 * 54059) ^ (cpl[(uint8_t) seq[i]] * 76963);
	  if (h1 > h2) {
	    unsigned tmp = h1;
	    h1 = h2;
	    h2 = tmp;
	  }
	  if ((bitH1[h1]) && (bitH2[h2]) && (hs.find(std::make_pair(h1, h2)) != hs.end())) {
	    filterSeq = false;
	    break;
	  }
	}
	if (filterSeq) ++filterCount;
	else {
	  ++passCount;
	  // Output FASTQ (read1 or read2)
	  std::string rname = bam_get_qname(rec);
	  if (rec->core.flag & BAM_FREAD2) {
	    dataOut << "@" << rname << "R2" << std::endl;
	    if (rec->core.flag & BAM_FREVERSE) {
	      for(int32_t i = rec->core.l_qseq - 1; i>=0; --i) dataOut << cpl[(uint8_t) seq[i]];
	      dataOut << std::endl;
	    } else {
	      dataOut << seq << std::endl;
	    }
	    dataOut << "+" << std::endl;
	    if (rec->core.flag & BAM_FREVERSE) {
	      for (int32_t i = rec->core.l_qseq - 1; i>=0; --i) dataOut << boost::lexical_cast<char>((uint8_t) (qualptr[i] + 33));
	    } else {
	      for (int32_t i = 0; i < rec->core.l_qseq; ++i) dataOut << boost::lexical_cast<char>((uint8_t) (qualptr[i] + 33));
	    }
	    dataOut << std::endl;
	  } else {
	    dataOut << "@" << rname << "R1" << std::endl;
	    if (rec->core.flag & BAM_FREVERSE) {
	      for(int32_t i = rec->core.l_qseq - 1; i>=0; --i) dataOut << cpl[(uint8_t) seq[i]];
	      dataOut << std::endl;
	    } else {
	      dataOut << seq << std::endl;
	    }
	    dataOut << "+" << std::endl;
	    if (rec->core.flag & BAM_FREVERSE) {
	      for (int32_t i = rec->core.l_qseq - 1; i>=0; --i) dataOut << boost::lexical_cast<char>((uint8_t) (qualptr[i] + 33));
	    } else {
	      for (int32_t i = 0; i < rec->core.l_qseq; ++i) dataOut << boost::lexical_cast<char>((uint8_t) (qualptr[i] + 33));
	    } 
	    dataOut << std::endl;
	  }
	}
      } else ++filterCount;
      ++lcount;
      if (lcount % RDXON_CHUNK_SIZE == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << lcount << " reads." << std::endl;
      }
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << lcount << " reads." << std::endl;
    std::cout << "Filtered reads: " << filterCount << " (" << ((double) filterCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    std::cout << "Passed reads: " << passCount << " (" << ((double) passCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    
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
  _filterForTheRarePE(TConfigStruct const& c, TBitSet const& bitH1, TBitSet const& bitH2, THashSet const& hs) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Paired-end FASTQ filtering." << std::endl;;

    // Data out
    boost::filesystem::path fq1 = c.outfile.string() + ".1.fq.gz";
    boost::iostreams::filtering_ostream dataOut1;
    dataOut1.push(boost::iostreams::gzip_compressor());
    dataOut1.push(boost::iostreams::file_sink(fq1.string().c_str(), std::ios_base::out | std::ios_base::binary));
    boost::filesystem::path fq2 = c.outfile.string() + ".2.fq.gz";
    boost::iostreams::filtering_ostream dataOut2;
    dataOut2.push(boost::iostreams::gzip_compressor());
    dataOut2.push(boost::iostreams::file_sink(fq2.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Data in
    std::ifstream file1(c.files[0].string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn1;
    dataIn1.push(boost::iostreams::gzip_decompressor());
    dataIn1.push(file1);
    std::ifstream file2(c.files[1].string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn2;
    dataIn2.push(boost::iostreams::gzip_decompressor());
    dataIn2.push(file2);

    // Streams
    std::istream instream1(&dataIn1);
    std::istream instream2(&dataIn2);
    std::string gline1;
    std::string header1;
    std::string seq1;
    std::string gline2;
    std::string header2;
    std::string seq2;
    uint64_t lcount = 0;
    uint32_t filterCount = 0;
    uint32_t passCount = 0;
    while(true) {
      bool getLine1 = std::getline(instream1, gline1);
      bool getLine2 = std::getline(instream2, gline2);
      if ((!getLine1) || (!getLine2)) break;
      if (lcount % 4 == 0) {
	header1 = gline1;
	header2 = gline2;
      } else if (lcount % 4 == 1) {
	seq1 = gline1;
	seq2 = gline2;
      }
      else if (lcount % 4 == 2) {} // Skip the spacing line
      else if (lcount % 4 == 3) {
	bool filterSeq = true;
	// Read1
	if (avgQual(gline1) >= c.minQual) filterSeq = _spotifyTheRare(c, seq1, gline1, bitH1, bitH2, hs);
	// Read2
	if ((filterSeq) && (avgQual(gline2) >= c.minQual)) filterSeq = _spotifyTheRare(c, seq2, gline2, bitH1, bitH2, hs);
	// Paired-end
	if (filterSeq) ++filterCount;
	else {
	  ++passCount;
	  // Output FASTQs
	  dataOut1 << header1 << std::endl;
	  dataOut1 << seq1 << std::endl;
	  dataOut1 << "+" << std::endl;
	  dataOut1 << gline1 << std::endl;
	  dataOut2 << header2 << std::endl;
	  dataOut2 << seq2 << std::endl;
	  dataOut2 << "+" << std::endl;
	  dataOut2 << gline2 << std::endl;
	}
      }
      ++lcount;
      if (lcount % (RDXON_CHUNK_SIZE * 4) == 0) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " read pairs." << std::endl;
      }
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processed " << (lcount / 4) << " read pairs." << std::endl;
    std::cout << "Filtered pairs: " << filterCount << " (" << ((double) filterCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    std::cout << "Passed pairs: " << passCount << " (" << ((double) passCount * 100.0 / (double) (filterCount + passCount)) << "%)" << std::endl;
    
    // Close input
    dataIn1.pop();
    dataIn1.pop();
    file1.close();
    dataIn2.pop();
    dataIn2.pop();
    file2.close();

    // Close output
    dataOut1.pop();
    dataOut1.pop();
    dataOut2.pop();
    dataOut2.pop();

    return true;
  }


}

#endif
