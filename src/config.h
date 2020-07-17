#ifndef CONFIG_H
#define CONFIG_H

namespace rdxon {

  struct SomaticConfig {
    uint16_t minFreq;
    uint16_t minQual;
    uint16_t minOccur;
    uint16_t minControlOccur;
    uint16_t maxOccur;
    uint16_t kmerLength;
    bool hasKmerTable;
    bool hasDumpFile;
    int32_t intype;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
    boost::filesystem::path kmerdb;
    boost::filesystem::path kmerX;
    boost::filesystem::path kmerY;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path dumpfile;
  };

  
  struct FilterConfig {
    uint16_t minFreq;
    uint16_t minQual;
    uint16_t minOccur;
    uint16_t maxOccur;
    uint16_t kmerLength;
    bool hasKmerTable;
    bool hasDumpFile;
    int32_t intype;
    boost::filesystem::path outfile;
    boost::filesystem::path kmerdb;
    boost::filesystem::path kmerX;
    boost::filesystem::path kmerY;
    boost::filesystem::path genome;
    std::vector<boost::filesystem::path> files;
    boost::filesystem::path dumpfile;

    FilterConfig() {}
    
    FilterConfig(SomaticConfig const& sc, bool const germ) {
      minFreq = sc.minFreq;
      minQual = sc.minQual;
      maxOccur = sc.maxOccur;
      kmerLength = sc.kmerLength;
      hasKmerTable = sc.hasKmerTable;
      hasDumpFile = sc.hasDumpFile;
      outfile = sc.outfile;
      genome = sc.genome;
      kmerdb = sc.kmerdb;
      kmerX = sc.kmerX;
      kmerY = sc.kmerY;
      dumpfile = sc.dumpfile;
      intype = sc.intype;
      
      if (germ) {
	// Germline
	minOccur = sc.minControlOccur;
	if (sc.files.size() == 2) {
	  files.push_back(sc.files[1]);
	} else if (sc.files.size() == 4) {
	  files.push_back(sc.files[2]);
	  files.push_back(sc.files[3]);
	}
      } else {
	// Tumor
	minOccur = sc.minOccur;
	if (sc.files.size() == 2) {
	  files.push_back(sc.files[0]);
	} else if (sc.files.size() == 4) {
	  files.push_back(sc.files[0]);
	  files.push_back(sc.files[1]);
	}
      }
    }
  };


  struct FetchConfig {
    uint16_t kmerLength;
    int32_t intype;
    bool hasHashTable;
    bool includeSeq;
    bool firstHitOnly;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
    boost::filesystem::path htable;
    boost::filesystem::path infile;

    FetchConfig() {}    
  };


}

#endif
