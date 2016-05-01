////////////////////////////////////////////////////////////////////////////////
//! \file
//! \author Adam M Phillippy
//! \date 03/26/2003
//!
//! \brief Class declarations for the handling of delta alignment files
//!
//! \see delta.cc
////////////////////////////////////////////////////////////////////////////////

#ifndef __DELTA_HH
#define __DELTA_HH

#include "tigrinc.hh"
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <map>


const std::string NUCMER_STRING = "NUCMER"; //!< nucmer id string
const std::string PROMER_STRING = "PROMER"; //!< promer id string

typedef char AlignmentType_t;               //!< type of alignment data
const AlignmentType_t NULL_DATA = 0;        //!< unknown alignment data type
const AlignmentType_t NUCMER_DATA = 'N';    //!< nucmer alignment data
const AlignmentType_t PROMER_DATA = 'P';    //!< promer alignment data

typedef unsigned char Dir_t;                //!< directional type
const Dir_t FORWARD_DIR = 0;                //!< forward direction
const Dir_t REVERSE_DIR = 1;                //!< reverse direction



//===================================================== DeltaAlignment_t =======
struct DeltaAlignment_t
//!< A single delta encoded alignment region
{
  long sR;    //!< start coordinate in the reference
  long eR;    //!< end coordinate in the reference
  long sQ;    //!< start coordinate in the reference
  long eQ;    //!< end coordinate in the reference
  long idyc;  //!< number of mismatches in the alignment
  long simc;  //!< number of similarity scores < 1 in the alignment
  long stpc;  //!< number of stop codons in the alignment

  float idy;               //!< percent identity [0 - 100]
  float sim;               //!< percent similarity [0 - 100]
  float stp;               //!< percent stop codon [0 - 100]

  std::vector<long> deltas;  //!< delta encoded alignment informaiton

  DeltaAlignment_t ( )
  {
    clear ( );
  }

  void clear ( )
  {
    sR = eR = sQ = eQ = 0;
    idy = sim = stp = 0;
    deltas.clear ( );
  }

  // Read one alignment
  inline bool read_nucmer(std::istream& is, const bool read_deltas = true) { return read(is, false, read_deltas); }
  inline bool read_promer(std::istream& is, const bool read_deltas = true) { return read(is, true, read_deltas); }
  bool read(std::istream& is, const bool promer, const bool read_deltas = true);
};
inline std::istream& operator>>(std::istream& is, DeltaAlignment_t& a) {
  a.read_nucmer(is);
  return is;
}
inline std::ostream& operator<<(std::ostream& os, const DeltaAlignment_t& a) {
  os << a.sR << ' ' << a.eR << ' ' << a.sQ << ' ' << a.eQ << ' ' << a.idyc << ' ' << a.simc << ' ' << a.stpc << '\n';
  for(auto d : a.deltas)
    os << d << '\n';
  return os;
}


//===================================================== DeltaRecord_t ==========
struct DeltaRecord_t
//!< A delta record representing the alignments between two sequences
{
  std::string idR;         //!< reference contig ID
  std::string idQ;         //!< query contig ID
  long lenR;  //!< length of the reference contig
  long lenQ;  //!< length of the query contig

  std::vector<DeltaAlignment_t> aligns; //!< alignments between the two seqs

  DeltaRecord_t ( )
  {
    clear ( );
  }

  void clear ( )
  {
    idR.erase ( );
    idQ.erase ( );
    lenR = lenQ = 0;
    aligns.clear ( );
  }

  bool read(std::istream& is);
};
inline std::istream& operator>>(std::istream& is, DeltaRecord_t& r) {
  r.read(is);
  return is;
}

inline std::ostream& operator<<(std::ostream& os, const DeltaRecord_t& r) {
  return os << '>' << r.idR << ' ' << r.idQ << ' ' << r.lenR << ' ' << r.lenQ;
}

//===================================================== DeltaReader_t ==========
//! \brief Delta encoded file reader class
//!
//! Handles the input of delta encoded alignment information for various MUMmer
//! utilities. Very basic functionality, can be expanded as necessary...
//!
//! \see DeltaRecord_t
//==============================================================================
class DeltaReader_t {

private:

  std::string delta_path_m;      //!< the name of the delta input file
  std::ifstream delta_stream_m;  //!< the delta file input stream
  std::string data_type_m;       //!< the type of alignment data
  std::string reference_path_m;  //!< the name of the reference file
  std::string query_path_m;      //!< the name of the query file
  DeltaRecord_t record_m;        //!< the current delta information record
  bool is_record_m;              //!< there is a valid record in record_m
  bool is_open_m;                //!< delta stream is open


  //--------------------------------------------------- readNextAlignment ------
  //! \brief Reads in a delta encoded alignment
  //!
  //! \param align read info into this structure
  //! \param read_deltas read delta information yes/no
  //! \pre delta_stream is positioned at the beginning of an alignment header
  //! \return void
  //!
  void readNextAlignment (DeltaAlignment_t & align, const bool read_deltas);


  //--------------------------------------------------- readNextRecord ---------
  //! \brief Reads in the next delta record from the delta file
  //!
  //! \param read_deltas read delta information yes/no
  //! \pre delta file must be open
  //! \return true on success, false on EOF
  //!
  bool readNextRecord (const bool read_deltas);


  //--------------------------------------------------- checkStream ------------
  //! \brief Check stream status and abort program if an error has occured
  //!
  //! \return void
  //!
  void checkStream ( )
  {
    if ( !delta_stream_m.good ( ) )
      {
	std::cerr << "ERROR: Could not parse delta file, "
		  << delta_path_m << std::endl;


        std::cerr << "error no: "
                  << int(delta_stream_m.rdstate() & std::ifstream::failbit)
                  << int(delta_stream_m.rdstate() & std::ifstream::badbit)
                  << int(delta_stream_m.rdstate() & std::ifstream::eofbit)
                  << std::endl;
	exit (-1);
      }
  }


public:
  //--------------------------------------------------- DeltaReader_t ----------
  //! \brief Default constructor
  //!
  //! \return void
  //!
  DeltaReader_t ( )
  {
    is_record_m = false;
    is_open_m = false;
  }


  //--------------------------------------------------- ~DeltaReader_t ---------
  //! \brief Default destructor
  //!
  //! \return void
  //!
  ~DeltaReader_t ( )
  {
    close ( );
  }


  //--------------------------------------------------- open -------------------
  //! \brief Opens a delta file by path
  //!
  //! \param delta file to open
  //! \return void
  //!
  void open (const std::string & delta_path);


  //--------------------------------------------------- close ------------------
  //! \brief Closes any open delta file and resets class members
  //!
  //! \return void
  //!
  void close ( )
  {
    delta_path_m.erase ( );
    delta_stream_m.close ( );
    data_type_m.erase ( );
    reference_path_m.erase ( );
    query_path_m.erase ( );
    record_m.clear ( );
    is_record_m = false;
    is_open_m = false;
  }


  //--------------------------------------------------- readNext --------------
  //! \brief Reads in the next delta record from the delta file
  //!
  //! \param read_deltas read delta information yes/no
  //! \pre delta file must be open
  //! \return true on success, false on EOF
  //!
  inline bool readNext (bool getdeltas = true)
  {
    return readNextRecord (getdeltas);
  }



  //--------------------------------------------------- readNextHeadersOnly ----
  //! \brief Reads in the next delta record from the delta file
  //!
  //! Only reads the alignment header information, does not read in the delta
  //! information and leaves each alignment structure's delta vector empty.
  //!
  //! \pre delta file must be open
  //! \return true on success, false on EOF
  //!
  inline bool readNextHeadersOnly ( )
  {
    return readNextRecord (false);
  }


  //--------------------------------------------------- getRecord --------------
  //! \brief Returns a reference to the current delta record
  //!
  //! \pre readNext( ) was successfully called in advance
  //! \return true on success, false on failure or end of file
  //!
  const DeltaRecord_t & getRecord ( ) const
  {
    assert (is_record_m);
    return record_m;
  }


  //--------------------------------------------------- getDeltaPath -----------
  //! \brief Get the path of the current delta file
  //!
  //! \pre delta file is open
  //! \return the path of the delta file
  //!
  const std::string & getDeltaPath ( ) const
  {
    assert (is_open_m);
    return delta_path_m;
  }


  //--------------------------------------------------- getDataType ------------
  //! \brief Get the data type of the current delta file
  //!
  //! \pre delta file is open
  //! \return the data type of the current file
  //!
  const std::string & getDataType ( ) const
  {
    assert (is_open_m);
    return data_type_m;
  }


  //--------------------------------------------------- getReferencePath -------
  //! \brief Get the path of the MUMmer reference file
  //!
  //! \pre delta file is open
  //! \return the path of the MUMmer reference file
  //!
  const std::string & getReferencePath ( ) const
  {
    assert (is_open_m);
    return reference_path_m;
  }


  //--------------------------------------------------- getQueryPath -----------
  //! \brief Get the path of the MUMmer query file
  //!
  //! \pre delta file is open
  //! \return the path of the MUMmer query file
  //!
  const std::string & getQueryPath ( ) const
  {
    assert (is_open_m);
    return query_path_m;
  }
};



//===================================================== SNP_t ==================
struct DeltaEdgelet_t;
struct DeltaEdge_t;
struct DeltaNode_t;
struct SNP_t
     //!< A single nuc/aa poly
{
  long buff;
  char cQ, cR;
  long pQ, pR;
  int conQ, conR;
  std::string ctxQ, ctxR;
  DeltaEdgelet_t * lp;
  DeltaEdge_t * ep;

  SNP_t ( )
  {
    cQ = cR = 0;
    buff = pQ = pR = 0;
    conQ = conR = 0;
  };
};



//===================================================== DeltaEdgelet_t =========
struct DeltaEdgelet_t
//!< A piece of a delta graph edge, a single alignment
{
  unsigned char isGOOD : 1;   //!< meets the requirements
  unsigned char isQLIS : 1;   //!< is part of the query's LIS
  unsigned char isRLIS : 1;   //!< is part of the reference's LIS
  unsigned char isGLIS : 1;   //!< is part of the reference/query LIS
  unsigned char dirR   : 1;   //!< reference match direction
  unsigned char dirQ   : 1;   //!< query match direction

  DeltaEdge_t * edge;
  float idy, sim, stp;        //!< percent identity [0 - 1]
  long idyc, simc, stpc;      //!< idy, sim, stp counts
  long loQ, hiQ, loR, hiR;    //!< alignment bounds
  int frmQ, frmR;             //!< reading frame

  std::string delta;          //!< delta information
  std::vector<SNP_t *> snps;  //!< snps for this edgelet

  DeltaEdgelet_t ( )
  {
    edge = NULL;
    isGOOD = true;
    isQLIS = isRLIS = isGLIS = false;
    dirR = dirQ = FORWARD_DIR;
    idy = sim = stp = 0;
    idyc = simc = stpc = 0;
    loQ = hiQ = loR = hiR = 0;
    frmQ = frmR = 1;
  }

  ~DeltaEdgelet_t ( )
  {
    std::vector<SNP_t *>::iterator i;
    for ( i = snps . begin( ); i != snps . end( ); ++ i )
      delete (*i);
  }

  bool isNegative() const
  { return ( dirR != dirQ ); }

  bool isPositive() const
  { return ( dirR == dirQ ); }

  int slope() const
  { return ( dirR == dirQ ? +1 : -1 ); }

  long loR2Q() const
  { return ( isPositive() ? loQ : hiQ ); }

  long hiR2Q() const
  { return ( isPositive() ? hiQ : loQ ); }

  long loQ2R() const
  { return ( isPositive() ? loR : hiR ); }

  long hiQ2R() const
  { return ( isPositive() ? hiR : loR ); }
};



//===================================================== DeltaEdge_t ============
struct DeltaEdge_t
//!< A delta graph edge, alignments between a single reference and query
{
  DeltaNode_t * refnode;      //!< the adjacent reference node
  DeltaNode_t * qrynode;      //!< the adjacent query node
  std::vector<DeltaEdgelet_t *> edgelets;  //!< the set of individual alignments

  DeltaEdge_t ( )
  { refnode = qrynode = NULL; }

  ~DeltaEdge_t ( )
  {
    std::vector<DeltaEdgelet_t *>::iterator i;
    for ( i = edgelets . begin( ); i != edgelets . end( ); ++ i )
      delete (*i);
  }

  void build (const DeltaRecord_t & rec);
};


//===================================================== DeltaNode_t ============
struct DeltaNode_t
//!< A delta graph node, contains the sequence information
{
  const std::string * id;             //!< the id of the sequence
  char * seq;                         //!< the DNA sequence
  long len;              //!< the length of the sequence
  std::vector<DeltaEdge_t *> edges;   //!< the set of related edges

  DeltaNode_t ( )
  { id = NULL; seq = NULL; len = 0; }

  ~DeltaNode_t ( )
  { free (seq); } // DeltaGraph_t will take care of destructing the edges
};



//===================================================== DeltaGraph_t ===========
//! \brief A graph of sequences (nodes) and their alignments (edges)
//!
//!  A bipartite graph with two partite sets, R and Q, where R is the set of
//!  reference sequences and Q is the set of query sequences. These nodes are
//!  named "DeltaNode_t". We connect a node in R to a node in Q if an alignment
//!  is present between the two sequences. The group of all alignments between
//!  the two is named "DeltaEdge_t" and a single alignment between the two is
//!  named a "DeltaEdgelet_t". Alignment coordinates reference the forward
//!  strand and are stored lo before hi.
//!
//==============================================================================
class DeltaGraph_t
{
public:

  std::map<std::string, DeltaNode_t> refnodes;
  //!< the reference graph nodes, 1 per aligned sequence

  std::map<std::string, DeltaNode_t> qrynodes;
  //!< the query graph nodes, 1 per aligned sequence

  std::string refpath;         //!< path of the reference FastA file
  std::string qrypath;         //!< path of the query FastA file
  AlignmentType_t datatype;    //!< alignment data type

  DeltaGraph_t()
  { datatype = NULL_DATA; }

  ~DeltaGraph_t ( )
  { clear(); }

  void build(const std::string & deltapath, bool getdeltas = true);
  void clean();
  void clear();
  long getNodeCount();
  long getEdgeCount();
  long getEdgeletCount();

  void flagGOOD();
  void flag1to1(float epsilon = -1, float maxolap = 100.0);
  void flagMtoM(float epsilon = -1, float maxolap = 100.0);
  void flagGLIS(float epsilon = -1);
  void flagQLIS(float epsilon = -1,
                float maxolap = 100.0,
                bool flagBad = true);
  void flagRLIS(float epsilon = -1,
                float maxolap = 100.0,
                bool flagbad = true);
  void flagScore(long minlen, float minidy);
  void flagUNIQ(float minuniq);

  void loadSequences();
  std::ostream & outputDelta(std::ostream & out);
};

#endif // #ifndef __DELTA_HH
