#ifndef JBPERF_H
#define JBPERF_H

//======================================================================
//
// File:        jbPerf.h
//
// Description: Simple higher level access to PAPI hardware counters
//
//----------------------------------------------------------------------
//
// Classes:     jbPerf
//
//----------------------------------------------------------------------
//
// Functions: + jbPerf ()
// Functions: + ~jbPerf ()
//
// Functions: + init ()  
//
// Functions: + mode (unsigned)
// Functions: + event (const char *)
//
// Functions: + begin ()
// Functions: + end ()
//
// Functions: + start (const char *)
// Functions: + stop (const char *)
// Functions: + next (const char *)
//
// Functions: + category (const char *)
// Functions: + advance ()
//
// Functions: + write ()
//
// Functions: - allocate (const char *)
// Functions: - getRegion (string, const char *)
// Functions: - include_category (string)
// Functions: - write_trace (string, string)
//
//----------------------------------------------------------------------
//
// Defines: USE_PAPI
// Defines: TRACE_PAPI
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================

#include <stdio.h>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <stack>

#ifdef USE_PAPI
#  include "papi.h"
#else
#  include <sys/time.h>
#  define PAPI_NULL 0
#endif

#include "jbPerf.def"
// *********************************************************************
namespace jb {
  // *********************************************************************

  typedef std::vector<long long> CountersType;
  typedef std::map<std::string,CountersType> RegionCountersType;
  typedef std::stack<std::string> FrameType;
  typedef std::map<std::string,int> StateType;

  //======================================================================
  // jbPerf base class
  //======================================================================

  class jbPerf {

  public:

    jbPerf();                      // Initialize the jbPerf object
    ~jbPerf();                     // Finalize the jbPerf object

    void init ();                  // Initialize PAPI, etc.
    void finalize ();              // Finalize PAPI, etc.

    void mode (unsigned);          // Set global mode 
    void event (const char *);     // Insert a PAPI event by name

    void begin ();                 // Begin PAPI hardware counters
    void end ();                   // End PAPI hardware counters

    //                                Start PAPI hardware counters for a region
    void start (const char *, long long = 0);
    //                                Stop PAPI hardware counters for a region
    void stop (const char *, long long = 0);      
    //                                Start a new counter for a region
    void next (const char *, bool = true); 
    void advance ();               // Start a new counter for all regions

    void category (const char *);  // Set a category for all regions

    void increment ();             // Go to next multiplexed event

    void write (bool isLast=false);// Dump hardware counters to files

    //--------------------------------------------------------------

  private:

    void allocate (const char * region); // Allocate storage for counters
    void getRegion (std::string &, const char *);
    inline void include_category (std::string &);
    void write_trace (std::string region);  // Write region name to trace file

    inline long long get_real (); // Return real (wallclock) time

    inline long long read_papi ();// Return current PAPI counter value
      
    //--------------------------------------------------------------

  private:

    int currentEvent_;         // current PAPI event
    int numPapiEvents_;        // number of events (size of event vectors)
    //                            vector of PAPI event codes
    std::vector<int> papiCodes_; 
    //                            vector of PAPI event names
    std::vector<std::string> papiNames_; 
    RegionCountersType papiCounters_;// Vector of PAPI event counters

    int papiEventSet_;         // handle to PAPI event set

    //                            Vector of counters for each region:
    RegionCountersType counters_;

    //                           index of last-written set of counters
    std::map<std::string,int> lastWrite_; 

    bool isBegun_;            // Whether we're between begin() and end()
    StateType state_;         // Current state for given region
    std::string category_;    // Current category
    long long initialTime_;   // Mode: see JB_MODE_* defines
    int mode_;                // Saved path to JBPERF directory
    char path_ [MAX_PATH_LENGTH]; // Original path to JBPERF directory
    FrameType frame_;         // Stack frame--required for exclusive counts

#ifdef USE_MPI
    int np_, ip_;
#endif
  };

  // *********************************************************************

  inline long long jbPerf::get_real ()
    {
#ifdef USE_PAPI
#  ifdef TRACE_PAPI
      printf ("%s:%d PAPI_get_real_usec ()\n",__FILE__,__LINE__);
      fflush(stdout);
#  endif
      return PAPI_get_real_usec();
#else
      struct timeval tv;
      struct timezone tz;
      gettimeofday(&tv,&tz);
      return (long long) (1000000)*tv.tv_sec + tv.tv_usec;
#endif
    };

  //----------------------------------------------------------------------
  
  inline void jbPerf::include_category (std::string & region)
    {
      if (category_!="") {
	region = category_ + "." + region;
      }
    }

  //----------------------------------------------------------------------

  inline long long jbPerf::read_papi ()
    {
      long long value[1] = {0};
#ifdef USE_PAPI
      if (numPapiEvents_ > 0) {
	_CALL_PAPI(PAPI_read(papiEventSet_,value),"read",PAPI_OK);
      }
#endif
      return value[0];
    }

  // *********************************************************************

  extern jbPerf perf;

  // *********************************************************************
} // namespace jb
// *********************************************************************

#endif
