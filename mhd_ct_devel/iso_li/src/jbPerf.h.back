/*****************************************************************************
 *                                                                           *
 * Copyright 2004 James Bordner                                              *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
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
// Functions: + papi (const char *)
//
// Functions: + begin ()
// Functions: + end ()
//
// Functions: + start (const char *)
// Functions: + stop (const char *)
// Functions: + next (const char *)
//
// Functions: + category (const char *, const char *)
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
#  define PAPI_NULL 0
#endif

#include <sys/time.h>
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
    void papi (const char *);     // Insert a PAPI event by name
    void user (const char *);      // Insert a user event by name

    void begin ();                 // Begin PAPI hardware counters
    void end ();                   // End PAPI hardware counters

    //                                Start PAPI hardware counters for a region
    void start (const char *);
    //                                Stop PAPI hardware counters for a region
    void stop (const char *);      
    //                                Start a new counter for a region
    void next (const char *, bool = true); 

    void increment (const char *, long long);
    //                                Update the user-counter for the region

    void advance ();               // Start a new counter for all regions

    void category (const char *, const char *);  
    //                                (re)set a category for all regions

    void write (bool isLast=false);// Dump hardware counters to files

    //--------------------------------------------------------------

  private:

    void allocate (const char * region); // Allocate storage for counters
    void getRegion (std::string &, const char *);
    inline void include_category (std::string &);
    void clear_category ();
    void write_trace (std::string region);  // Write region name to trace file

    inline long long get_real (); // Return real (wallclock) time

    inline long long read_papi ();// Return current PAPI counter value
      
    //--------------------------------------------------------------

  private:

    // PAPI event

    int papiCode_;             // PAPI event code
    std::string papiName_;     // PAPI event name
    std::string PAPIName_;     // PAPI event name
    int papiEventSet_;         // handle to PAPI event set

    // User event

    int numUserEvents_;
    long long userEvent_[MAX_USER_EVENTS];  // counters for user events
    std::map<std::string,int> userIndex_;   // region user event
    std::string userEventName_[MAX_USER_EVENTS]; // inverse userIndex_ mapping

#ifdef OVERHEAD
    long long overhead_start_; // initial time
    long long overhead_write_; // estimate of overhead amount due to write()
    long long overhead_;       // estimate of total overhead amount
#endif

    int call_depth_;           // depth of jbPerf function calls for overhead_

    //                            Vector of counters for each region:
    RegionCountersType counters_;

    //                           indices of last-written set of counters
    std::map<std::string,int> lastWrite_; 

    bool isBegun_;            // Whether we're between begin() and end()
    StateType state_;         // Current state for given region
    std::string category_;    // Current category
    long long initialTime_;   // Mode: see JB_MODE_* defines
    int mode_;                // Saved path to JBPERF directory
    char path_ [MAX_PATH_LENGTH]; // Original path to JBPERF directory
    FrameType frame_;         // Stack frame--required for exclusive counts

    int np_, ip_;
  };

  // *********************************************************************

  inline long long jbPerf::get_real ()
    {

#  ifdef TRACE_PAPI
      printf ("%s:%d PAPI_get_real_usec ()\n",__FILE__,__LINE__);
      fflush(stdout);
#  endif

      struct timeval tv;
      struct timezone tz;
      gettimeofday(&tv,&tz);
      return (long long) (1000000) * tv.tv_sec + tv.tv_usec;

      // The following led to negative values on copper.ncsa.uiuc.edu!
      //      return PAPI_get_real_usec();

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
      long long value = 0;
#ifdef USE_PAPI
      if (papiCode_) {
	int retval;
	_CALL_PAPI(PAPI_read(papiEventSet_,&value),"read",PAPI_OK,retval);
      }
#endif
      return value;
    }

  // *********************************************************************

  extern jbPerf perf;

  // *********************************************************************
} // namespace jb
// *********************************************************************

#endif
