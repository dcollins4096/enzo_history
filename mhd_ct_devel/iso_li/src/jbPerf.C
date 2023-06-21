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
//======================================================================
//
// File:        jbPerf.C
//
// Description: Performance tool built on PAPI
//
// See jbPerf.h for list of functions
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================

// #define TRACE_JBPERF
// #define TRACE_STATE
// #define TRACE_PAPI
// #define TRACE_FRAME
#define CHECK_PAPI
// #define DEBUG

#include <assert.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#ifdef USE_MPI
#  include <mpi.h>
#endif
#include "jbPerf.h"

#define DUMMY "PERF"
#define DIR_MODE 0755
#define MAX_EVENT_NAME_LENGTH 20

// *********************************************************************
namespace jb {
// *********************************************************************

//----------------------------------------------------------------------
// jbPerf(): Initialize the jbPerf object
//----------------------------------------------------------------------
  
jbPerf::jbPerf ()
  : papiCode_(PAPI_NULL),
    papiName_(),
    PAPIName_(),
    papiEventSet_(PAPI_NULL),
    numUserEvents_(0),
    userIndex_(),
    userEventName_(),
#ifdef OVERHEAD
    overhead_(0),
    overhead_write_(0),
#endif
    call_depth_(0),
    counters_(),
    isBegun_(false),
    state_(),
    category_(),
    frame_(),
    initialTime_(-1),
    mode_(JB_MODE_ALL - JB_MODE_TRACE - JB_MODE_NEWWRITE)
{
  _TRACE_JBPERF("jbPerf");

  // Clear user counters

  for (int i=0; i<MAX_USER_EVENTS; i++) userEvent_[i] = 0;

  // Create jbPerf directory

  mkdir (JBPERF_DIR,DIR_MODE);

  // Save path to jbPerf directory

  getcwd(path_, MAX_PATH_LENGTH);
  strcat (path_,"/");
  strcat (path_,JBPERF_DIR);
  //  _OVERHEAD_STOP;

}

//----------------------------------------------------------------------
// ~jbPerf(): Finalize the jbPerf object
//----------------------------------------------------------------------
  
jbPerf::~jbPerf ()
{
  _OVERHEAD_START;
  _TRACE_JBPERF("~jbPerf");

#ifdef USE_PAPI
  int retval;
#ifdef PAPI3
  _CALL_PAPI(PAPI_cleanup_eventset(papiEventSet_),"cleanup_eventset",
             PAPI_OK,retval);
#else
  _CALL_PAPI(PAPI_cleanup_eventset(&papiEventSet_),"cleanup_eventset",
             PAPI_OK,retval);
#endif
  _CALL_PAPI(PAPI_destroy_eventset(&papiEventSet_),"destroy_eventset",
             PAPI_OK,retval);
#endif
  _OVERHEAD_STOP;
#ifdef OVERHEAD

  if (ip_==0) {
    long long total_time = get_real() - overhead_start_;
    double total_ratio = (double)(overhead_+overhead_write_)/total_time;
    printf ("jbPerf overhead (us) = %lld (%4.2f%%)\n",
	    overhead_+overhead_write_,100.0*total_ratio);
  }
#endif
}

//----------------------------------------------------------------------
// init(): Initialize PAPI, etc.
//----------------------------------------------------------------------

void jbPerf::init ()
{
  _TRACE_JBPERF("init");

#ifdef USE_PAPI
  int retval;
  _CALL_PAPI(PAPI_library_init(PAPI_VER_CURRENT),"library_init",
	     PAPI_VER_CURRENT,retval);
  _CALL_PAPI(PAPI_set_debug(PAPI_VERB_ECONT),"set_debug",PAPI_OK,retval);
//  _CALL_PAPI(PAPI_thread_init((unsigned long (*)(void))(pthread_self), 0),
//             "thread_init",PAPI_OK,retval),
  _CALL_PAPI(PAPI_create_eventset(&papiEventSet_),"create_eventset",
	     PAPI_OK,retval);
#endif

#ifdef OVERHEAD
  overhead_start_ = get_real();
#endif
  _OVERHEAD_START;

#ifdef USE_MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &ip_);
  MPI_Comm_size (MPI_COMM_WORLD, &np_);
#else
  ip_ = 0;
  np_ = 1;
#endif
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// finalize(): Finalize PAPI, etc.
//----------------------------------------------------------------------

void jbPerf::finalize ()
{
  _OVERHEAD_START;
  _TRACE_JBPERF("finalize");
#ifdef USE_PAPI
#endif
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// mode(): Set global mode 
//----------------------------------------------------------------------

void jbPerf::mode (unsigned mode)
{
  _OVERHEAD_START;
  _TRACE_JBPERF("mode");
  mode_ = mode;
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// papi(): Insert a PAPI event by name
//----------------------------------------------------------------------

void jbPerf::papi (const char * eventName)
{
  _TRACE_JBPERF1("papi",eventName);
  _OVERHEAD_START;
#ifdef USE_PAPI
  assert (!isBegun_);
  if (!isBegun_) {

    // Check if an event is already inserted

    if (papiCode_ != PAPI_NULL) {
      fprintf (stderr,"%s:%d: Warning: jbPerf::papi(%s) called with "
	       "previously inserted event %s\n",__FILE__,__LINE__,
	       eventName,PAPIName_.c_str());
      fflush(stderr);
    } else {

      // Update PAPIName_ (e.g."PAPI_FP_INS") and papiName_ (e.g."fp-ins")

      //   (convert, e.g., "fp-ins" to "PAPI_FP_INS")

      papiName_ = eventName;
      PAPIName_ = eventName;
      PAPIName_.replace(PAPIName_.find("-"),1,"_");
      PAPIName_ = "PAPI_" + PAPIName_;

      // Update papiCode_

      int retval;
      _CALL_PAPI(PAPI_event_name_to_code
		 ((char *)(PAPIName_.c_str()),&papiCode_),
		 "event_name_to_code",PAPI_OK,retval);
      _CALL_PAPI(PAPI_query_event (papiCode_),"query_event",PAPI_OK,retval);

    }
  } else {
    fprintf (stderr,"%s:%d: Warning: jbPerf::papi() called after "
             "event counting started\n",__FILE__,__LINE__);
    fflush(stderr);
  }
#endif
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// user(): Insert a user event by name
//----------------------------------------------------------------------

void jbPerf::user (const char * eventName)
{
  _OVERHEAD_START;
  _TRACE_JBPERF1("user",eventName);

  if (numUserEvents_ < MAX_USER_EVENTS) {

    // INSERT USER EVENT

    userIndex_[eventName] = numUserEvents_;
    userEventName_[numUserEvents_] = eventName;
    ++ numUserEvents_;

  } else {

    // TOO MANY USER EVENTS

    fprintf (stderr,"%s:%d: Warning: jbPerf::event(%s) called with "
	     "more than %d user events\n",__FILE__,__LINE__,
	     eventName,MAX_USER_EVENTS);
    fflush(stderr);
  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// begin(): Begin PAPI hardware counters
//----------------------------------------------------------------------

void jbPerf::begin ()
{
  _OVERHEAD_START;
  _TRACE_JBPERF("begin");

  assert (!isBegun_);

  if (!isBegun_) {

#ifdef USE_PAPI
    if (papiCode_ != PAPI_NULL) {
      int retval;
#ifdef PAPI3
      _CALL_PAPI(PAPI_add_event (papiEventSet_,papiCode_),
		 "add_event",PAPI_OK,retval);
#else
      _CALL_PAPI(PAPI_add_event (&papiEventSet_,papiCode_),
		 "add_event",PAPI_OK,retval);
#endif
      _CALL_PAPI(PAPI_start(papiEventSet_),"start",PAPI_OK,retval);
    }
#endif

    isBegun_ = true;

#ifdef USE_MPI
    // Barrier to attempt to synchronize starting times
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    initialTime_ = get_real();

    // Start a top-level "dummy" region to simplify updating exclusize counts

    start (DUMMY);
    
  } else {
    fprintf (stderr,"%s:%d: Warning: jbPerf::begin() called after "
             "event counting started\n",__FILE__,__LINE__);
    fflush(stderr);
  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// end(): End PAPI hardware counters
//----------------------------------------------------------------------

void jbPerf::end ()
{
  _OVERHEAD_START;
  _TRACE_JBPERF("end");

  // Stop the top blank region (making sure there's no categories in effect)
  clear_category();
  stop (DUMMY);
  if (!frame_.empty()) {
    fprintf (stderr,"%s:%d: Warning: function stack not empty\n"
             "  Possible performance.start(region)/stop(region) mismatch?\n",
             __FILE__,__LINE__);
    fflush(stderr);
  }

  if (isBegun_) {

    // Stop PAPI counters

#ifdef USE_PAPI
    if (papiCode_ != PAPI_NULL) {
      long long value = 0;
      int retval;
      _CALL_PAPI(PAPI_stop(papiEventSet_, &value),"stop",PAPI_OK,retval);
    }
#endif

    // Mark PAPI counters as stopped

    isBegun_ = false;

    // Write un-written events, including last

    write (true);

  } else {
    fprintf (stderr,"%s:%d: Warning: jbPerf::end() called after "
             "event counting stopped\n",__FILE__,__LINE__);
    fflush(stderr);
  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// start(): Start PAPI hardware counters for a region
//----------------------------------------------------------------------

void jbPerf::start (const char * region_base)
{
  _OVERHEAD_START;
  _TRACE_JBPERF1("start",region_base);

  assert (isBegun_);

  std::string region = region_base;

  // Include categories if needed

  include_category (region);

  // Write region name to trace file if JB_MODE_TRACE mode is set

  if (mode_ & JB_MODE_TRACE) write_trace (region);

  // Save region on stack

  frame_.push(region);
#ifdef TRACE_FRAME
  printf ("TRACE_FRAME: push(%s)\n",region.c_str());
#endif

  // Update state

  state_[region] = JB_STATE_STARTED;
  _TRACE_STATE;

  // Adjust stored counters using current values

  CountersType &c = counters_[region];

  // Allocate space if this is the first call to region
  if (c.size()==0) {
    allocate (region.c_str());
    c[IND_BEGIN + IND_START] = get_real();
  }

  // Initialize times and counters

  c[IND_THIS + IND_REAL] = get_real();
#ifdef USE_PAPI
  c[IND_THIS + IND_PAPI] = read_papi();
  c[IND_THIS + IND_VIRT] = PAPI_get_virt_usec();
#endif

  std::map<std::string,int>::iterator pIndex;
  for (pIndex = userIndex_.begin(); 
       pIndex != userIndex_.end(); ++pIndex) {
    int index = (*pIndex).second;
    c[IND_THIS + IND_USER + index] = userEvent_[index];
  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// stop(): Stop PAPI hardware counters for a region
//----------------------------------------------------------------------

void jbPerf::stop (const char * region_base)
{
  _OVERHEAD_START;
  _TRACE_JBPERF1("stop",region_base);
  assert (isBegun_);

  std::string region = region_base;

  // Include categories if needed

  include_category (region);

  // Write region name to trace file if JB_MODE_TRACE mode is set

  if (mode_ & JB_MODE_TRACE) write_trace ("-" + region);

  // Check region against stack

  if (frame_.top() != region) {
    fprintf (stderr,"DEBUG_JBPERF: start(%s) - stop(%s) mismatch\n",
             frame_.top().c_str(),region_base);
    fflush(stderr);
  }

  // Pop region stack 

#ifdef TRACE_FRAME
  printf ("TRACE_FRAME: pop(%s)\n",frame_.top().c_str());
#endif
  frame_.pop();

  // Update state

  assert (state_[region] == JB_STATE_STARTED);
  state_[region] = JB_STATE_STOPPED;
  _TRACE_STATE;

  // Define short aliases for counters

  CountersType &c = counters_[region];

  // Get counter values

  c[IND_THIS+IND_REAL] =           get_real() - c[IND_THIS+IND_REAL];
#ifdef USE_PAPI
  c[IND_THIS+IND_PAPI] =          read_papi() - c[IND_THIS+IND_PAPI];
  c[IND_THIS+IND_VIRT] = PAPI_get_virt_usec() - c[IND_THIS+IND_VIRT];
#endif

  int i;
  for (i=0; i<MAX_USER_EVENTS; i++) {
    c[IND_THIS+IND_USER+i] = userEvent_[i] - c[IND_THIS+IND_USER+i];
  }

  c[IND_THIS+IND_CALLS] =                   0 - c[IND_THIS+IND_CALLS];

  // Adjust stored counters using current values

  int indLast;
  
  indLast = c.size()-2*NUM_EVENTS;
  for (i=0; i<NUM_EVENTS; i++) {
    c[IND_INCL+indLast+i] += c[IND_THIS+i];
    c[IND_EXCL+indLast+i] += c[IND_THIS+i] + c[IND_ADJUST+i];
    c[IND_ADJUST+i] = 0;
  }

  ++ c[IND_INCL+indLast+IND_CALLS];

  // Update parent's exclusive counts

  if (frame_.size() > 0) {
    CountersType &parent = counters_[frame_.top()];
    for (i=0; i<NUM_EVENTS; i++) {
      parent[IND_ADJUST+i] -= c[IND_THIS+i];
    }

  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// next(): Start a new counter for a region
//----------------------------------------------------------------------

void jbPerf::next (const char * region_base, bool useCategory)
{
  _OVERHEAD_START;
  _TRACE_JBPERF1("next",region_base);

  if (isBegun_ && (mode_ & JB_MODE_NEXT)) {

    std::string region = region_base;

    // Include categories if needed

    if (useCategory) include_category (region);

    // Update state

    if (state_[region] != JB_STATE_STOPPED) {
      printf ("error in next(%s)\n",region.c_str());
    }

    assert (state_[region] == JB_STATE_STOPPED);
    state_[region] = JB_STATE_NEW;
    _TRACE_STATE;

    // Push 0's onto counters vector for next iteration

    CountersType &c = counters_[region];
    c.reserve(c.size() + 2*NUM_EVENTS);
    for (int i=0; i<2*NUM_EVENTS; i++) {
      c.push_back(0);
    }
    int indLast = c.size()-2*NUM_EVENTS;
    c[IND_START+indLast] = get_real();
  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// count(): Update the user-counter for the region
//----------------------------------------------------------------------

void jbPerf::increment (const char *counter, 
			long long value)
{
  _OVERHEAD_START;
  _TRACE_JBPERF("increment");

  userEvent_[userIndex_[counter]] += value;
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// advance(): Start a new counter for all regions
//----------------------------------------------------------------------

void jbPerf::advance ()
{
  _OVERHEAD_START;
  _TRACE_JBPERF("advance");

  if (isBegun_ && (mode_ & JB_MODE_ADVANCE)) {

    RegionCountersType::iterator itRegion;

    // Advance each region 
    for (itRegion  = counters_.begin(); 
         itRegion != counters_.end(); ++itRegion) {

      const std::string &region = (*itRegion).first;

      // Call next() for stopped regions

      if (state_[region] == JB_STATE_STOPPED) next (region.c_str(), false);

    }
  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// category(): Set a category for all regions
//----------------------------------------------------------------------

void jbPerf::category (const char * category, const char * value)
{
  _OVERHEAD_START;
  _TRACE_JBPERF1("category",value);

  if (mode_ & JB_MODE_CATEGORY) category_ = value;
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// write(): Dump hardware counters to files
//----------------------------------------------------------------------

void jbPerf::write (bool isLast)
{
  _OVERHEAD_WRITE_START;
  _TRACE_JBPERF("write");

  // Save current path and change to the jbPerf directory

  char pathOrig[MAX_PATH_LENGTH];
  getcwd(pathOrig, MAX_PATH_LENGTH);
  chdir (path_);

  static bool isFirst = true;

  // *** BEGIN NEW WRITE ***

  if (mode_ & JB_MODE_NEWWRITE) {

    static FILE *fp;

    // Open file if first call

    if (isFirst) {
      char filename[5];
      snprintf (filename,5,"%04d",ip_);
      fp = fopen (filename, "w");
    }

    // ======================================
    // Write regions
    // ======================================

    RegionCountersType::iterator pRegion;

    for (pRegion  = counters_.begin(); 
         pRegion != counters_.end(); ++pRegion) {

      std::string region = (*pRegion).first;

      // (skip dummy region)

      if (region != DUMMY) {

	CountersType &c = counters_[region];

        // Search for end of previous write, if any

	int iBegin = lastWrite_[region] + 2*NUM_EVENTS;
        int iEnd = c.size()-4*NUM_EVENTS;
	// only print last after end() is called
        if (isLast) iEnd+=2*NUM_EVENTS;

        if (isLast) {
          // Stop region if it's not stopped and it's the last call to write
          if (state_[region] == JB_STATE_STARTED) {
            fprintf (stderr,"%s:%d: Warning: region %s active when write(true)"
                     "called: stopping\n",__FILE__,__LINE__,region.c_str());
            fflush(stderr);
            stop (region.c_str());
          }
          // Don't print last set of events if they're new
          if (state_[region] == JB_STATE_NEW) {
            iEnd -= NUM_EVENTS;
          }
        }

	printf ("region %s\n",region.c_str());
        for (int i=iBegin; i<= iEnd; i+=2*NUM_EVENTS) {
          if (c[IND_INCL+IND_REAL+i] != 0) {

	printf ("DEBUG %d %d\n",__LINE__,i);
	    fprintf (fp,"region %s\n",region.c_str());
	printf ("DEBUG %d\n",__LINE__);
	    fprintf (fp,"time-stamp %lld\n", 
		     c[IND_INCL+IND_START+i]-initialTime_);
            fprintf (fp,"call-count %lld\n",c[IND_CALLS+i]);
            fprintf (fp,"time-real-excl %lld\n",c[IND_EXCL+IND_REAL+i]);
            fprintf (fp,"time-real-incl %lld\n",c[IND_INCL+IND_REAL+i]);
            fprintf (fp,"time-virt-excl %lld\n",c[IND_EXCL+IND_VIRT+i]);
            fprintf (fp,"time-virt-incl %lld\n",c[IND_INCL+IND_VIRT+i]);
#ifdef USE_PAPI
	    fprintf (fp,"%s-excl %lld\n",
		     papiName_.c_str(),c[IND_EXCL+IND_PAPI+i]);
	    fprintf (fp,"%s-incl %lld\n",
		     papiName_.c_str(),c[IND_INCL+IND_PAPI+i]);
#endif
	    for (int k=0; k<numUserEvents_; k++) {
	      fprintf (fp,"%s %lld\n",userEventName_[k].c_str(),
		       c[IND_EXCL+IND_USER+k+i]);
	      fprintf (fp,"%s %lld\n",userEventName_[k].c_str(),
		       c[IND_INCL+IND_USER+k+i]);
	    }
	    fprintf (fp,"\n");
          }
        }

        // Mark end of write with END_MARKER

        lastWrite_[region] = iEnd;

	// Close file if last call
	if (isLast) {
	  fclose (fp);
	}

      }
    }


  // *** END NEW WRITE ***

  } else {

  // *** BEGIN OLD WRITE ***

    // ======================================
    // Write start()/stop() regions
    // ======================================

    RegionCountersType::iterator pRegion;

    for (pRegion  = counters_.begin(); 
         pRegion != counters_.end(); ++pRegion) {

      std::string region = (*pRegion).first;
      std::string filename = region;

      // (skip dummy region)

      if (region != DUMMY) {

	CountersType &c = counters_[region];

        // Search for end of previous write, if any

	int iBegin = lastWrite_[region] + 2*NUM_EVENTS;
        int iEnd = c.size()-4*NUM_EVENTS;
	// only print last after end() is called
        if (isLast) iEnd+=2*NUM_EVENTS;

        if (isLast) {
          // Stop region if it's not stopped and it's the last call to write
          if (state_[region] == JB_STATE_STARTED) {
            fprintf (stderr,"%s:%d: Warning: region %s active when write(true)"
                     "called: stopping\n",__FILE__,__LINE__,region.c_str());
            fflush(stderr);
            stop (region.c_str());
          }
          // Don't print last set of events if they're new
          if (state_[region] == JB_STATE_NEW) {
            iEnd -= NUM_EVENTS;
          }
        }

        FILE *fp = fopen (filename.c_str(),iBegin ? "a" : "w");

        // Print header line

        if (isFirst) {
          fprintf (fp,"time-start");
          fprintf (fp," time-real-excl");
          fprintf (fp," time-real-incl");
          fprintf (fp," call-count");
#ifdef USE_PAPI
          fprintf (fp," time-virt-excl");
          fprintf (fp," time-virt-incl");

	  fprintf (fp," %s-excl",papiName_.c_str());
	  fprintf (fp," %s-incl",papiName_.c_str());
#endif
	  for (int k=0; k<numUserEvents_; k++) {
	    fprintf (fp," %s-excl",userEventName_[k].c_str());
	    fprintf (fp," %s-incl",userEventName_[k].c_str());
	  }
	    
          fprintf (fp,"\n");
          fflush (fp);
        }

        for (int i=iBegin; i<= iEnd; i+=2*NUM_EVENTS) {
          if (c[IND_INCL+IND_REAL+i] != 0) {
	    //                                      START TIME
            fprintf (fp,"%lld", c[IND_INCL+IND_START+i]-initialTime_);
            fprintf (fp," %lld",c[IND_EXCL+IND_REAL+i]); // time-real-excl
            fprintf (fp," %lld",c[IND_INCL+IND_REAL+i]); // time-real-incl
            fprintf (fp," %lld",c[IND_CALLS+i]); // call-count
#ifdef USE_PAPI
            fprintf (fp," %lld",c[IND_EXCL+IND_VIRT+i]); // time-virt-excl
            fprintf (fp," %lld",c[IND_INCL+IND_VIRT+i]); // time-virt-incl
	    fprintf (fp," %lld",c[IND_EXCL+IND_PAPI+i]);
	    fprintf (fp," %lld",c[IND_INCL+IND_PAPI+i]);
#endif
	    for (int k=0; k<numUserEvents_; k++) {
	      fprintf (fp," %lld",c[IND_EXCL+IND_USER+k+i]);
	      fprintf (fp," %lld",c[IND_INCL+IND_USER+k+i]);
	    }
	    fprintf (fp,"\n");
            fflush (fp);
          }
        }

        // Mark end of write with END_MARKER

        lastWrite_[region] = iEnd;

        fclose (fp);

      }
    }
  // *** END OLD WRITE ***
  }

  // Change back to the original path

  chdir (pathOrig);

  // Update isFirst

  isFirst = false;

  _OVERHEAD_WRITE_STOP;
}

//----------------------------------------------------------------------
// (): 
//----------------------------------------------------------------------

void jbPerf::allocate (const char * region)
{
  _OVERHEAD_START;
  _TRACE_JBPERF("allocate");

  // Allocate timers

  CountersType &c = counters_[region];

  c.reserve(4*NUM_EVENTS); // ADJUST, THIS, FIRST+INCL, FIRST+EXCL
  for (int i=0; i<4*NUM_EVENTS; i++) {
    c.push_back(0);
  }
  _OVERHEAD_STOP;
}

//----------------------------------------------------------------------
// (): 
//----------------------------------------------------------------------

void jbPerf::write_trace (std::string region_base)
{
  _OVERHEAD_START;

  // Save current path and change to the jbPerf directory

  char pathOrig[MAX_PATH_LENGTH];
  getcwd(pathOrig, MAX_PATH_LENGTH);
  chdir (path_);

  // Write the file

  std::string filename = std::string("TRACE");
  if (category_!="") {
    filename = filename + "." + category_;
  }
  FILE *fp = fopen (filename.c_str(),"a");
  fprintf (fp,"%lld %s\n",get_real()-initialTime_,region_base.c_str());
  fclose (fp);

  // Change back to the original path

  chdir (pathOrig);
  _OVERHEAD_STOP;
}

void jbPerf::clear_category()
{
  _OVERHEAD_START;
  category("DEBUG","");
  _OVERHEAD_STOP;
}

void copy_string (char *str1, char *str2, int len)
{
  assert (abs(len+1) < MAX_F77_REGION_LENGTH);
  for (int i=0; i<len; i++) str1[i] = str2[i];
  str1[len] = 0;
}

jbPerf perf;

// *********************************************************************
} // namespace jb
// *********************************************************************

static char cstr1[MAX_F77_REGION_LENGTH];
static char cstr2[MAX_F77_REGION_LENGTH];

/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(init) ()) 
{ 
  jb::perf.init(); 
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(finalize) ()) 
{ 
  jb::perf.finalize(); 
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(mode) (unsigned *mode)) 
{ 
  jb::perf.mode(*mode); 
}
/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(papi) (char *fregion, int *len))
{
  jb::copy_string(cstr1,fregion,*len);
  jb::perf.papi(cstr1);
}
/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(user) (char *fregion, int *len))
{
  jb::copy_string(cstr1,fregion,*len);
  jb::perf.user(cstr1);
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(begin) ()) 
{ 
  jb::perf.begin(); 
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(end) ()) 
{ 
  jb::perf.end(); 
}
/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(start) (char *fregion, int *len))
{
  jb::copy_string(cstr1,fregion,*len);
  jb::perf.start(cstr1);
}

/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(stop) (char *fregion, int *len))
{
  jb::copy_string(cstr1,fregion,*len);
  jb::perf.stop(cstr1);
}

/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(next) (char *fregion, int *len))
{
  jb::copy_string(cstr1,fregion,*len);
  jb::perf.next(cstr1);
}
/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(increment) (char *fevent, int *elen,
				 int *value))
{
  jb::copy_string(cstr1,fevent,*elen);
  jb::perf.increment(cstr1,(long long) *value);
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(advance) ()) 
{ 
  jb::perf.advance(); 
}
/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(category) (char *fcategory, int *clen,
				char *fvalue,    int *vlen))
{
  int i;
  jb::copy_string(cstr1,fcategory,*clen);
  jb::copy_string(cstr2,fvalue,*vlen);
  jb::perf.category(cstr1,cstr2);
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(write) ()) 
{ 
  jb::perf.write(); 
}

/*------------------------------------------------------------------*/

