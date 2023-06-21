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
// jbPerf
//----------------------------------------------------------------------
  
jbPerf::jbPerf ()
  : currentEvent_(0),
    numPapiEvents_ (0),
    papiCodes_(),
    papiNames_(),
    papiEventSet_(PAPI_NULL),
    papiCounters_(),
    counters_(),
    isBegun_(false),
    state_(),
    category_(),
    frame_(),
    initialTime_(-1),
    mode_(JB_MODE_ALL - JB_MODE_TRACE)
{
  _TRACE_JBPERF("jbPerf");

  // Create jbPerf directory

  mkdir (JBPERF_DIR,DIR_MODE);

  // Save path to jbPerf directory

  getcwd(path_, MAX_PATH_LENGTH);
  strcat (path_,"/");
  strcat (path_,JBPERF_DIR);

}

//----------------------------------------------------------------------
  
jbPerf::~jbPerf ()
{
  _TRACE_JBPERF("~jbPerf");

#ifdef USE_PAPI
  _CALL_PAPI(PAPI_cleanup_eventset(&papiEventSet_),"cleanup_eventset",PAPI_OK);
  _CALL_PAPI(PAPI_destroy_eventset(&papiEventSet_),"destroy_eventset",PAPI_OK);
#endif
}

//----------------------------------------------------------------------

void jbPerf::init ()
{
  _TRACE_JBPERF("init");

#ifdef USE_PAPI
  _CALL_PAPI(PAPI_library_init(PAPI_VER_CURRENT),"library_init",PAPI_VER_CURRENT);
  _CALL_PAPI(PAPI_set_debug(PAPI_VERB_ECONT),"set_debug",PAPI_OK);
//  _CALL_PAPI(PAPI_thread_init((unsigned long (*)(void))(pthread_self), 0),
//             "thread_init");
  _CALL_PAPI(PAPI_create_eventset(&papiEventSet_),"create_eventset",PAPI_OK);
#endif

#ifdef USE_MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &ip_);
  MPI_Comm_size (MPI_COMM_WORLD, &np_);
#endif
}

//----------------------------------------------------------------------

void jbPerf::finalize ()
{
  _TRACE_JBPERF("finalize");
#ifdef USE_PAPI
#endif
}

//----------------------------------------------------------------------

void jbPerf::mode (unsigned mode)
{
  _TRACE_JBPERF("mode");
  mode_ = mode;
}

//----------------------------------------------------------------------

void jbPerf::event (const char * eventName)
{
  _TRACE_JBPERF1("event",eventName);

#ifdef USE_PAPI
  assert (!isBegun_);
  if (!isBegun_) {

    // Check if event is already inserted

    int i;
    for (i=0; i<numPapiEvents_; i++) {
      if (papiNames_[i] == eventName) {
	fprintf (stderr,"%s:%d: Warning: jbPerf::event() called with "
		 "previously inserted event %s %s\n",__FILE__,__LINE__,
		 papiNames_[i].c_str(),eventName);
	fflush(stderr);
	break;
      }
    }
    if (i==numPapiEvents_) {

      // Update papiNames_

      //   (convert, e.g., "fp-ins" to "PAPI_FP_INS")

      std::string papiName = eventName;
      papiName.replace(papiName.find("-"),1,"_");
      papiName = "PAPI_" + papiName;
      papiNames_.push_back(eventName);

      // Update papiCodes_

      int eventCode;

      _CALL_PAPI(PAPI_event_name_to_code
		 ((char *)(papiName.c_str()),&eventCode),
		 "event_name_to_code",PAPI_OK);
      _CALL_PAPI(PAPI_query_event (eventCode),"query_event",PAPI_OK);
      papiCodes_.push_back(eventCode);

      // update numPapiEvents_

      ++ numPapiEvents_ ;
    }
  } else {
    fprintf (stderr,"%s:%d: Warning: jbPerf::event() called after "
             "event counting started\n",__FILE__,__LINE__);
    fflush(stderr);
  }
#endif
}

//----------------------------------------------------------------------

void jbPerf::begin ()
{
  _TRACE_JBPERF("begin");

  assert (!isBegun_);

  if (!isBegun_) {

#ifdef USE_PAPI
    if (numPapiEvents_ > 0) {
      _CALL_PAPI(PAPI_add_event (&papiEventSet_,papiCodes_[currentEvent_]),
		 "add_event",PAPI_OK);
      _CALL_PAPI(PAPI_start(papiEventSet_),"start",PAPI_OK);
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
}

//----------------------------------------------------------------------

void jbPerf::end ()
{
  _TRACE_JBPERF("end");

  // Stop the top blank region (making sure there's no category in effect)
  category ("");
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
    if (numPapiEvents_ > 0) {
      long long value[1] = {0};
      _CALL_PAPI(PAPI_stop(papiEventSet_, value),"stop",PAPI_OK);
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
}

//----------------------------------------------------------------------

void jbPerf::start (const char * region_base, long long log)
{

  _TRACE_JBPERF1("start",region_base);

  assert (isBegun_);

  std::string region = region_base;

  // Prepend category if needed

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

#ifdef USE_PAPI
  CountersType &p = papiCounters_[region];
  p[currentEvent_] = read_papi();
  c[IND_THIS + IND_VIRT] = PAPI_get_virt_usec();
#endif

  c[IND_THIS + IND_REAL] = get_real();

}

//----------------------------------------------------------------------

void jbPerf::stop (const char * region_base, long long log)
{
  _TRACE_JBPERF1("stop",region_base);
  assert (isBegun_);

  std::string region = region_base;

  // Prepend category if needed

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

#ifdef USE_PAPI
  c[IND_THIS+IND_VIRT] = PAPI_get_virt_usec() - c[IND_THIS+IND_VIRT];
  CountersType &p = papiCounters_[region];
  p[currentEvent_] = read_papi() - p[currentEvent_];
#endif
  c[IND_THIS+IND_REAL] = get_real() - c[IND_THIS+IND_REAL];
  c[IND_THIS+IND_CALLS] = 0 - c[IND_THIS+IND_CALLS];
  c[IND_THIS+IND_LOG] = 0 - c[IND_THIS+IND_LOG];

  // Adjust stored counters using current values

  int indLast;

  indLast = c.size()-2*NUM_EVENTS;
  for (int i=0; i<NUM_EVENTS; i++) {
    c[IND_INCL+indLast+i] += c[IND_THIS+i];
    c[IND_EXCL+indLast+i] += c[IND_THIS+i] + c[IND_ADJUST+i];
    c[IND_ADJUST+i] = 0;
  }

  ++ c[IND_INCL+indLast+IND_CALLS];
  c[IND_INCL+indLast+IND_LOG] += log;

#ifdef USE_PAPI
  indLast = p.size()-2*numPapiEvents_;
  for (int i=0; i<numPapiEvents_; i++) {
    p[indLast+2*i  ] += p[i];                         // Update inclusive
    p[indLast+2*i+1] += p[i] + p[numPapiEvents_ + i]; // Update exclusive
    p[numPapiEvents_ + i] = 0;                        // Clear excl. adjusts
  }
#endif

  // Update parent's exclusive counts

  if (frame_.size() > 0) {
    CountersType &parent = counters_[frame_.top()];
    for (int i=0; i<NUM_EVENTS; i++) {
      parent[IND_ADJUST+i] -= c[IND_THIS+i];
    }
#ifdef USE_PAPI
    CountersType &papiParent = papiCounters_[frame_.top()];
    for (int i=0; i<numPapiEvents_; i++) {
      papiParent[numPapiEvents_+i] -= p[i];
    }
#endif
  }
}

//----------------------------------------------------------------------

void jbPerf::next (const char * region_base, bool useCategory)
{
  _TRACE_JBPERF1("next",region_base);

  if (isBegun_ && (mode_ & JB_MODE_NEXT)) {

    std::string region = region_base;

    // Prepend category if needed

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

#ifdef USE_PAPI

    // Push 0's onto papiCounters vector for next iteration

    CountersType &p = papiCounters_[region];
    p.reserve(p.size() + 2*numPapiEvents_);
    for (int i=0; i<2*numPapiEvents_; i++) {
      p.push_back(0);
    }

#endif

  }
}

//----------------------------------------------------------------------

void jbPerf::advance ()
{
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
}

//----------------------------------------------------------------------

void jbPerf::category (const char * category)
{
  _TRACE_JBPERF1("category",category);

  if (mode_ & JB_MODE_CATEGORY) category_ = category;
}

//----------------------------------------------------------------------

void jbPerf::increment ()
{
  _TRACE_JBPERF("increment");
}

//----------------------------------------------------------------------

void jbPerf::write (bool isLast)
{
  _TRACE_JBPERF("write");

    // Save current path and change to the jbPerf directory

    char pathOrig[MAX_PATH_LENGTH];
    getcwd(pathOrig, MAX_PATH_LENGTH);
    chdir (path_);

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
#ifdef USE_PAPI
	CountersType &p = papiCounters_[region];
#endif

        // Search for end of previous write, if any

	int iBegin = lastWrite_[region] + 2*NUM_EVENTS;
        int iEnd = c.size()-4*NUM_EVENTS;
	// only print last after end() is called
        if (isLast) iEnd+=2*NUM_EVENTS;
	bool isFirst = (lastWrite_[region]==0 && iEnd >= iBegin);

        if (isLast) {
          // Stop region if it's not stopped and it's the last call to write
          if (state_[region] == JB_STATE_STARTED) {
            fprintf (stderr,"%s:%d: Warning: region %s active when write(true)"
                     "called: stopping\n",region.c_str(),__FILE__,__LINE__);
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
          fprintf (fp," user-event");
#ifdef USE_PAPI
          fprintf (fp," time-virt-excl");
          fprintf (fp," time-virt-incl");

          for (int i=0; i<numPapiEvents_; i++) {
            fprintf (fp," %s-excl",papiNames_[i].c_str());
            fprintf (fp," %s-incl",papiNames_[i].c_str());
          }
#endif
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
            fprintf (fp," %lld",c[IND_LOG+i]);  // user-event
#ifdef USE_PAPI
            fprintf (fp," %lld",c[IND_EXCL+IND_VIRT+i]); // time-virt-excl
            fprintf (fp," %lld",c[IND_INCL+IND_VIRT+i]); // time-virt-incl
            for (int k=0; k<numPapiEvents_; k++) {
	      // @@@ MULTIPLE PAPI EVENTS NOT IMPLEMENTED YET!
              fprintf (fp," %lld",p[(i/(NUM_EVENTS))*numPapiEvents_+2*k+1]);
              fprintf (fp," %lld",p[(i/(NUM_EVENTS))*numPapiEvents_+2*k+0]);
            }
#endif
            fprintf (fp,"\n");
            fflush (fp);
          }
        }

        // Mark end of write with END_MARKER

        lastWrite_[region] = iEnd;

        fclose (fp);

      }
    }

  // Change back to the original path

  chdir (pathOrig);

}

//----------------------------------------------------------------------

void jbPerf::allocate (const char * region)
{
  _TRACE_JBPERF("allocate");

  // Allocate timers

  CountersType &c = counters_[region];

  c.reserve(4*NUM_EVENTS); // ADJUST, THIS, FIRST+INCL, FIRST+EXCL
  for (int i=0; i<4*NUM_EVENTS; i++) {
    c.push_back(0);
  }

#ifdef USE_PAPI
  // Allocate PAPI counters

  CountersType &p = papiCounters_[region];
 
  p.reserve(4*numPapiEvents_);
  for (int i=0; i<4*numPapiEvents_; i++) {
    p.push_back(0);
  }

#endif

  
}

//----------------------------------------------------------------------

void jbPerf::write_trace (std::string region_base)
{

  // Save current path and change to the jbPerf directory

  char pathOrig[MAX_PATH_LENGTH];
  getcwd(pathOrig, MAX_PATH_LENGTH);
  chdir (path_);

  // Write the file

  std::string filename = std::string("TRACE");
#ifdef USE_MPI
  // ip is already in category...
  //  char ipstring[10];
  //  sprintf (ipstring,"%d.",ip_);
  //  filename = std::string (ipstring) + filename;
#endif
  if (category_!="") {
    filename = filename + "." + category_;
  }
  FILE *fp = fopen (filename.c_str(),"a");
  fprintf (fp,"%lld %s\n",get_real()-initialTime_,region_base.c_str());
  fclose (fp);

  // Change back to the original path

  chdir (pathOrig);
}

jbPerf perf;

// *********************************************************************
} // namespace jb
// *********************************************************************

static char jbstr[MAX_F77_REGION_LENGTH];

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
HEADER(void XINC_FUN(event) (char *str, int *len))
{
  assert (*len+1 < MAX_F77_REGION_LENGTH);
  for (int i=0; i<*len; i++) jbstr[i] = str[i];
  jbstr[*len] = 0;
  jb::perf.event(jbstr);
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
HEADER(void XINC_FUN(start) (char *str, int *len, int *update))
{
  assert (*len+1 < MAX_F77_REGION_LENGTH);
  for (int i=0; i<*len; i++)  jbstr[i] = str[i];
  jbstr[*len] = 0;
  jb::perf.start(jbstr,*update);
}

/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(stop) (char *str, int *len, int *update))
{
  assert (*len+1 < MAX_F77_REGION_LENGTH);
  for (int i=0; i<*len; i++) jbstr[i] = str[i];
  jbstr[*len] = 0;
  jb::perf.stop(jbstr,*update);
}

/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(next) (char *str, int *len))
{
  assert (*len+1 < MAX_F77_REGION_LENGTH);
  for (int i=0; i<*len; i++) jbstr[i] = str[i];
  jbstr[*len] = 0;
  jb::perf.next(jbstr);
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(advance) ()) 
{ 
  jb::perf.advance(); 
}
/*------------------------------------------------------------------*/
HEADER(void XINC_FUN(category) (char *str, int *len))
{
  assert (*len+1 < MAX_F77_REGION_LENGTH);
  for (int i=0; i<*len; i++) jbstr[i] = str[i];
  jbstr[*len] = 0;
  jb::perf.category(jbstr);
}
/*------------------------------------------------------------------*/
HEADER(void FINC_FUN(write) ()) 
{ 
  jb::perf.write(); 
}

/*------------------------------------------------------------------*/

