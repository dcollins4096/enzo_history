
//
// GlobalStats constructor and destructor and clear function.
// Clear is separate from the constructor; since the object is declared as static, the constructor
// is only called once.
//

#include <stdio.h>
#include "GlobalStats.h"

GlobalStats::GlobalStats(){
  //fprintf(stderr,"FOO. Construct\n");
}
GlobalStats::~GlobalStats(){
  //fprintf(stderr,"FOO, destruct\n");
}
GlobalStats::Clear(){
