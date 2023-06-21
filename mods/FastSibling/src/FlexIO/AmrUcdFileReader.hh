/*
	Code objects with attibutes for behavior on instantiation
	For example
	Array
		fast indexing
		slow copy
	button
		texturemapped
		
*/

#ifndef __AMRUCDFILEREADER_HH_
#define __AMRUCDFILEREADER_HH_
#include <IO.hh>
#include "AmrFileReader.hh"
#include "AmrUcdGridHierarchy.hh"

/*
	It appears that AmrUcdGrid should not be an 
	external data structure.
	
	It should be hidden inside of the AmrGridHierarchy
*/
class AmrUcdFileReader : public AmrFileReader {
protected:
	AmrUcdGridHierarchy genUCD;
	// feed the AmrUcdGrid's to the hierarchy
	// but reclaimation of those grids is questionable.
	
	// However, the amrgridhierarchy should accept a
	// reference to AmrGrid as input.
public:
	AmrUcdFileReader(IObase &f):AmrFileReader(f){
	}
	void getUcd(FlexArray<AmrNode*> &nodes, FlexArray<int> &cells);
};

#endif
