#include <stdio.h>
#include <stdlib.h>
#include "AmrUcdGridHierarchy.hh"

void AmrUcdGridHierarchy::AmrUcdGrid::copyFromAmrGrid(AmrGrid &g){
  level = g.level;
  dataveclen=g.dataveclen;
  rank = g.rank;
  datatype = IObase::Int2DataType(g.datatype);
  // Perhaps should do floating point only
  vdata = g.data;
  // copy to float datatype
  printf("AmrUcdGridHierarchy:: dataveclen=%u\n",dataveclen);
  if(g.data){
    int nelements = IObase::nElements(g.rank,g.dims);
    register int i;
    //if(datatype!=IObase::Float32)
      data = new float[nelements];
    switch(datatype){
    case IObase::Float32:
      printf("AmrUcdGridHierarchy:: Original data is float use directly\n");
      printf("\tpointerbase=%lu last=%lu\n",data,(data+nelements));
	     //data=(float*)g.data;
	     for(i=0;i<nelements;i++) data[i]=(float)((float*)(g.data))[i];
      break;
    case IObase::Float64:
      printf("AmrUcdGridHierarchy:: Original data is double... copy\n");
      for(i=0;i<nelements;i++) data[i]=(float)((double*)(g.data))[i];
      break;
    case IObase::Int32: 
    case IObase::uInt32:
      for(i=0;i<nelements;i++) data[i]=(float)((int*)(g.data))[i];
      break;
    case IObase::Int16:
    case IObase::uInt16:
      for(i=0;i<nelements;i++) data[i]=(float)((short*)(g.data))[i];
      break;
    case IObase::Int8:
    case IObase::uInt8:
    default:
      for(i=0;i<nelements;i++) data[i]=(float)((char*)(g.data))[i];
      break;
    };
  }
  else
    data=0; // danger
  register int i;
  for(i=0;i<3;i++){
    dims[i]=g.dims[i];
    origin[i]=g.origin[i];
    delta[i]=g.delta[i];
  }
  bounds.setFromExtents(g.minext,g.maxext);
  for(i=0,nelements=1;i<rank;i++) nelements*=dims[i];
}

void AmrUcdGridHierarchy::AmrUcdGrid::regenPoints(){
  int i;
  //bounds.setFromOriginDx(origin,delta,dims); // already set from ext
  // should set from extents
  // create points
  points.setSize(nelements);
  points.setSize(nelements);
  // for each point, lets initialize
  long index,dataindex;
  int mindex[3]={0,0,0};
  //----------Node Initialization Cycle-----------------
  printf("AmrUcdGridHierarchy:: dataveclen=%u\n",dataveclen);
  for(index=0,dataindex=0;index<nelements;index++,dataindex+=dataveclen,(mindex[0])++){
    AmrNode *n = points.getData()+index;
    for(i=0;mindex[i]>=dims[i] && i<2;i++){
      mindex[i]=0;
      (mindex[i+1])++;
    }
    // assign using correct type (could put outside of loop to optimize)
    n->data = data+dataindex;
    n->level = level; // semiredunant, but necessary for AVS
    /*
    if(datatype==IObase::Float32)
      n->fdata=fdata+dataindex;
    else if(IObase::Float64)
      n->ddata=ddata+dataindex;
    else if(IObase::Int32 || IObase::uInt32)
      n->idata=idata+dataindex;
    else if(IObase::Int16 || IObase::uInt16)
      n->sdata=sdata+dataindex;
    else if(IObase::Int8 || IObase::uInt8)
      n->cdata=cdata+dataindex;
      */
    n->child=n->parent=0; // for now (will fix this up later)
    n->index=-1; // nil index for now
    n->gridID = gridID;
    for(i=0;i<3;i++)
      n->location.array[i] = origin[i] + delta[i]*(double)(mindex[i]);
    // assumes f77 order for arrays
  }
}


AmrUcdGridHierarchy::~AmrUcdGridHierarchy() { purge(); }

void AmrUcdGridHierarchy::purge(){
  for(int l=0;l<grids.getSize();l++){
    AmrUcdGrid *grid;
    FOREACH(grid,grids[l]){
      delete grid;
    }
  }
  grids.purge();
  levelinfo.purge();
}

void AmrUcdGridHierarchy::addGrid(AmrGrid &g){
  if(grids.getSize()<g.level+1)
    grids.setSize(g.level+1);
  AmrUcdGrid *grid,*parent;
  grid = new AmrUcdGrid(g);
  //--------------Grid Parenting Cycle----------------
  if(g.level>0){ // make grid hierarchy (eg. assign parents) (perhaps doubly-link)
    // compare the bounding box to all grids in parent level
    FOREACH(parent,grids[g.level-1]){ // OrderedList needs copy constructor
      if(parent->bounds.contains(grid->bounds)){  // compared 
		grid->parent=parent;
		parent->children.append(grid);// do forward ref
		break;
      }
    }
  }
  grids[g.level].add(grid);
}

void AmrUcdGridHierarchy::buildNodeHierarchy(){
  int i;
  levelinfo.setSize(grids.getSize());
  for(i=0;i<levelinfo.getSize();i++){
    grids[i].reset();  // restart
    AmrUcdGrid *g = grids[i].getNext();
    if(g)
      levelinfo[i].setDelta(g->delta);
    else{
      double fakedelta[5]={1,1,1,1,1};
      printf("AmrUcdGridHierarchy::buildNodeHierarchy():  Using Fake delta for level %u (shouldn't do that)\n",i);
      levelinfo[i].setDelta(fakedelta);
    }
  }
  for(i=0;i<(levelinfo.getSize()-1);i++)
    levelinfo[i].setStride(levelinfo[i+1].delta);
  // cycle through each level and find parent nodes for each node
  //-------------Node Parenting Cycle----------------------------
  puts("node parenting cycle");
  for(i=1;i<levelinfo.getSize();i++){
    AmrUcdGrid *grid;
    printf("\tL%u\n",i);
    FOREACH(grid,grids[i]){ // level 2 grids should not be here!!!
      // First find offset
      int offset[3]; // we'll skip that for now and assume aligned
      int origin[3]; // origin with respect to parent
      int stride[3]; // just here for convenience
      // now stride through the parent grid doing parenting replacement
      AmrUcdGrid *parent=grid->parent;
      if(!parent) continue; // failsafe for NULL parent
      // stride through the parent doing replacement
      // first compute the strides through the parent grid
      // as well as the origin in parent coordinates
      for(int j=0;j<3;j++){
	register double dp=parent->origin[j];
	register double dc=grid->origin[j];
	register double dx=dc-dp; // difference between origins
	dx /= parent->delta[j]; // if diff > .5 then use it? 
	if(j==0) printf("dparent=%lf dchild=%lf dx=%lf\n",dp,dc,dx);
	// (no just need to round it)
	origin[j]=(int)(dx+0.5); // this is a marginally bogus way to round up.
	stride[j]=levelinfo[i-1].childstride[j];
      }
      printf("NodeParenting L%u: origin=%u:%u:%u stride=%u\n",i,origin[0],origin[1],origin[2],stride[0]);
      int px,x,py,y,pz,z;
      // typedef DATATYPE AmrNode;
      FlexArray<AmrNode> *pnodes=&(parent->points);
      FlexArray<AmrNode> *nodes=&(grid->points);
      for(z=0,pz=origin[2];z<grid->dims[2];z+=stride[2],pz++){
	int zi=grid->dims[1]*z;
	int pzi=parent->dims[1]*pz;
	for(y=0,py=origin[1];y<grid->dims[1];y+=stride[1],py++){
	  int yi=(zi+y)*grid->dims[0];
	  int pyi=(pzi+py)*parent->dims[0];
	  for(x=0,px=origin[0];x<grid->dims[0];x+=stride[0],px++){
	    register int pindex=px+pyi;
	    register int cindex=x+yi;
	    // ugh! finally we get to assign parents to the nodes
	    (*pnodes)[pindex].child=&((*nodes)[cindex]);
	    (*nodes)[cindex].parent=&((*pnodes)[pindex]);
	  }
	}
      }
    }
  }
}

void AmrUcdGridHierarchy::print(){
  //************Print Some Results***********************
  for(int i=0;i<levelinfo.getSize();i++){
    AmrUcdGrid *grid;
    printf("Level %u\n",i);
    FOREACH(grid,grids[i]){ // level 2 grids should not be here!!!
      printf("-------------Slice grid %lu ID=%u---------\n",(unsigned long)grid,grid->gridID);
      int x,y,index;
      
      printf("=====: ");
      for(x=0;x<grid->dims[0];x++){
	// print a row
	printf("%5u ",x);
      }
      
      //puts("-----------------------------")
      for(y=0,index=0;y<grid->dims[1];y++,index+=grid->dims[0]){
	printf("\n[%3u]X: ",y);
	for(x=0;x<grid->dims[0];x++){
	  // print a row
	  printf("%5.3f ",grid->points[index+x].location.cartesian.x);
	}
	printf("\n     Y: ");
	for(x=0;x<grid->dims[0];x++){
	  // print a row
	  printf("%5.3f ",grid->points[index+x].location.cartesian.y);
	}	
	printf("\n     Z: ");
	for(x=0;x<grid->dims[0];x++){
	  // print a row
	  printf("%5.3f ",grid->points[index+x].location.cartesian.z);
	}
	printf("\n     P: ");
	for(x=0;x<grid->dims[0];x++){
	  if(grid->points[index+x].parent)
	    printf("%5u ",grid->points[index+x].parent->gridID);
	  else
	    printf("None  ");
	}
	printf("\n     C: ");
	for(x=0;x<grid->dims[0];x++){
	  if(grid->points[index+x].child)
	    printf("%5u ",grid->points[index+x].child->gridID);
	  else
	    printf(" None ");
	}
      }
      printf("\n");
    }
  }
}

void AmrUcdGridHierarchy::buildUCD(FlexArray<AmrNode*> &nodes,FlexArray<int> &cells){
  // OK, first build pointlist by adding only cells without children
  // for each one added to list, set its index member
  // also set index for all parents recursively. (have recursive setindex in 
  // as a static member function of each node)
  int ctr,i;
  //**** Init all node indexes to -1
  printf("AmrUcdGridHierarchy::buildUCD(): Initialize UCD\n");
  for(i=0;i<levelinfo.getSize();i++){
    AmrUcdGrid *grid;
    FOREACH(grid,grids[i]){
      AmrNode *points=grid->points.getData();
      for(long idx=0;idx<grid->nelements;idx++)
	points[i].index=-1;
    }
  }
  // print();
  //**** Top-down assignment of indices for each node without child
  // If a node has a child, then the child will be superceding it 
  // for indices
  printf("\ttop-down assign indices\n");
  for(i=levelinfo.getSize()-1;i>=0;--i){
    printf("\tLevel[%u]\n",i);
    AmrUcdGrid *grid;
    FOREACH(grid,grids[i]){
      AmrNode *points=grid->points.getData();
      for(long idx=0;idx<grid->nelements;idx++){      
	if(!points[idx].child){ // if no child
	  // then assign unique index
	  points[idx].setIndex(nodes.getSize());
	  nodes.append(points+idx); // append to nodelist
	}
      }
    }
  }
  printf("Number of nodes in all = %u\n",nodes.getSize());
  // OK, nows the time in schprockets when we dance...
  //**** Top-down building of Hexahedral cells
  //printf("now build hexahedral cells\n");
  for(i=levelinfo.getSize()-1;i>=0;i--){
    printf("\tL%u\n",i);
    AmrUcdGrid *grid;
    int nodeindices[8],mapped[8];
    FOREACH(grid,grids[i]){ 
      AmrNode *points=grid->points.getData();
      // verify that at least one point doesn't have a child
      // if all have children then continue
      int x,y,z;
      int *dims=grid->dims;
      int zstep=dims[1]*dims[0];
      int ystep=dims[0]; // xstep=1
      int idx;
      // setup the cell indices
      for(idx=0;idx<8;idx++){
		nodeindices[idx]=0;
      }
      for(idx=0;idx<4;idx++) {
		nodeindices[idx]+=zstep;
	  }
      for(idx=0;idx<=4;idx+=4){
		nodeindices[idx+2]+=1;
		nodeindices[idx+3]+=1;
		nodeindices[idx]+=ystep;
		nodeindices[idx+3]+=ystep;
      }
      //***** Build the Cell list 
      for(z=0;z<(dims[2]-1);z++){
	//printf("\t\tZ%u\n",z);
	for(y=0;y<(dims[1]-1);y++){
	  for(x=0;x<(dims[0]-1);x++){
	    // OK, now we find the real indices and 
	    for(idx=0;idx<8;idx++){ 
	      mapped[idx]=points[nodeindices[idx]].index;
	      if(mapped[idx]>nodes.getSize()){
			printf("indexing error at grid %u, position %u,%u,%u idx[%u]=%u\n",
		       grid->gridID,x,y,z,idx,mapped[idx]);
	      }
	    }
	    // append to cells list
	    cells.append(mapped,8); // copies the mapped array?
	    // now we save this sequence;
	    for(idx=0;idx<8;idx++) (nodeindices[idx])++;
	  }
	  for(idx=0;idx<8;idx++) (nodeindices[idx])++;
	}
	for(idx=0;idx<8;idx++) (nodeindices[idx])+=ystep;
      }
    }
  }
}
