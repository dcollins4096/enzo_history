#ifndef PROJECTION_H
#define PROJECTION_H

class Projection : public AnalysisBaseClass {

 public:
  Projection( TopGridData *mdata, HierarchyEntry *topgrid, 
	      int num_fields, int *fields, int axis = 0, 
	      int max_level = 0, char * ProjectionTitle = NULL,
	      FLOAT *left = NULL, FLOAT *right = NULL);
  ~Projection();

  void Project();

  // create a new file with name = outputname
  void Write(char *outputname = NULL);
  void Plot(char *outputname = NULL, int color_table = 0);

  int GetDimension(int dim){return ProjectionDims[dim];};

  EnzoArrayFloat *CreateFieldArrayFloat(field_type_int field);
  EnzoArrayFloat *CreateFieldArrayFloat(char *field_name);

 private:
  char * ProjectionTitle;
  int ProjectionAxis;
  int ProjectionDims[2];
  FLOAT ProjectionDX[2];

  // these the mapping of x, y, z in to the height and 
  // width of the projection
  // Dim1 is width, Dim2 is height
  // Dim0 is the projection_axis defined in the
  // parameters
  int Dim1;
  int Dim2;

  // ProjectionDims[0]*ProjectionDims[1]
  // for the lazy compiler
  int BufferSize;
 
  void LoopProject( HierarchyEntry *Grid,
		    int current_level);

  ealInt *LevelField;
  ealFloat **ProjectionField;
  int *FieldType;

  int NumberOfFields;

  void ProjectGrid(grid *Grid, int CurrentLevel);
  void ProjectGridParticles(grid *Grid);
    
};

#endif /* PROJECTION_H */
