struct PhotonBuffer {
  float		Photons;
  int		Type;
  float		Energy;                   
  float		ColumnDensity;
  double	CrossSection;
  FLOAT		EmissionTimeInterval;
  FLOAT		EmissionTime;
  FLOAT		CurrentTime;
  FLOAT		Radius;
  long		ipix;
  int		level;
  FLOAT		SourcePosition[3];
  float		SourcePositionDiff;
  int		SuperSourceID;
};

struct GroupPhotonList {
  int ToLevel;
  int ToGrid;
  //  int FromLevel;
  //  int FromGrid;
  PhotonBuffer buffer;
};

