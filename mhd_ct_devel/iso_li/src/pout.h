void Pout( char * string, int i1 = -12345, int i2 = -12345,
           int i3 = -12345, int i4 = -12345, int i5 = -12345);

void PoutF( char * string, float f1 = -12345.6,float f2 = -12345.6,float f3 = -12345.6,
	    float f4 = -12345.6,float f5 = -12345.6,float f6 = -12345.6 );

void dump(float *A, int nx, int ny, int nz, int nb, char * filename);
int WriteInThisF(int flag);
//<dcc>
void WriteCube(float * array, int Dims[], char* string, int dNum, int gNum);
//</dcc>
double wall_time (char * string);
