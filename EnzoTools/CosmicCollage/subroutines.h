
/* in file Utilities.C */
int chooseaxis(int proj_number);
void chooseoffset(double *xoffset, double *yoffset, int proj_number);

/* in file ProjectionRoutines.C */
void InitializeProjectionArray(void);
void EraseProjectionArray(void);
void ReadInProjectionArray(int axis, int boxnumber);
void RebinProjectionArray(int thisprojection, double xoffset, double yoffset);
void ApplyCorrectionsToRebinnedArray(int thisprojection);
void AddRebinnedProjectionToMaster(void);
void DeleteLoopArrays(void);
void CreateAndZeroLoopArrays(int thisprojection);
void WriteOutProjectionArray(void);
void ReadInUserData(void);
void WriteOutRebinnedProjectionArray(int thisprojection);

void DeleteClusterArrays(void);
void ReadClusterFile(int axis, int boxnumber);
void ModifyArraysWithClusterInfo(int axis, int boxnumber);

/* in file UserDefines.C */
void UserDefines(void);
