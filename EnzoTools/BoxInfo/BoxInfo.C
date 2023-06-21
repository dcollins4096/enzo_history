/* ------------------------------- BoxInfo.C ------------------------- *

   This tool goes through an Enzo AMR hierarchy and calculates a huge
   number of things about the entire box.  A separate tool (based on
   this one) will calculate info for just the most massive halo.

   This tool no longer makes assumptions on what MultiSpecies is set
   to.  It currently supports 0,1, and 2.  Note that any distribution
   functions which are added must be done so in the proper place based on
   the value of MultiSpecies.

   Quantities calculated for a given data dump:
   - maximum baryon density in the box (and its position)
   - Temperature and H2 fraction at that point
   - clumping factor of entire volume
   - Mean temperature, entropy, H2 fraction, e- fraction, H- fraction 
     of entire volume
   - maximum, minimum  temperature, density, H2 fraction, 
     e- fraction, H- fraction
   - Metal Mass/Volume Filling Fractions
   - Mean Metallicity Mass and Volume weighted
   - 2D distribution functions (both mass and volume-weighted), including:
       - Temperature vs. overdensity
       - Entropy vs. overdensity
       - Temperature vs. entropy
       - Jeans mass vs. overdensity
       - H2 fraction (rho_H2/(0.76*rho_tot)) vs. overdensity
       - f_H- vs. overdensity
       - f_e- vs. overdensity
       - f_e- vs. f_H2
       - f_H- vs. f_H2
       - f_e- vs. f_H-
   - 1D distribution functions (both mass and volume-weighted), including:
       - Temperature
       - overdensity
       - entropy
       - H2 fraction
       - e- fraction
       - H- fraction

 * ------------------------------------------------------------------- */

#include "BoxInfo.h"
/*------------------------------------------------------------------
 *MAKE SURE TO ADD NEW 1D and 2D Distributions for multispecies = 1!
 *----------------------------------------------------------------*/



/*--------------------------- main -----------------------------
 *
 *
 *--------------------------------------------------------------*/
int main(int argc, char *argv[]){
  char *inputfilename=NULL;
  char *halofilename=NULL;
  int junk,i,return_value;
  inputfilename=argv[1];


  boxsize = omegamatter = omegalambda = hubble = -1.0;
  rootgridsize = multispecies = -1;

  // read parameter file, get some values
  return_value =  ReadParameterFile(inputfilename);

  if(DEBUG) fprintf(stderr,"%lf %lf %lf %lf %lf %d %d %lf\n", boxsize,omegamatter,
		    omegalambda,redshift,hubble,rootgridsize,multispecies, currenttime);

  // this code assumes that we have molecular hydrogen, etc.  Hence multispecies > 1
  if(multispecies >2){
    fprintf(stderr,"MultiSpecies = %d. Should be <= 2\n",multispecies);
    exit(-1);
  }
  fprintf(stderr, "Multispecies = %d. \n", multispecies);

  if(DEBUG) fprintf(stderr,"Main: about to get hierarchy info\n");

  // get hierarchy file name
  int inlen = (int) strlen(inputfilename);
  char *hierfile = new char[inlen+6];
  strcpy(hierfile,inputfilename);
  hierfile=strcat(hierfile,".hierarchy");

  if(DEBUG) fprintf(stderr,"Main: hierarchy file name is:  %s\n",hierfile);

  // get total number of grids - to generate grid information
  total_number_grids = NumberOfGrids(hierfile);

  if(DEBUG) fprintf(stderr,"Main: there are %i grids\n",total_number_grids);

  // declares grid arrays - basically just uses new and sets things to zero
  DeclareGridArrays();

  if(DEBUG) fprintf(stderr,"Main: entering GetGridInfo\n");

  // get grid information from hierarchy file
  return_value = GetGridInfo(total_number_grids,hierfile);
  assert( return_value == total_number_grids );

  delete [] hierfile;

  // set all of the distribution function arrays to zero (or whatever
  // they should be) and also set mean values, min, max values to
  // something appropriate.
  if(DEBUG) fprintf(stderr,"Main: entering SetAllArraysAndValues\n");
  SetAllArraysAndValues();

  // Loop through grids and get information for each grid, making sure
  // to zero under subgrids
  if(DEBUG) fprintf(stderr,"Main: entering GetCellInformation\n");
  for(i=0; i<total_number_grids; i++){
    return_value=GetCellInformationFromGrid(i,total_number_grids);
    assert( return_value != FAILURE );    
  }

  // This cleans up all of the mean values by dividing through by
  // total mass, etc. and then cleans up the 1D distribution functions.
  // It also calculates the cumulative distribution functions in 1D.
  if(DEBUG) fprintf(stderr,"Main: entering CleanUpArraysAndValues\n");
  CleanUpArraysAndValues();

  // Actually write everything out - this creates a huge number of files.
  if(DEBUG) fprintf(stderr,"Main: entering OutputAllInformation\n");
  OutputAllInformation(argv[1]);

  // Out of here.
  if(DEBUG) fprintf(stderr,"Main: exiting...\n");
  exit(0);
}


/*------------------ OutputAllInformation -----------------------
 *
 *  Pretty straightforward, really.  Writes all of the information
 *  out to files.
 *
 *--------------------------------------------------------------*/
void OutputAllInformation(char *infilename){
  if(DEBUG) fprintf(stderr,"In OutputAllInformation  \n");
  int i, numpixels;
  FILE *output;

  // misc. information, mean information
  output = fopen("misc_and_meaninfo.dat","w");
  fprintf(output,"# mean information file\n");
  fprintf(output,"InputFileName       = %s\n",infilename);
  fprintf(output,"Redshift            = %lf\n",redshift);
  fprintf(output,"CurrentTime         = %lf\n",currenttime);
  fprintf(output,"BoxSize             = %lf\n",boxsize);
  fprintf(output,"RootGridSize        = %d\n",rootgridsize);
  fprintf(output,"TotalNumberOfGrids  = %d\n",total_number_grids);
  fprintf(output,"TotalMassInSim      = %e\n",totalmass);
  fprintf(output,"TotalVolumeInSim    = %e\n",totalvolume);
  fprintf(output,"\n");
  fprintf(output,"MaxBaryonDensity      = %e\n",maxbaryondens);
  fprintf(output,"MaxBaryonTemperature  = %e\n",maxbar_temp);
  if(multispecies == 2){
    fprintf(output,"MaxBaryonH2fraction   = %e\n",maxbar_H2frac);
  }
  fprintf(output,"MaxBaryonPosition     = %.12lf %.12lf %.12lf\n",maxbar_xpos, maxbar_ypos, maxbar_zpos);
  fprintf(output,"\n");
  fprintf(output,"RhoMeanMassWeighted           = %e\n",rhobar_mean_mw);
  fprintf(output,"RhoMeanVolumeWeighted         = %e\n",rhobar_mean_vw);
  fprintf(output,"RhoSquaredMeanMassWeighted    = %e\n",rhobarsquared_mean_mw);
  fprintf(output,"RhoSquaredMeanVolumeWeighted  = %e\n",rhobarsquared_mean_vw);
  fprintf(output,"ClumpingFactorMassWeighted    = %e\n",clumping_factor_mw);
  fprintf(output,"ClumpingFactorVolumeWeighted  = %e\n",clumping_factor_vw);
  fprintf(output,"\n");
  fprintf(output,"MeanTemperatureMassWeighted          = %e\n",meantemp_mw);
  fprintf(output,"MeanTemperatureVolumeWeighted        = %e\n",meantemp_vw);
  fprintf(output,"MeanEntropyMassWeighted              = %e\n",meanent_mw);
  fprintf(output,"MeanEntropyVolumeWeighted            = %e\n",meanent_vw);
  if(multispecies == 2){ 
    fprintf(output,"MeanH2FractionMassWeighted           = %e\n",meanH2frac_mw);
    fprintf(output,"MeanH2FractionVolumeWeighted         = %e\n",meanH2frac_vw);
    fprintf(output,"MeanHMinusMassWeighted               = %e\n",meanHMfrac_mw);
    fprintf(output,"MeanHMinusVolumeWeighted             = %e\n",meanHMfrac_vw);
  }
  if(multispecies > 0){
    fprintf(output,"MeanElectronFractionMassWeighted     = %e\n",meanelecfrac_mw);
    fprintf(output,"MeanElectronFractionVolumeWeighted   = %e\n",meanelecfrac_vw);
    fprintf(output,"MeanHIIFractionMassWeighted          = %e\n",meanHIIfrac_mw);
    fprintf(output,"MeanHIIFractionVolumeWeighted        = %e\n",meanHIIfrac_vw);
  }
  fprintf(output,"\n");
  fprintf(output,"MaximumTemperatureInVolume        = %e\n",maxtemp);
  fprintf(output,"MinimumTemperatureInVolume        = %e\n",mintemp);
  fprintf(output,"MaximumOverDensityInVolume        = %e\n",maxdens);
  fprintf(output,"MinimumOverDensityInVolume        = %e\n",mindens);
  fprintf(output,"MaximumEntropyInVolume            = %e\n",maxent);
  fprintf(output,"MinimumEntropyInVolume            = %e\n",minent);
  if(multispecies == 2){
  fprintf(output,"MaximumH2FractionInVolume         = %e\n",maxH2frac);
  fprintf(output,"MinimumH2FractionInVolume         = %e\n",minH2frac);
  fprintf(output,"MaximumElectronFractionInVolume   = %e\n",maxelecfrac);
  fprintf(output,"MinimumElectronFractionInVolume   = %e\n",minelecfrac);
  fprintf(output,"MaximumHMinusFractionInVolume     = %e\n",maxHMfrac);
  fprintf(output,"MinimumHMinusFractionInVolume     = %e\n",minHMfrac);
  }
  if(multispecies > 0){
    fprintf(output,"MaximumElectronFractionInVolume   = %e\n",maxelecfrac);
    fprintf(output,"MinimumElectronFractionInVolume   = %e\n",minelecfrac);
    fprintf(output,"MaximumHIIFractionInVolume        = %e\n",maxHIIfrac);
    fprintf(output,"MinimumHIIFractionInVolume        = %e\n",minHIIfrac);
  }
  fprintf(output,"MaximumJeansMassInVolume          = %e\n",maxjeans);
  fprintf(output,"MinimumJeansMassInVolume          = %e\n",minjeans);

  fprintf(output,"\n");
  fprintf(output,"Metal mass filling fraction    = %e\n", metalmff/totalmass);
  fprintf(output,"Metal volume filling fraction  = %e\n", metalvff/totalvolume);

  if(metalmff > 0){
      fprintf(output,"Mean Metallicity of Gas MW     = %e\n", metalmass/metalmff/0.022);
      fprintf(output,"Mean Metallicity of Gas VW     = %e\n", metalvol/metalvff/0.022);
  }
  if(metalmff == 0){
      fprintf(output,"Mean Metallicity of Gas MW     = %e\n", 0.0);
      fprintf(output,"Mean Metallicity of Gas VW     = %e\n", 0.0);
  }

  fprintf(output,"\n");
  fprintf(output,"DensityLogMinAndMaxLimis        = %lf %lf\n",DMIN,DMAX);
  fprintf(output,"EntropyLogMinAndMaxLimits       = %lf %lf\n",SMIN,SMAX);
  fprintf(output,"TemperatureLogMinAndMaxLimits   = %lf %lf\n",TMIN,TMAX);
  if(multispecies == 2){
    fprintf(output,"H2FracLogMinAndMaxLimits        = %lf %lf\n",H2MIN,H2MAX);
    fprintf(output,"ElectronFracLogMinAndMaxLimits  = %lf %lf\n",EMIN,EMAX);
    fprintf(output,"HMinusFracLogMinAndMaxLimits    = %lf %lf\n",HMMIN,HMMAX);
  }
  if(multispecies > 0){
   fprintf(output,"ElectronFracLogMinAndMaxLimits  = %lf %lf\n",EMIN,EMAX);
   fprintf(output,"HIIFracLogMinAndMaxLimits       = %lf %lf\n",HIIMIN,HIIMAX);
  }
  fprintf(output,"JeansMassMinAndMaxLimits        = %lf %lf\n",JMIN,JMAX);
  fprintf(output,"UserSpecifiedNumberOfBins       = %d\n",NUMBINS);
  fprintf(output,"UserSpecifiedOmegaBaryon        = %lf\n",OMEGA_BARYON);
  fprintf(output,"UserSpecifiedOmegaMatter        = %lf\n",OMEGA_MATTER);
  fprintf(output,"UserSpecifiedDebugFlag          = %d\n",DEBUG);
  fprintf(output,"UserSpecifiedVerboseDebugFlag   = %d\n",VERBOSEDEBUG);
  fprintf(output,"\n");  // add a newline just in case


  fclose(output);

  // 1D distribution functions
  output = fopen("temperature_1ddist.dat","w");
  for(i=0; i < NUMBINS; i++)
    fprintf(output,"%lf %e %e %e %e\n", temp_1ddist_value[i],
	    temp_1ddist_mw[i], temp_1ddist_mw_cum[i],
	    temp_1ddist_vw[i], temp_1ddist_vw_cum[i]);
  fclose(output);

  output = fopen("overdensity_1ddist.dat","w");
  for(i=0; i < NUMBINS; i++)
    fprintf(output,"%lf %e %e %e %e\n", od_1ddist_value[i],
	    od_1ddist_mw[i], od_1ddist_mw_cum[i],
	    od_1ddist_vw[i], od_1ddist_vw_cum[i]);
  fclose(output);

  output = fopen("entropy_1ddist.dat","w");
  for(i=0; i < NUMBINS; i++)
    fprintf(output,"%lf %e %e %e %e\n", ent_1ddist_value[i],
	    ent_1ddist_mw[i], ent_1ddist_mw_cum[i],
	    ent_1ddist_vw[i], ent_1ddist_vw_cum[i]);
  fclose(output);


  output = fopen("jeansmass_1ddist.dat","w");
  for(i=0; i < NUMBINS; i++)
    fprintf(output,"%lf %e %e %e %e\n", jeans_1ddist_value[i],
	    jeans_1ddist_mw[i], jeans_1ddist_mw_cum[i],
	    jeans_1ddist_vw[i], jeans_1ddist_vw_cum[i]);
  fclose(output);

  if(multispecies == 2){
    output = fopen("H2fraction_1ddist.dat","w");
    for(i=0; i < NUMBINS; i++)
      fprintf(output,"%lf %e %e %e %e\n", H2frac_1ddist_value[i],
	      H2frac_1ddist_mw[i], H2frac_1ddist_mw_cum[i],
	      H2frac_1ddist_vw[i], H2frac_1ddist_vw_cum[i]);
    fclose(output);
    
    output = fopen("Hminusfraction_1ddist.dat","w");
    for(i=0; i < NUMBINS; i++)
      fprintf(output,"%lf %e %e %e %e\n", HMfrac_1ddist_value[i],
	      HMfrac_1ddist_mw[i], HMfrac_1ddist_mw_cum[i],
	      HMfrac_1ddist_vw[i], HMfrac_1ddist_vw_cum[i]);
    fclose(output);
  }
  if(multispecies > 0){
    output = fopen("electronfraction_1ddist.dat","w");
    for(i=0; i < NUMBINS; i++)
      fprintf(output,"%lf %e %e %e %e\n", efrac_1ddist_value[i],
	      efrac_1ddist_mw[i], efrac_1ddist_mw_cum[i],
	      efrac_1ddist_vw[i], efrac_1ddist_vw_cum[i]);
    fclose(output);
    /*--------------ADD HII distribution Functions when neeed----*/
  }
  // ------------ 2D distribution functions
  numpixels = NUMBINS;

  // temperature vs. overdensity
  output = fopen("tempvsrho_mw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(temp_od_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);

  output = fopen("tempvsrho_vw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(temp_od_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);

  // entropy vs. overdensity
  output = fopen("entvsrho_mw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(ent_od_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);

  output = fopen("entvsrho_vw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(ent_od_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);

  // temperature vs. entropy
  output = fopen("tempvsent_mw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(temp_ent_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);

  output = fopen("tempvsent_vw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(temp_ent_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);

  // jeans masss vs. overdensity
  output = fopen("jeansvsrho_mw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(jeans_od_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);

  output = fopen("jeansvsrho_vw.dat","w");
  fwrite(&numpixels,sizeof(int),1,output);
  fwrite(jeans_od_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
  fclose(output);
  
  if(multispecies == 2){
    // H2 fraction vs. overdensity
    output = fopen("H2fracvsrho_mw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(H2frac_od_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    output = fopen("H2fracvsrho_vw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(H2frac_od_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    // H minus fraction vs. overdensity
    output = fopen("Hmfracvsrho_mw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(HMfrac_od_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    output = fopen("Hmfracvsrho_vw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(HMfrac_od_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
   
    // electron fraction vs. H2fraction
    output = fopen("efracvsH2frac_mw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(efrac_H2frac_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    output = fopen("efracvsH2frac_vw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(efrac_H2frac_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    // H minus fraction vs. H2fraction
    output = fopen("HMfracvsH2frac_mw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(HMfrac_H2frac_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    output = fopen("HMfracvsH2frac_vw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(HMfrac_H2frac_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    // electron fraction vs. H minus fraction
    output = fopen("efracvsHMfrac_mw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(efrac_HMfrac_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    output = fopen("efracvsHMfrac_vw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(efrac_HMfrac_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    
  }
  if(multispecies > 0){
    // electron fraction vs. overdensity
    output = fopen("efracvsrho_mw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(efrac_od_2ddist_mw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    
    output = fopen("efracvsrho_vw.dat","w");
    fwrite(&numpixels,sizeof(int),1,output);
    fwrite(efrac_od_2ddist_vw,sizeof(double),NUMBINS*NUMBINS,output);
    fclose(output);
    /*------------------ADD HII 2d distributions------------------*/
  }
 
  if(DEBUG) fprintf(stderr,"Leaving OutputAllInformation  \n");
  return;
}


/*------------------ CleanUpArraysAndValues ----------------------
 *
 *  This function takes all of the mean weighted quantities and 
 *  divides through by the appropriate weight.  It also calculates
 *  cumulative distribution functions and normalizes all of the 1D
 *  distribution functions correctly.
 *
 *----------------------------------------------------------------*/
void CleanUpArraysAndValues(void){

  int i;
  double temptotal_mw, odtotal_mw, enttotal_mw, H2fractotal_mw, 
    efractotal_mw, HMfractotal_mw, jeanstotal_mw, HIIfractotal_mw,
    temptotal_vw, odtotal_vw, enttotal_vw, H2fractotal_vw, 
    efractotal_vw, HMfractotal_vw, jeanstotal_vw, HIIfractotal_vw;

  // zero out all of the total values - these are used in the 
  // cumulative dist'n function stuff 
  temptotal_mw = odtotal_mw = enttotal_mw = H2fractotal_mw = 
    efractotal_mw = HMfractotal_mw = jeanstotal_mw = HIIfractotal_mw =
    temptotal_vw = odtotal_vw = enttotal_vw = H2fractotal_vw = 
    efractotal_vw = HMfractotal_vw = jeanstotal_vw = HIIfractotal_vw = 0.0;

  // divide all mean quantities through by the appropriate weights
  rhobar_mean_mw /= totalmass;
  rhobar_mean_vw /= totalvolume;

  rhobarsquared_mean_mw /= totalmass;
  rhobarsquared_mean_vw /= totalvolume;

  // calculate clumping factor two different ways - really,
  // traditionally clumping factor is calculated in a volume-weighted
  // way, so that's the one I should use.  Still, it's interesting to see.
  clumping_factor_mw = rhobarsquared_mean_mw / rhobar_mean_mw / rhobar_mean_mw;
  clumping_factor_vw = rhobarsquared_mean_vw / rhobar_mean_vw / rhobar_mean_vw;

  meantemp_mw /= totalmass;
  meantemp_vw /= totalvolume;

  meanent_mw /= totalmass;
  meanent_vw /= totalvolume;

  if(multispecies == 2){
    meanH2frac_mw /= totalmass;
    meanH2frac_vw /= totalvolume;
    
    meanHMfrac_mw /= totalmass;
    meanHMfrac_vw /= totalvolume;
  }
  if(multispecies > 0){
    meanelecfrac_mw /= totalmass;
    meanelecfrac_vw /= totalvolume;
    
    meanHIIfrac_mw /= totalmass;
    meanHIIfrac_vw /= totalvolume;
  }
  
  // calculate cumulative distribution function, counting backwards 
  // from the highest bin.  This is done by keeping track of the total
  // value starting at the highest bin and summing backwards.  We'll then
  // normalize to one later.
  for(i=NUMBINS-1; i>= 0; i--){
    
    // calculate total cumulative (working backwards)
    temptotal_mw += temp_1ddist_mw[i];
    temptotal_vw += temp_1ddist_vw[i];

    odtotal_mw += od_1ddist_mw[i];
    odtotal_vw += od_1ddist_vw[i];

    enttotal_mw += ent_1ddist_mw[i];
    enttotal_vw += ent_1ddist_vw[i];

    if(multispecies == 2){
      H2fractotal_mw += H2frac_1ddist_mw[i];
      H2fractotal_vw += H2frac_1ddist_vw[i];
      
      efractotal_mw += efrac_1ddist_mw[i];
      efractotal_vw += efrac_1ddist_vw[i];
      
      HMfractotal_mw += HMfrac_1ddist_mw[i];
      HMfractotal_vw += HMfrac_1ddist_vw[i];
    }
    if(multispecies > 0){
      efractotal_mw += efrac_1ddist_mw[i];
      efractotal_vw += efrac_1ddist_vw[i];
    
      HIIfractotal_mw += HIIfrac_1ddist_mw[i];
      HIIfractotal_vw += HIIfrac_1ddist_vw[i];
    }
    
    jeanstotal_mw += jeans_1ddist_mw[i];
    jeanstotal_vw += jeans_1ddist_vw[i];

    // actually set this value in the arrays that hold
    // the cumulative distribution functions
    temp_1ddist_mw_cum[i] = temptotal_mw;
    temp_1ddist_vw_cum[i] = temptotal_vw;

    od_1ddist_mw_cum[i] = odtotal_mw;
    od_1ddist_vw_cum[i] = odtotal_vw;

    ent_1ddist_mw_cum[i] = enttotal_mw;
    ent_1ddist_vw_cum[i] = enttotal_vw;
    
    if(multispecies == 2){
    H2frac_1ddist_mw_cum[i] = H2fractotal_mw;
    H2frac_1ddist_vw_cum[i] = H2fractotal_vw;

    HMfrac_1ddist_mw_cum[i] = HMfractotal_mw;
    HMfrac_1ddist_vw_cum[i] = HMfractotal_vw;
    }
    if(multispecies > 0){
      efrac_1ddist_mw_cum[i] = efractotal_mw;
      efrac_1ddist_vw_cum[i] = efractotal_vw;
      
      HIIfrac_1ddist_mw_cum[i] = HIIfractotal_mw;
      HIIfrac_1ddist_vw_cum[i] = HIIfractotal_vw; 
    }
    
    jeans_1ddist_mw_cum[i] = jeanstotal_mw;
    jeans_1ddist_vw_cum[i] = jeanstotal_vw;

  }


  // normalize all distribution functions to binsize and total mass, etc. 
  // do both cumulative and differential
  for(i=0; i<NUMBINS; i++){
    // differential: divide through by total mass
    // and binsize to make this comparable between data files with
    // different values of NUMBINS, etc.
    temp_1ddist_mw[i] /= (totalmass*tbin);
    temp_1ddist_vw[i] /= (totalvolume*tbin);

    od_1ddist_mw[i] /= (totalmass*dbin);
    od_1ddist_vw[i] /= (totalvolume*dbin);

    ent_1ddist_mw[i] /= (totalmass*sbin);
    ent_1ddist_vw[i] /= (totalvolume*sbin);
    
    if(multispecies == 2){
    H2frac_1ddist_mw[i] /= (totalmass*H2bin);
    H2frac_1ddist_vw[i] /= (totalvolume*H2bin);

    HMfrac_1ddist_mw[i] /= (totalmass*HMbin);
    HMfrac_1ddist_vw[i] /= (totalvolume*HMbin);
    }
    if(multispecies > 0){
      HIIfrac_1ddist_mw[i] /= (totalmass*HIIbin);
      HIIfrac_1ddist_vw[i] /= (totalvolume*HIIbin);
      
      efrac_1ddist_mw[i] /= (totalmass*ebin);
      efrac_1ddist_vw[i] /= (totalvolume*ebin);
    }
    
    jeans_1ddist_mw[i] /= (totalmass*jbin);
    jeans_1ddist_vw[i] /= (totalvolume*jbin);

    // cumulative: just divide by the total value (which
    // corresponds to the value in bin zero)
    temp_1ddist_mw_cum[i] /= temptotal_mw;
    temp_1ddist_vw_cum[i] /= temptotal_vw;

    od_1ddist_mw_cum[i] /= odtotal_mw;
    od_1ddist_vw_cum[i] /= odtotal_vw;

    ent_1ddist_mw_cum[i] /= enttotal_mw;
    ent_1ddist_vw_cum[i] /= enttotal_vw;

    if(multispecies == 2){
    H2frac_1ddist_mw_cum[i] /= H2fractotal_mw;
    H2frac_1ddist_vw_cum[i] /= H2fractotal_vw;

    efrac_1ddist_mw_cum[i] /= efractotal_mw;
    efrac_1ddist_vw_cum[i] /= efractotal_vw;

    HMfrac_1ddist_mw_cum[i] /= HMfractotal_mw;
    HMfrac_1ddist_vw_cum[i] /= HMfractotal_vw;
    }
    if(multispecies > 0){
      HIIfrac_1ddist_mw_cum[i] /= HIIfractotal_mw;
      HIIfrac_1ddist_vw_cum[i] /= HIIfractotal_vw;
      
      efrac_1ddist_mw_cum[i] /= efractotal_mw;
      efrac_1ddist_vw_cum[i] /= efractotal_vw;
    }
      
    jeans_1ddist_mw_cum[i] /= jeanstotal_mw;
    jeans_1ddist_vw_cum[i] /= jeanstotal_vw;
  }

  return;
}


/*-------------------- SetAllArraysAndValues -------------------
 *
 *  This routine zeros out various arrays, and also sets the correct
 *  values for the 1D distribution bings.  It sets mean values to 
 *  zero (which are summed over later) and sets max and min values
 *  for the box to totally ridiculous values so we'll write over it.
 *
 *--------------------------------------------------------------*/
void SetAllArraysAndValues(void){
  if(DEBUG) fprintf(stderr,"In SetAllArraysAndValues\n");
  int i,j;
  double enzo_densconversion;
  
  maxbaryondens = maxbar_temp  
    = maxbar_xpos = maxbar_ypos = maxbar_zpos 
    = -1.0;
  
  if(multispecies == 2)
    maxbar_H2frac = -1.0;
  
  
  rhobar_mean_mw = rhobarsquared_mean_mw 
    = rhobar_mean_vw = rhobarsquared_mean_vw 
    = clumping_factor_mw = clumping_factor_vw = meantemp_mw = meantemp_vw
    = meanent_mw = meanent_vw = 0.0;
  if(multispecies == 2){
    meanH2frac_mw = meanH2frac_vw
      = meanelecfrac_mw = meanelecfrac_vw
      = meanHMfrac_mw = meanHMfrac_vw 
      = meanHIIfrac_mw = meanHIIfrac_vw = 0.0;
  }
  if(multispecies >0){
    meanHIIfrac_mw = meanHIIfrac_vw
      = meanelecfrac_mw = meanelecfrac_vw = 0.0;
  }

  metalmff = metalvff = metalmass = metalvol = 0.0;

  maxtemp = maxdens = maxent = maxH2frac = maxelecfrac = maxHMfrac = maxHIIfrac = maxjeans = -1.0;
  
  mintemp = mindens = minent = minH2frac = minelecfrac = minHMfrac = minHIIfrac = minjeans = HUGE_NUMBER;
  
  totalmass = totalvolume = 0.0;
  
  // calculate bin sizes (in log10 space)
  dbin = (DMAX-DMIN)/NUMBINS;
  sbin = (SMAX-SMIN)/NUMBINS;
  tbin = (TMAX-TMIN)/NUMBINS;
  jbin = (JMAX-JMIN)/NUMBINS;
  if(multispecies == 2){
    H2bin = (H2MAX-H2MIN)/NUMBINS;
    HMbin = (HMMAX-HMMIN)/NUMBINS;
    HIIbin = (HIIMAX-HIIMIN)/NUMBINS;
  }
  if(multispecies > 0){
    H2bin = (H2MAX-H2MIN)/NUMBINS;
    ebin = (EMAX-EMIN)/NUMBINS;
  }

  // zero out 2D distribution functions
  for(i=0; i < NUMBINS; i++)
    for(j=0; j < NUMBINS; j++){
      temp_od_2ddist_mw[i][j] = temp_od_2ddist_vw[i][j]
	= ent_od_2ddist_mw[i][j] = ent_od_2ddist_vw[i][j]
	= temp_ent_2ddist_mw[i][j] = temp_ent_2ddist_vw[i][j]
	= jeans_od_2ddist_mw[i][j] = jeans_od_2ddist_vw[i][j]	
	= 0.0;
      if(multispecies == 2)
	H2frac_od_2ddist_mw[i][j] = H2frac_od_2ddist_vw[i][j]
	  = HMfrac_od_2ddist_mw[i][j] = HMfrac_od_2ddist_vw[i][j]
	  = efrac_od_2ddist_mw[i][j] = efrac_od_2ddist_vw[i][j]
	  = efrac_H2frac_2ddist_mw[i][j] = efrac_H2frac_2ddist_vw[i][j]
	  = HMfrac_H2frac_2ddist_mw[i][j] = HMfrac_H2frac_2ddist_vw[i][j]
	  = efrac_HMfrac_2ddist_mw[i][j] = efrac_HMfrac_2ddist_vw[i][j]
	  = 0.0;
      /* ADD CASE FOR MULTISPECIES = 1*/
    }
  /* Need to zero out any other 2D distribution functions added!*/
    
  // zero out 1D distribution functions
  for(i=0; i<NUMBINS; i++){
    temp_1ddist_mw[i] = temp_1ddist_vw[i]
      = od_1ddist_mw[i] = od_1ddist_vw[i]
      = ent_1ddist_mw[i] = ent_1ddist_vw[i]
      = jeans_1ddist_mw[i] = jeans_1ddist_vw[i]
      = 0.0;
    if(multispecies == 2)
      H2frac_1ddist_mw[i] = H2frac_1ddist_vw[i]
	= HMfrac_1ddist_mw[i] = HMfrac_1ddist_vw[i]
	= 0.0;
    if(multispecies > 0)
      efrac_1ddist_mw[i] = efrac_1ddist_vw[i]
	= HIIfrac_1ddist_mw[i] = HIIfrac_1ddist_vw[i]
	= 0.0;
  }
  /* Need to zero out any other 1D distribution functions added!*/


  // set bin values for 1D dist'n fctns - calculate vals
  // in center of bin so that the plots look nice 'n stuff
  for(i=0; i<NUMBINS; i++){

    if(i==0){  // first bin

      temp_1ddist_value[i] = tbin/2.0 + TMIN;
      od_1ddist_value[i] = dbin/2.0 + DMIN;
      ent_1ddist_value[i] = sbin/2.0 + SMIN;
      jeans_1ddist_value[i] = jbin/2.0 + JMIN;
      if(multispecies == 2){
	H2frac_1ddist_value[i] = H2bin/2.0 + H2MIN;	
	HMfrac_1ddist_value[i] = HMbin/2.0 + HMMIN;
      }
      if(multispecies > 0){
	efrac_1ddist_value[i] = ebin/2.0 + EMIN;
	HIIfrac_1ddist_value[i] = HIIbin/2.0 +HIIMIN;
      }
      
    } else {  // all other bins

      temp_1ddist_value[i] = (((double) i) + 0.5)*tbin + TMIN;
      od_1ddist_value[i] = (((double) i) + 0.5)*dbin + DMIN ;
      ent_1ddist_value[i] = (((double) i) + 0.5)*sbin + SMIN ;
      jeans_1ddist_value[i] = (((double) i) + 0.5)*jbin + JMIN ;
      if(multispecies == 2){
	H2frac_1ddist_value[i] = (((double) i) + 0.5)*H2bin + H2MIN ;
	HMfrac_1ddist_value[i] = (((double) i) + 0.5)*HMbin + HMMIN ;
      }
      if(multispecies > 0){
	efrac_1ddist_value[i] = (((double) i) + 0.5)*ebin + EMIN ;
	HIIfrac_1ddist_value[i] = (((double) i) + 0.5)*HIIbin + HIIMIN ;
      }
    }

  } // for(i=0; i<NUMBINS; i++){


  // density constant for cell mass calculations
  densconstant = omegamatter * 2.78 * 100000000000.0 * hubble * hubble; // msolar/mpc^3
  if(DEBUG) fprintf(stderr,"densconstant:  %e\n",densconstant);


  // need this for jeans mass prefactor
  enzo_densconversion = omegamatter * 1.8788e-29 * hubble* hubble* pow( (1.0+redshift), 3.0);

  /* jeans mass prefactor:  this is set up so that all the user has to
     do is multiply it by T^3/2 / rho^0.5 (where rho is in Enzo code units)
     And then you get the jeans mass in solar masses.

     Mjeans = (3/(4pi))^(1/2) * (5*kb/(G*mu*m_h))^(3/2) * T^(3/2) / rho^(1/2)

     kb = boltz. constant
     mu = mean molecular weight (assuming 1.22)
     m_h = hydrogen atom mass      */

  jeansmass_prefactor = sqrt( 3.0/4.0/3.14159/enzo_densconversion )
    * pow( (5.0*1.38e-16/6.67e-8/1.22/1.67e-24), 1.5  )
    / 1.989e33;

  if(DEBUG) fprintf(stderr,"Leaving SetAllArraysAndValues\n");
  return;
}


/*--------------------- DeclareGridArrays ----------------------
 *
 *  This simple function just uses new to declare the arrays for
 *  grid information like position, level, etc.  Also set various
 *  values to absurd values so we can error-check later.
 *
 *--------------------------------------------------------------*/
void DeclareGridArrays(void){
  if(DEBUG) fprintf(stderr,"In DeclareGridArrays\n");

  int i;

  /* initialize arrays which grid information is stored in
     (information such as level, grid bounds and number of
     cells, etc.) */
  gridlevel = new int[total_number_grids];  // level info

  // grid number of cells info
  griddx = new int[total_number_grids];  
  griddy = new int[total_number_grids];
  griddz = new int[total_number_grids];

  // grid bounds info
  gridlex = new double[total_number_grids];
  gridley = new double[total_number_grids];
  gridlez = new double[total_number_grids];
  gridrex = new double[total_number_grids];
  gridrey = new double[total_number_grids];
  gridrez = new double[total_number_grids];

  // grid file name info
  gridfilenames = new char*[total_number_grids];

  // number of particles info
  gridnump = new int[total_number_grids];

  /* initialize all of the arrays to negative values 
     (in the case of the int and double arrays) or
     to a string in the case of the filename array. */
  for(i=0; i<total_number_grids; i++){
    gridlevel[i]=griddx[i]=griddy[i]=griddz[i]=gridnump[i]=-1;

    gridlex[i]=gridley[i]=gridlez[i]=
      gridrex[i]=gridrey[i]=gridrez[i]=-1.0;

    gridfilenames[i]= new char[MAX_LINE_LENGTH];
  }

  if(DEBUG) fprintf(stderr,"Leaving DeclareGridArrays  \n");
  return;
}

/*---------------------- CleanGridArrays -----------------------
 *
 *  Cleans up the memory consumed by the grid arrays created in 
 *  DeclareGridArrays.  Nothing exciting here.
 *  
 *--------------------------------------------------------------*/
void CleanGridArrays(void){
  if(DEBUG) fprintf(stderr,"In CleanGridArrays  \n");

  int i;

  // clean up dynamically allocated stuff
  for(i=0; i<total_number_grids; i++)
    delete [] gridfilenames[i];
  
  delete [] gridfilenames;
  delete [] gridlevel;
  delete [] griddx;
  delete [] griddy;
  delete [] griddz;
  delete [] gridlex;
  delete [] gridley;
  delete [] gridlez;
  delete [] gridrex;  
  delete [] gridrey;  
  delete [] gridrez;  
  delete [] gridnump;

  if(DEBUG) fprintf(stderr,"Leaving CleanGridArrays\n");
  return;
}


/*------------------ GetCellInformationFromGrid ---------------------------
 *
 *  This function opens up the grid, extracts all of the baryon quantities
 *  we care about.  It then calls FlagGridCells and GetDensityInfo, which is
 *  the routine that actually calculates all of the quantities that we care
 *  about.
 *
 *-------------------------------------------------------------------------*/
int GetCellInformationFromGrid(int gridnum,int total_number_grids){

  if(DEBUG) fprintf(stderr,"in GetCellInformationFromGrid %d\n",gridnum);


  // we only need one of all of the following
  hid_t file_id;
  hid_t       mem_type_id;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       dsp_id;
  hid_t       typ_id;

  hsize_t     size;
  hsize_t     dims[4];
  hsize_t     xdims[4];
  hsize_t     maxdims[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  // need one of these per dataset!
  hid_t dens_dset_id,temp_dset_id, H2I_dset_id, 
    HM_dset_id, electron_dset_id,HII_dset_id, metal_dset_id;

  int i, ndims;

  // open grid file and extract various quantities of interest

  // open file - only once
  fprintf(stderr,"GetCellInformationFromGrid %s\n",gridfilenames[gridnum]);
  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  assert( file_id != h5_error );
  
  // open density dataset
  dens_dset_id = H5Dopen(file_id, "Density");
  assert( dens_dset_id != h5_error );

  // open temperature dataset
  temp_dset_id = H5Dopen(file_id,"Temperature");
  assert( temp_dset_id != h5_error );

  if(multispecies == 2){ 
      // open H2 density dataset
      H2I_dset_id = H5Dopen(file_id,"H2I_Density");
      assert( H2I_dset_id != h5_error );

      // open H minus density dataset
      HM_dset_id = H5Dopen(file_id,"HM_Density");
      assert( HM_dset_id != h5_error );
  }
  if(multispecies > 0){
      // open electron density dataset
      electron_dset_id = H5Dopen(file_id,"Electron_Density");
      assert( electron_dset_id != h5_error );

      // open H+ density dataset
      HII_dset_id = H5Dopen(file_id,"HII_Density");
      assert( HII_dset_id != h5_error );
  }

#ifdef USE_METALS
  // open metal density dataset
  metal_dset_id = H5Dopen(file_id, "Metal_Density");
  assert( metal_dset_id != h5_error );
#endif

  // open density dataspace (to get dimensions) 
  // only once!
  dsp_id = H5Dget_space(dens_dset_id);
  assert( dsp_id != h5_error );

  // get data type (only once!)
  typ_id = H5Dget_type(dens_dset_id);
  assert( typ_id != h5_error );
  

  // get dimensional information from dataspace (only once)
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  // from the dimensional information, calculate the size of the buffer.
  // only once!
  size = 1;
  if(DEBUG) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(DEBUG) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }
  if(DEBUG) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

  // get arrays
  if ( H5Tequal( typ_id, H5T_IEEE_F32BE ) )
    {
      // allocate buffers - one per projection buffer!
      densitybuff = new float[(int) size];
      temperaturebuff = new float[(int) size];
      H2Ibuff = new float[(int) size];
      Hminusbuff = new float[(int) size];
      electronbuff = new float[(int) size];
      HIIbuff = new float[(int) size];

#ifdef USE_METALS
      metalbuff = new float[(int) size];
#endif

      mem_type_id = H5T_NATIVE_FLOAT;

      // read density field into an array
      h5_status = H5Dread(dens_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, densitybuff);
      if(DEBUG) fprintf(stderr,"double read status %d for density field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      // read temperature field into an array
      h5_status = H5Dread(temp_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, temperaturebuff);
      if(DEBUG) fprintf(stderr,"double read status %d for temperature field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      if(multispecies == 2){
	  // read H2 density field into an array
	  h5_status = H5Dread(H2I_dset_id, mem_type_id, 
			      mem_dsp_id, file_dsp_id, 
			      H5P_DEFAULT, H2Ibuff);
	  if(DEBUG) fprintf(stderr,"double read status %d for H2 density field\n", 
			    (int) h5_status);
	  assert( h5_status != h5_error ); 


	  // read H minus density field into an array
	  h5_status = H5Dread(HM_dset_id, mem_type_id, 
			      mem_dsp_id, file_dsp_id, 
			      H5P_DEFAULT, Hminusbuff);
	  if(DEBUG) fprintf(stderr,"double read status %d for H minus density field\n", 
			    (int) h5_status);
	  assert( h5_status != h5_error ); 
      }
      
      if(multispecies > 0){
	  // read electron density field into an array
	  h5_status = H5Dread(electron_dset_id, mem_type_id, 
			      mem_dsp_id, file_dsp_id, 
			      H5P_DEFAULT, electronbuff);
	  if(DEBUG) fprintf(stderr,"double read status %d for electron density field\n", 
			    (int) h5_status);
	  assert( h5_status != h5_error );
	  
	  // read H+ density field into an array
	  h5_status = H5Dread(HII_dset_id, mem_type_id, 
			      mem_dsp_id, file_dsp_id, 
			      H5P_DEFAULT, HIIbuff);
	  if(DEBUG) fprintf(stderr,"double read status %d for H+ density field\n", 
			    (int) h5_status);
	  assert( h5_status != h5_error ); 
      }  

#ifdef USE_METALS
      // read metal density field into an array
      h5_status = H5Dread(metal_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT, metalbuff);
      if(DEBUG) fprintf(stderr,"double read status %d for metal density field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 
#endif


    }

  // ---------- close hdf5 file, doing appropriate error checking
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );


  // ---------- must close each dataset - one per buffer extracted!
  h5_status = H5Dclose(dens_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(temp_dset_id);
  assert( h5_status != h5_error );

  if(multispecies == 2){
      h5_status = H5Dclose(H2I_dset_id);
      assert( h5_status != h5_error );
      
      h5_status = H5Dclose(HM_dset_id);
      assert( h5_status != h5_error );
  }
  if(multispecies > 0){
      h5_status = H5Dclose(electron_dset_id);
      assert( h5_status != h5_error );
      
      h5_status = H5Dclose(HII_dset_id);
      assert( h5_status != h5_error );
  }

#ifdef USE_METALS
  h5_status = H5Dclose(metal_dset_id);
  assert( h5_status != h5_error );
#endif

  // ---------- close file
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  // ---------- create flag field and zero it out
  flagbuffer = new int[(int) size];
  for(i=0; i < ((int) size) ; i++) flagbuffer[i] = 0;

  // ---------- check density buffer as a debug operation
  for(i=0; i < ((int) size) ; i++){
    if(densitybuff[i]<=0.0) 
      fprintf(stderr,"GetCellInformationFromGrid: negative or zero density value!\n");
  }

  // flag grid cells here!
  FlagGridCells(gridnum,total_number_grids);

  // CALL THE ROUTINE THAT ACTUALLY GETS THE INFORMATION! 
  GetDensityInfo(gridnum);


  // ------------------- erase buffers - once per buffer!
  delete [] flagbuffer;
  delete [] densitybuff;
  delete [] temperaturebuff;
  delete [] H2Ibuff;
  delete [] Hminusbuff;
  delete [] electronbuff;
  delete [] HIIbuff;

#ifdef USE_METALS
  delete [] metalbuff;
#endif

  if(DEBUG) fprintf(stderr,"exiting GetCellInformationFromGrid\n");
  return SUCCESS;
}



/*------------------------- GetDensityInfo ----------------------------
 *
 *  This is the routine that actually calculates all of the quantites
 *  we care about.  Loop through grid cells and for all cells that are
 *  the most highly-resolved, we calculate all of the quantites shown
 *  at the top of the file.
 * 
 *---------------------------------------------------------------------*/
int GetDensityInfo(int gridnum){
  if(DEBUG) fprintf(stderr,"in GetDensityInfo\n");

  double cellsize,cellvol,cellmass, ccx,ccy,ccz, H2frac, efrac, HMfrac,
    overdensity,entropy,jeansmass,HIIfrac;

  int i,j,k,cellindex,totalcells, logodbin, logtempbin, logjeansbin, 
    logefracbin, logH2fracbin, logHminusfracbin,logentbin, logHIIfracbin;


  totalcells = griddx[gridnum]*griddy[gridnum]*griddz[gridnum];

  // size of grid cell, in code units
  cellsize = (gridrex[gridnum] - gridlex[gridnum]) / ((double) griddx[gridnum] );

  // cell volume in code units
  cellvol = cellsize*cellsize*cellsize;

  // density -> mass conversion, cell volume included
  massconv = densconstant * cellvol * boxsize * boxsize * boxsize 
    / hubble / hubble / hubble;
  
  // loop over all cells
  for(k=0; k < griddz[gridnum]; k++)
    for(j=0; j < griddy[gridnum]; j++)
      for(i=0; i < griddx[gridnum]; i++){

	// cell index (3d -> 1d conversion - USES FORTRAN ORDERING)
	cellindex = k*griddx[gridnum]*griddy[gridnum] + j*griddx[gridnum] + i;

	if(flagbuffer[cellindex] != -1){ // read just non-overlapped grid cells

	  // calculate center of cell (i,j,k)
	  ccx = gridlex[gridnum] + ( (gridrex[gridnum]-gridlex[gridnum]) / 
				     ((double) griddx[gridnum]) ) 
	    * ( ((double) i) + 0.5 );
	  ccy = gridley[gridnum] + ( (gridrey[gridnum]-gridley[gridnum]) / 
				     ((double) griddy[gridnum]) ) 
	    * ( ((double) j) + 0.5 );
	  ccz = gridlez[gridnum] + ( (gridrez[gridnum]-gridlez[gridnum]) / 
				     ((double) griddz[gridnum]) ) 
	    * ( ((double) k) + 0.5 );

	  // calculate mass of gas in cells
	  cellmass = ((double) densitybuff[cellindex])*cellvol;  // cell mass in code units

	  // calculate overdensity - this is what we're really interested
	  // in, typically.
	  overdensity = ((double) densitybuff[cellindex]) / (OMEGA_BARYON/OMEGA_MATTER);

	  if(multispecies == 2){
	    // H2, Hminus, electron fractions (note funky electron density units)
	    H2frac = ((double) H2Ibuff[cellindex]) / ((double) densitybuff[cellindex]) / 0.76;
	    HMfrac = ((double) Hminusbuff[cellindex]) / ((double) densitybuff[cellindex]) / 0.76;
	    
	  }
	  if(multispecies > 0){
	    efrac = ((double) electronbuff[cellindex]) / ((double) densitybuff[cellindex]);
	    HIIfrac = ((double) HIIbuff[cellindex]) / ((double) densitybuff[cellindex]) / 0.76;
	  }

#ifdef USE_METALS
	  if(metalbuff[cellindex] / densitybuff[cellindex] >= 1.0e-8){
	    metalmff += cellmass;
	    metalvff += cellvol;
	    metalmass += ((double) metalbuff[cellindex])*cellvol;
	    metalvol += ((double) metalbuff[cellindex])*cellvol/((double) densitybuff[cellindex]);
	  }
#endif

	  // calculate entropy assuming gamma = 5/3
	  entropy = ((double) temperaturebuff[i]) * 
	    (
	     pow( ((double) (densitybuff[i]*omegamatter*RHOCRIT_CGS*hubble*hubble)), 
			  ((double) (-2./3.)))
	     / pow( (1.0+redshift),2.0)
	     );

	  // calculate jeans mass in solar masses
	  jeansmass = jeansmass_prefactor * pow( ((double) temperaturebuff[cellindex]), 1.5)
	    / sqrt( ((double) densitybuff[cellindex]) );

	  // Is this the most dense cell we've looked at so far?  If so,
	  // keep track of some information.
	  if(overdensity > maxbaryondens){
	    maxbaryondens = overdensity;
	    maxbar_temp = ((double) temperaturebuff[cellindex]);
	    if(multispecies == 2)
	      maxbar_H2frac = H2frac;
	    if(multispecies > 0)
	      maxbar_HIIfrac = HIIfrac;
	    maxbar_xpos = ccx;
	    maxbar_ypos = ccy;
	    maxbar_zpos = ccz;
	  }

	  // increment total mass, volume: these are the weights for
	  // most of the quantities below.
	  totalmass += cellmass;
	  totalvolume += cellvol;

	  // calculate mass and volume-weighted mean quantities.
	  rhobar_mean_mw += overdensity * cellmass;
	  rhobar_mean_vw += overdensity * cellvol; 

	  rhobarsquared_mean_mw += overdensity*overdensity * cellmass;
	  rhobarsquared_mean_vw += overdensity*overdensity * cellvol;
	  
	  meantemp_mw += ((double) temperaturebuff[cellindex]) * cellmass;
	  meantemp_vw += ((double) temperaturebuff[cellindex]) * cellvol;

	  meanent_mw += entropy*cellmass;
	  meanent_vw += entropy*cellvol;
	  if(multispecies == 2){
	    meanH2frac_mw += H2frac*cellmass;
	    meanH2frac_vw += H2frac*cellvol;

	    meanHMfrac_mw += HMfrac*cellmass;
	    meanHMfrac_vw += HMfrac*cellvol;	  
	  }
	  if(multispecies > 0){
	    meanelecfrac_mw += efrac*cellmass;
	    meanelecfrac_vw += efrac*cellvol;

	    meanHIIfrac_mw += HIIfrac*cellmass;
	    meanHIIfrac_vw += HIIfrac*cellvol;
	  }
	  if( ((double) temperaturebuff[cellindex]) > maxtemp ) 
	    maxtemp = ((double) temperaturebuff[cellindex]);

	  if( ((double) temperaturebuff[cellindex]) < mintemp ) 
	    mintemp = ((double) temperaturebuff[cellindex]);
								  
	  if( overdensity > maxdens) maxdens = overdensity;
	  if( overdensity < mindens) mindens = overdensity;

	  if( entropy > maxent) maxent = entropy;
	  if( entropy < minent) minent = entropy;

	  if(multispecies == 2){
	    if(H2frac > maxH2frac) maxH2frac = H2frac;
	    if(H2frac < minH2frac) minH2frac = H2frac;
	    
	    if(HMfrac > maxHMfrac) maxHMfrac = HMfrac;
	    if(HMfrac < minHMfrac) minHMfrac = HMfrac;
	  }
	  if(multispecies > 0){
	    if(efrac > maxelecfrac) maxelecfrac = efrac;
	    if(efrac < minelecfrac) minelecfrac = efrac;
	    
	    if(HIIfrac > maxHIIfrac) maxHIIfrac = HIIfrac;
	    if(HIIfrac < minHIIfrac) minHIIfrac = HIIfrac;
	  }
	  if(jeansmass > maxjeans) maxjeans = jeansmass;
	  if(jeansmass < minjeans) minjeans = jeansmass;
	  

	  // calculate indices for 1D and 2D distribution functions
	  logodbin = (int) ( (log10(overdensity) - DMIN)/dbin);
	  logtempbin = (int) ((log10( ((double) temperaturebuff[cellindex])) - TMIN)/tbin);
	  logentbin = (int) ((log10(entropy) - SMIN)/sbin);
	  logjeansbin = (int) ((log10(jeansmass) - JMIN)/jbin);
	  if(multispecies > 0){
	    logefracbin = (int) ((log10(efrac) - EMIN)/ebin);
	    logHIIfracbin = (int) ((log10(HIIfrac) - HIIMIN)/HIIbin);
	  }
	  if(multispecies == 2){
	    logH2fracbin = (int) ((log10(H2frac) - H2MIN)/H2bin);
	    logHminusfracbin = (int) ((log10(HMfrac) - HMMIN)/HMbin);
	  }
	  // bounds check indices -- if they're too high/low, just set to the max.
	  if(logodbin < 0) logodbin = 0;
	  if(logodbin > NUMBINS-1) logodbin = NUMBINS-1;

	  if(logtempbin < 0) logtempbin = 0;
	  if(logtempbin > NUMBINS-1) logtempbin = NUMBINS-1;

	  if(logentbin < 0) logentbin = 0;
	  if(logentbin > NUMBINS-1) logentbin = NUMBINS-1;

	  if(logjeansbin < 0) logjeansbin = 0;
	  if(logjeansbin > NUMBINS-1) logjeansbin = NUMBINS-1;
	  
	  if(multispecies > 0){
	    if(logefracbin < 0) logefracbin = 0;
	    if(logefracbin > NUMBINS-1) logefracbin = NUMBINS-1;
	    
	    if(logHIIfracbin < 0) logHIIfracbin = 0;
	    if(logHIIfracbin > NUMBINS-1) logHIIfracbin = NUMBINS-1;
	  }
	  if(multispecies == 2){
	    if(logH2fracbin < 0) logH2fracbin = 0;
	    if(logH2fracbin > NUMBINS-1) logH2fracbin = NUMBINS-1;

	    if(logHminusfracbin < 0) logHminusfracbin = 0;
	    if(logHminusfracbin > NUMBINS-1) logHminusfracbin = NUMBINS-1;
	  }
	  // set all of the 2D distribution functions
	  temp_od_2ddist_mw[logtempbin][logodbin] += cellmass;
	  temp_od_2ddist_vw[logtempbin][logodbin] += cellvol;

	  ent_od_2ddist_mw[logentbin][logodbin] += cellmass;
	  ent_od_2ddist_vw[logentbin][logodbin] += cellvol;

	  temp_ent_2ddist_mw[logtempbin][logentbin] += cellmass;
	  temp_ent_2ddist_vw[logtempbin][logentbin] += cellvol;

	  jeans_od_2ddist_mw[logjeansbin][logodbin] += cellmass;
	  jeans_od_2ddist_vw[logjeansbin][logodbin] += cellvol;

	  if(multispecies == 2){
	    H2frac_od_2ddist_mw[logH2fracbin][logodbin] += cellmass;
	    H2frac_od_2ddist_vw[logH2fracbin][logodbin] += cellvol;

	    HMfrac_od_2ddist_mw[logHminusfracbin][logodbin] += cellmass;
	    HMfrac_od_2ddist_vw[logHminusfracbin][logodbin] += cellvol;

	  

	    efrac_H2frac_2ddist_mw[logefracbin][logH2fracbin] += cellmass;
	    efrac_H2frac_2ddist_vw[logefracbin][logH2fracbin] += cellvol;
	    
	    HMfrac_H2frac_2ddist_mw[logHminusfracbin][logH2fracbin] += cellmass;
	    HMfrac_H2frac_2ddist_vw[logHminusfracbin][logH2fracbin] += cellvol;
	    
	    efrac_HMfrac_2ddist_mw[logefracbin][logHminusfracbin] += cellmass;
	    efrac_HMfrac_2ddist_vw[logefracbin][logHminusfracbin] += cellvol;
	    
	    H2frac_1ddist_mw[logH2fracbin] += cellmass;
	    H2frac_1ddist_vw[logH2fracbin] += cellvol;
	    
	    HMfrac_1ddist_mw[logHminusfracbin] += cellmass;
	    HMfrac_1ddist_vw[logHminusfracbin] += cellvol;
	  }
	  if(multispecies > 0){
	    efrac_od_2ddist_mw[logefracbin][logodbin] += cellmass;
	    efrac_od_2ddist_vw[logefracbin][logodbin] += cellvol;
	  }

	  temp_1ddist_mw[logtempbin] += cellmass;
	  temp_1ddist_vw[logtempbin] += cellvol;

	  od_1ddist_mw[logodbin] += cellmass;
	  od_1ddist_vw[logodbin] += cellvol;

	  ent_1ddist_mw[logentbin] += cellmass;
	  ent_1ddist_vw[logentbin] += cellvol;
	  
	  efrac_1ddist_mw[logefracbin] += cellmass;
	  efrac_1ddist_vw[logefracbin] += cellvol;
	  
	  jeans_1ddist_mw[logjeansbin] += cellmass;
	  jeans_1ddist_vw[logjeansbin] += cellvol;

	} // if(flagbuffer[cellindex] != -1){

      } // end of triply-nested for loop

  if(DEBUG) fprintf(stderr,"leaving GetDensityInfo\n");
  return SUCCESS;
}


/*--------------------- FlagGridCells --------------------------
 *
 *  This function is handed the number of a grid, and there is a
 *  buffer for flags.  We go through each of the OTHER grids, 
 *  only looking for grids that are at level L+1 (assuming this
 *  grid is on level L) and that is overlapping the grid in
 *  question in some way.  If this is true, we go over all of the
 *  cells in the grid and flag each one individually if that cell
 *  happens to be either a) covering a more refined region or
 *  b) if that cell is completely outside the bounds of the
 *  desired projection volume.
 *
 *  The algorithm chosen is probably not the most effective one,
 *  though I can't think of a better one offhand.  If you've
 *  actually taken the time to read the source and can
 *  significantly improve on this, please email me your idea at
 *  bwoshea@cosmos.ucsd.edu, and I'll be very, very grateful.
 *
 *--------------------------------------------------------------*/
int FlagGridCells(int gridnum, int total_number_grids){

  if(DEBUG) fprintf(stderr,"in FlagGridCells: %i %i\n",gridnum,total_number_grids);

  int counter, flagged_cells_this_grid,i,j,k,cellindex;

  double clex, cley, clez, crex, crey, crez;  // cell left and right edges (x,y,z)

  flagged_cells_this_grid=0;  // reset counter to zero

  /* loop over all grids and flag cells that are covered by a grid of
     higher resolution.  We can reject a lot of grids offhand because
     they are:
     1) not of level L+1 (where the grid of interest is level L) or
     2) outside of the bounds of the grid of interest
   */
  for(counter=0; counter < total_number_grids; counter++){
    
    // if the grid of interest is on level l, we only want to
    // look at grids on level l+1 to avoid working too hard
    // (and to make our lives simpler)
    // also takes care of getting rid of grid[gridnum]
    // so we don't flag ALL of the cells!
    if(gridlevel[counter] != (gridlevel[gridnum]+1)) continue;

    // does grid i overlap our grid of interest?  If not, skip it.
    if( !((gridrex[counter] > gridlex[gridnum]) && (gridlex[counter] < gridrex[gridnum]) &&
	  (gridrey[counter] > gridley[gridnum]) && (gridley[counter] < gridrey[gridnum]) &&
	  (gridrez[counter] > gridlez[gridnum]) && (gridlez[counter] < gridrez[gridnum])) 
	) continue;

    // if we've passed all of those tests, grid[counter] is
    //   1)  on level L+1 (where grid[gridnum] is on L)
    //   2)  guaranteed to overlap grid[gridnum] in some way
    // so we loop over every cell in grid[gridnum] and flag all cells that
    // grid[counter] overlaps, and also if there are any cells that are
    // outside the bounds.

    for(i=0; i < griddx[gridnum]; i++){
      for(j=0; j < griddy[gridnum]; j++){
	for(k=0; k < griddz[gridnum]; k++){

	  // calculate cell index
	  cellindex = k*griddx[gridnum]*griddy[gridnum] + j*griddx[gridnum] + i;

	  clex = gridlex[gridnum] + (( (double) i ) * (gridrex[gridnum] - gridlex[gridnum])/
				     ((double) griddx[gridnum]));

	  crex = gridlex[gridnum] + (( (double) (i+1) ) * (gridrex[gridnum] - gridlex[gridnum])/
				     ((double) griddx[gridnum]));

	  cley = gridley[gridnum] + (( (double) j ) * (gridrey[gridnum] - gridley[gridnum])/
				     ((double) griddy[gridnum]));

	  crey = gridley[gridnum] + (( (double) (j+1) ) * (gridrey[gridnum] - gridley[gridnum])/
				     ((double) griddy[gridnum]));

	  clez = gridlez[gridnum] + (( (double) k ) * (gridrez[gridnum] - gridlez[gridnum])/
				     ((double) griddz[gridnum]));

	  crez = gridlez[gridnum] + (( (double) (k+1) ) * (gridrez[gridnum] - gridlez[gridnum])/
				     ((double) griddz[gridnum]));

	  /*  grid[counter] (ie, grid on level L+1) can either only partially fill the 
	     cell in question on grid[gridnum] (grid on level L) or can totally encompass
	     it (much less likely - it would have to be a 2x2x2 grid!) - therefore, see if
	     grid[counter] overlaps grid[gridnum] in ANY way, and flag the cell if it hasn't been
	     flagged before! */

	  if( (gridlex[counter] < crex) && (gridrex[counter] > clex) &&
	      (gridley[counter] < crey) && (gridrey[counter] > cley) &&
	      (gridlez[counter] < crez) && (gridrez[counter] > clez)
	      )
	    // flag cell if it hasn't yet been flagged, and increment counter
	    if(flagbuffer[cellindex] != -1){
	      flagbuffer[cellindex]=-1;
	      flagged_cells_this_grid++;
	    }
	  
	  /* we're going with the approach where, if ANY part of the grid cell overlaps
	     the user-defined projection boundaries, we want to include it - that way
	     there aren't any funky edge effects.  
	     but, only flag if the cell hasn't been flagged already (so as to not screw
	     up the flagged_cells counter)
	  */

	  //if(flagbuffer[cellindex] != -1){
	  //  flagbuffer[cellindex]=-1;
	  //  flagged_cells_this_grid++;
	  //}
	  
	}  // for(k=...
      } // for(j=0...
    }  // for(i=0...

    if(DEBUG) fprintf(stderr,"FlagGridCells:  In Grid: %i (level %i) Grid Looked At: %i (level %i) cells flagged this grid (ttl): %i\n",
		      gridnum,gridlevel[gridnum],counter,gridlevel[counter],flagged_cells_this_grid);
    if(DEBUG) fprintf(stderr,"FlagGridCells: %lf %lf %lf   %lf %lf %lf\n",gridlex[counter],gridley[counter],gridlez[counter],
		      gridrex[counter],gridrey[counter],gridrez[counter]);

  }  // for(counter=0...

  if(DEBUG) fprintf(stderr,"exiting FlagGridCells...\n");
  return SUCCESS;
}


/*----------------- GetGridInfo(int,char) ------------------------
 *
 *   Reads through grid file and retrieves all important information -
 *   grid dimensions, bounds, level, and file name, and puts them all 
 *   into arrays which were allocated in main().
 *
 *---------------------------------------------------------------- */
int GetGridInfo(int numberofgrids,char *hierfilename){
  if(DEBUG) fprintf(stderr,"in GetGridInfo\n");
  FILE *hierfile;
  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];
  int grid=0,i;

  int griddimx,griddimy,griddimz,gridsex,gridsey,gridsez,
    grideex,grideey,grideez,nbfields;
  
  if(DEBUG) fprintf(stderr,"hierarchy file: %s\n",hierfilename);

  // open hierarchy file
  hierfile = fopen(hierfilename,"r");

  // read through hierarchy file, get grid info
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){

    // lines that start with "Grid =" indicate a new set of grid info
    if(strncmp(line,"Grid = ",7)==0){

      fgets(line, MAX_LINE_LENGTH, hierfile);  // junk line
      fgets(line, MAX_LINE_LENGTH, hierfile);  // grid dim
      sscanf(line,"GridDimension     = %d %d %d",
	     &griddimx,&griddimy,&griddimz);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // start index
      sscanf(line,"GridStartIndex    = %d %d %d",
	     &gridsex,&gridsey,&gridsez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // end index
      sscanf(line,"GridEndIndex      = %d %d %d",
	     &grideex,&grideey,&grideez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  // left edge
      sscanf(line,"GridLeftEdge      = %lf %lf %lf",
	     &gridlex[grid],&gridley[grid],&gridlez[grid]);

      fgets(line, MAX_LINE_LENGTH, hierfile);  // right edge
      sscanf(line,"GridRightEdge     = %lf %lf %lf",
	     &gridrex[grid],&gridrey[grid],&gridrez[grid]);
      
      for(i=0;i<3;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
      sscanf(line,"NumberOfBaryonFields = %d",&nbfields);

      if(nbfields==0){
	particleonly=1;

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  // number of particles on this grid

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"ParticleFileName = %s",gridfilenames[grid]);

	//fprintf(stderr,"GetGridInfo:  No baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      } else {
	// "junk" lines
	for(i=0;i<2;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
	
	// grid file name
	sscanf(line,"BaryonFileName = %s",gridfilenames[grid]);

	// "junk" lines
	for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  // number of particles on this grid

	//fprintf(stderr,"GetGridInfo: WITH baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      }

      // grid dims - buffers stripped
      griddx[grid] = 1+grideex-gridsex;
      griddy[grid] = 1+grideey-gridsey;
      griddz[grid] = 1+grideez-gridsez;

      // calculate level from grid bounds, etc.
      gridlevel[grid] = (int) (log10(
			      (1.0 / ((double) rootgridsize)) / 
			      ( (gridrex[grid]-gridlex[grid]) / 
				( (double) griddx[grid] ) )
			      ) / log10(2.0)); 

      grid++;  // increment grid number!
    }
  }

  fclose(hierfile);

  // clean up dynamically allocated stuff
  delete [] line;
  delete [] gfname;

  if(DEBUG) fprintf(stderr,"exiting GetGridInfo\n");

  // return code - number of grids!
  return grid;
}


/*------------------------- NumberOfGrids(char) --------------- 
 *
 *   Reads the hierarchy file and counts the number of grids.
 *   This will be used shortly thereafter to get all of the 
 *   grid info.
 *
 *-------------------------------------------------------------*/
int NumberOfGrids(char *hierfilename){
  if(DEBUG) fprintf(stderr,"in NumberOfGrids\n");

  FILE *hierfile;
  char *line = new char[MAX_LINE_LENGTH];
  int numgrids=0;

  if(DEBUG) fprintf(stderr,"hierarchy file: %s\n",hierfilename);
    
  hierfile = fopen(hierfilename,"r");

  // read through hierarchy file, counting # of grids
  // lines that start with "Grid = " are the beginning of a new
  // piece of grid info
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)
    if(strncmp(line,"Grid = ",7)==0) numgrids++;
    
  fclose(hierfile);

  if(DEBUG) fprintf(stderr,"NumberOfGrids:  there are %i grids\n",numgrids);
  
  // clean up dynamically allocated stuff
  delete [] line;

  if(DEBUG) fprintf(stderr,"exiting NumberOfGrids\n");

  // return # of grids
  return numgrids;  
}



/* ---------------------------- ReadParameterFile -----------------------------
 * 
 *  Reads the parameter file, extracts stuff needed in the rest of the program.
 *
 * ---------------------------------------------------------------------------- */

int ReadParameterFile(char *filename){
  FILE *headerfile;
  char *line = new char[MAX_LINE_LENGTH];
  int ijunk;

  if(DEBUG) fprintf(stderr,"ReadParameterFile:  reading %s\n",filename);

  // open up header file and go through line by line
  headerfile=fopen(filename,"r");

  while( fgets(line, MAX_LINE_LENGTH, headerfile) != NULL){

    // comoving box size (in Mpc/h)
    if(strncmp(line,"CosmologyComovingBoxSize",24)==0){
      sscanf(line,"CosmologyComovingBoxSize   = %lf",
	     &boxsize);
      if(DEBUG){
	fprintf(stderr,"found boxsize:  %s\n",line);
	fprintf(stderr,"box size is: %g\n",boxsize);
      }
    }
    
    // top grid dimension
    if(strncmp(line,"TopGridDimensions",17)==0)
      sscanf(line,"TopGridDimensions   = %d %d %d",
	     &rootgridsize,&ijunk,&ijunk);

    // multispecies
    if(strncmp(line,"MultiSpecies",12)==0)
      sscanf(line,"MultiSpecies                   = %d",
	     &multispecies);

    // omega lambda
    if(strncmp(line,"CosmologyOmegaLambdaNow",23)==0)
      sscanf(line,"CosmologyOmegaLambdaNow    = %lf",
	     &omegalambda);

    // omega matter
    if(strncmp(line,"CosmologyOmegaMatterNow",23)==0)
      sscanf(line,"CosmologyOmegaMatterNow    = %lf",
	     &omegamatter);

    // current redshift
    if(strncmp(line,"CosmologyCurrentRedshift",24)==0)
      sscanf(line,"CosmologyCurrentRedshift   = %lf",
	     &redshift);

    // current time
    if(strncmp(line,"InitialTime",11)==0)
      sscanf(line,"InitialTime         = %lf",
	     &currenttime);


    // hubble constant (in units of 100 km/s/Mpc)
    if(strncmp(line,"CosmologyHubbleConstantNow",26)==0)
      sscanf(line,"CosmologyHubbleConstantNow = %lf",
	     &hubble);


  } // end of while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL)

  fclose(headerfile);

  if(DEBUG) fprintf(stderr,"ReadParameterFile:  root grid is %d\n",rootgridsize);

  return 1;

}

// end of the file, hooray
