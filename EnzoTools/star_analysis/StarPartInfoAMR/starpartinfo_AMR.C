/*-----------------------------------------------------------

starinfo.C

1) dm/dt vs. z
2) mtotal vs. z
3) creation time vs. metallicity (scatter plot)
4) dynamical time vs. metallicity  (scatter plot)
5) mass vs. metallicity  (scatter plot)

-----------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

char *starfilename = "starinfo.dat";
char *timefilename = "time_vs_redshift.txt";

#define NUMREDSHIFTBINS 199

float redshift[NUMREDSHIFTBINS];
float time[NUMREDSHIFTBINS];

float stellarmass[NUMREDSHIFTBINS-1];
float dmdt[NUMREDSHIFTBINS-1];
float cum_mass[NUMREDSHIFTBINS-1];

#define SPEEDBINS 1000   // each speed bin corresponds to 1 km/s

float speeddisp[SPEEDBINS];

double boxsize = 30.0,    // ComovingBoxSize, in Mpc/h 
  hubble = 0.7,
  gridsize = 128.0,       // number of root grid cells!
  enzo_redshift = 6.52,
  initialredshift = 60.0,
  omegamatter = 0.3;
  




int main(int argc, char *argv[]){

  FILE *starfile, *timefile,*outfile;
  int i,numpart,thisspeedbin;
  float xpos, ypos, zpos, xvel, yvel, zvel,
    crtime, dyntime, metallicity, mass,maxmass,
    minmass,maxmetal,minmetal,speed,minspeed,maxspeed;

  double lengthconv, densconv, timeconv, velconv, massconv;

  lengthconv = boxsize / hubble * 3.08e+24 / (1.0 + enzo_redshift);  // length conversion to proper cm
  
  densconv = omegamatter * 1.8788e-29 * hubble * hubble * pow( (1.0+enzo_redshift), 3.0);  // density in proper g/cm^3

  timeconv = pow( (4.0*3.14159*(6.67259e-8)*omegamatter * 1.8788e-29 * hubble * hubble*
		   pow( (1.0+initialredshift) , 3.0) )
		  , -0.5);   // time conversion in seconds

  velconv = lengthconv/timeconv*(1.0+enzo_redshift)/(1.0+initialredshift);  // velocity in proper cm/s

  massconv = densconv * pow( (lengthconv), 3.0)/1.989e33;  // mass in solar masses  

  for(i=0; i < NUMREDSHIFTBINS; i++) redshift[i] = time[i] = -1.0;
  
  for(i=0; i < NUMREDSHIFTBINS-1; i++) stellarmass[i] = dmdt[i] = cum_mass[i] = 0.0;

  for(i=0; i < SPEEDBINS; i++) speeddisp[i] = 0.0;

  timefile=fopen(timefilename,"r");
  
  i=0;

  while(fscanf(timefile,"%f %f",&time[i],&redshift[i])==2) i++;

  fclose(timefile);

  for(i=0; i<NUMREDSHIFTBINS; i++) printf("%f %f\n",time[i],redshift[i]);

  fprintf(stderr,"boo! (6)\n");

  printf("lengthconv:  %e\n",lengthconv);
  printf("densconv:  %e\n",densconv);
  printf("timeconv:  %e\n",timeconv);
  printf("velconv:  %e\n",velconv);
  printf("massconv:  %e\n",massconv);

  starfile = fopen(starfilename,"r");

  numpart=0;

  while(fscanf(starfile,"%f %f %f %e %e %e %e %f %f %e",&xpos,&ypos,&zpos,
	       &xvel,&yvel,&zvel,
	       &mass,&crtime,&dyntime,&metallicity)==10){

    speed = sqrt( xvel*xvel + yvel*yvel + zvel*zvel);

    if(numpart==0){
      maxmass = minmass = mass;
      maxmetal = minmetal = metallicity;
      minspeed = maxspeed = speed;
    }

    if(speed < minspeed) minspeed = speed;
    if(speed > maxspeed) maxspeed = speed;

    thisspeedbin = ((int) (speed*velconv/1.0e+5));

    if(thisspeedbin < 0) thisspeedbin = 0;
    if(thisspeedbin >= SPEEDBINS) thisspeedbin = SPEEDBINS;

    speeddisp[thisspeedbin] += 1.0;

    if(metallicity < minmetal) minmetal = metallicity;
    if(metallicity > maxmetal) maxmetal = metallicity;

    if(mass < minmass) minmass = mass;
    if(mass > maxmass) maxmass = mass;

    int i=0;

    for(i=0; i<NUMREDSHIFTBINS-2; i++){
      if(time[i] < crtime && crtime <= time[i+1]) stellarmass[i] += mass;
    }

    numpart++;

    if(numpart%100000==0) fprintf(stderr,"%i ptcles read in\n",numpart);
  }

  fclose(starfile);



  printf("%i particles total\n",numpart);
  printf("max, min masses:  %f %f\n",maxmass*massconv,minmass*massconv);
  printf("max, min metallicity:  %e %e\n",maxmetal/0.022,minmetal/0.022);
  printf("max, min speed:  %e %e\n",maxspeed*velconv/1.0e5,minspeed*velconv/1.0e5);
  fflush(stdout);

  outfile = fopen("starinfo_sfr.dat","w");

  for(i=0; i<NUMREDSHIFTBINS-1; i++){
    if(i==0) cum_mass[i] = stellarmass[i]*((float) massconv);
    if(i != 0){
      cum_mass[i] = cum_mass[i-1] + stellarmass[i]* ((float) massconv);
    }

    dmdt[i] = stellarmass[i] / (time[i+1]-time[i]);

    dmdt[i] *= (float) ( (massconv/timeconv) * (365.25*24.0*3600.0) );

    fprintf(outfile,"%f\t%f\t%e\t%e\t%e\t%e\n",
	   (redshift[i+1]+redshift[i])/2.0,
	   (time[i+1]+time[i])/2.0 * (timeconv /(365.25*24.0*3600.0*1.0e6)) ,
	   stellarmass[i]*((float) massconv),
	   cum_mass[i],
	   dmdt[i],
	   dmdt[i]/((float) pow((boxsize/hubble), 3.0)));
  }

  fclose(outfile);

  outfile = fopen("speeddisp.dat","w");
  for(i=0; i<SPEEDBINS; i++)
    fprintf(outfile,"%f %f\n",(float) i,(int) speeddisp[i]);

  fclose(outfile);
}
