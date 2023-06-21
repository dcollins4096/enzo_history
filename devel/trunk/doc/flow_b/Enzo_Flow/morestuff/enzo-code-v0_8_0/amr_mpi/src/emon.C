// ====================================================================
//
// This program reads the enzo output file "OutputLevelInformation.out"
// and displays the information in a user-readable format.
//
// Author: James Bordner (jbordner@cosmos.ucsd.edu)
// Date: 2003-08-22
//
// ====================================================================

#if defined(USE_NCURSES) || defined(USE_CURSES)
#define HAVE_CURSES
#endif

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <vector>
#ifdef USE_NCURSES
#include <ncurses.h>
#endif
#ifdef USE_CURSES
#include <curses.h>
#endif

#ifdef HAVE_CURSES
#  define PRINT printw
#else
#  define PRINT printf
#endif

#define MAX_LINE_LENGTH 1024

//-------------------------------------------------------------
// function prototypes
//-------------------------------------------------------------

void get_input();
void get_line (FILE*);
int get_levels (char *);
bool next (char *buffer, int &i);
float avg_size(int all_cells, int active_cells, int grids);

//-------------------------------------------------------------
// Global (ack!) data
//-------------------------------------------------------------

int g_curr_line = -1;
int g_prev_line = -2;
std::vector<char *> g_lines;
bool g_redraw = true;
int g_dimension = 3;

//-------------------------------------------------------------
// MAIN
//-------------------------------------------------------------

main(int argc, char **argv)
{

  // Get hierarchy depth from the command line

  if (argc != 1) {
    fprintf (stderr,"usage: %s\n",argv[0]);
    exit(1);
  }

  // Open the "OutputLevelInformation.out" file

  FILE *fp = fopen ("OutputLevelInformation.out","r");
  if (fp==NULL) {
    fprintf (stderr,"\n%s: failed to open 'OutputLevelInformation.out'!\n\n",
	     argv[0]);
    exit(1);
  }

  // Initialize curses

#ifdef HAVE_CURSES
  initscr();
  cbreak();
  noecho();
  nodelay (stdscr, TRUE);
  keypad (stdscr, TRUE);
#endif

  int scan_return; 

#ifdef HAVE_CURSES
  werase (stdscr);
#endif

  int i;

  while (true) {

    // Clear screen

    // Read global data

    double time,memory,ratio,sgimem1,sgimem2;
    int depth,grids;

    get_line (fp);

    get_input();

    char buffer[MAX_LINE_LENGTH];

    if (g_curr_line != g_prev_line || g_redraw) {

      g_redraw = false;

      char *buffer = g_lines[g_curr_line];
      int levels = get_levels (buffer);

      i = 0;
#ifdef HAVE_CURSES
  werase (stdscr);
#endif
      sscanf (buffer+i,"%lf",&time);    next(buffer,i);
      sscanf (buffer+i,"%d",&depth);    next(buffer,i);
      sscanf (buffer+i,"%d",&grids);    next(buffer,i);
      sscanf (buffer+i,"%lf",&memory);  next(buffer,i);
      sscanf (buffer+i,"%lf",&ratio);   next(buffer,i);
      sscanf (buffer+i,"%lf",&sgimem1); next(buffer,i);
      sscanf (buffer+i,"%lf",&sgimem2); next(buffer,i);

      //      scan_return = sscanf (buffer,"%lf %d %d %lf %lf %lf %lf",
      //			    &time,&depth,&grids,&memory,&ratio,
      //			    &sgimem1,&sgimem2);

      PRINT ((char *)"========================================================================\n");
      PRINT ((char *)"                               ENZO MONITOR\n");
      PRINT ((char *)"========================================================================\n\n");
      PRINT ((char *)"Step: %d\n",g_curr_line+1);
      PRINT ((char *)"Time: %f\n",time);

      PRINT ((char *)"                         total    active    avg.\n");
      PRINT ((char *)" level  grids  memory    zones     zones    size\n");
      PRINT ((char *)"--------------------------------------------------------------\n");


      int grids_sum=0, all_cells_sum=0, active_cells_sum=0;
      double memory_sum=0.0,coverage_sum=0.0,ratio=0.0;
      double wasted_sum=0.0,average_sum=0.0;

      for (int level = 0; level<levels; level++) {

	int grids, all_cells, active_cells;
	double memory,coverage,ratio,wasted,average;

	sscanf (buffer+i,"%d",&grids);       next(buffer,i);
	sscanf (buffer+i,"%lf",&memory);     next(buffer,i);
	sscanf (buffer+i,"%lf",&coverage);   next(buffer,i);
	sscanf (buffer+i,"%lf",&ratio);      next(buffer,i);
	sscanf (buffer+i,"%d",&all_cells); next(buffer,i);
	sscanf (buffer+i,"%d",&active_cells);next(buffer,i);
	wasted = (float) all_cells / active_cells;
	average = avg_size(all_cells,active_cells,grids);

	grids_sum += grids;
	memory_sum += memory;
	all_cells_sum += all_cells;
	active_cells_sum += active_cells;

	if (grids) {
	  PRINT ((char *)" %3d    %5d  %6.1f  %8d  %8d %6.1f\n",
		 level, grids, memory,all_cells,active_cells,
		 average);
	} else {
	  break;
	}
      }

      wasted_sum  = (float) all_cells_sum / active_cells_sum;
      average_sum = avg_size(all_cells_sum,active_cells_sum,grids_sum);

      PRINT ((char *)"--------------------------------------------------------------\n");
      PRINT ((char *)"Summary %5d  %6.1f  %8d  %8d %6.1f\n",
	     grids_sum, memory_sum,all_cells_sum,active_cells_sum,
	     average_sum);

      g_prev_line = g_curr_line;
#ifdef HAVE_CURSES
      wrefresh (stdscr);
      mvprintw (0,0,"");  // Move cursor to top
#endif
    }

  }
#ifdef HAVE_CURSES
  endwin();
#endif
}

void get_input ()
{
#ifdef HAVE_CURSES
  int ch = wgetch (stdscr);
  if (ch != ERR) {
    switch (ch) {
    case 'q': // EXIT
      endwin();
      exit (0);
    case KEY_LEFT: // previous line
      g_curr_line = g_curr_line - 1 >= 0 ? g_curr_line -1 : 0;
      break;
    case KEY_RIGHT: // next line
      g_curr_line = g_curr_line + 1 < g_lines.size() ? g_curr_line + 1 : g_lines.size()-1;
      break;
    case 'h': // help!
      PRINT ((char *)"(q)uit, (h)elp, (<--) previous, (-->) next,  (2) 2D, (3) 3D\n");
      break;
    case '2': // assume 2D problem
      g_dimension = 2;
      g_redraw = true;
      break;
    case '3': // assume 2D problem
      g_dimension = 3;
      g_redraw = true;
      break;
    default:
      wrefresh (stdscr); 
      break;
    }
  }
#endif
}

void get_line (FILE *fp)
{
  char buffer[MAX_LINE_LENGTH];

  // read line

  int i=0,ch;
  while ((ch = fgetc(fp)) != EOF && ch != '\n') {
    buffer[i++] = ch;
  }
  buffer[i]=0;

  if (strlen(buffer)>0) {

    // store buffer in list of lines

    int s = g_lines.size();
    g_lines.reserve (s+1);
    g_lines.push_back(strdup(buffer));
    if (g_curr_line == s-1) ++g_curr_line;
  }
}

int get_levels (char *buffer)
{
  if (strlen(buffer)>0) {
    int count=0;
    int i=0;
    while (next (buffer,i)) 
      ++ count;
    return (count - 7)/6;

  } else {
    return 0;
  }

}

bool next (char *line, int &i)
  // advance to next field in line
  // return false if we're at the end of the line
{
  int n = strlen(line);
  while (i<n && isspace(line[i])) i++;
  if (i>=n) {
    i = n;
    return false;
  }
  while (i<n && ! isspace(line[i])) i++;
  i = (i < n ? i : n);
  return true;
}

float avg_size(int all_cells, int active_cells, int grids)
{
  // Compute average size X of grids, assuming
  // (X+6)^3 - X^3 = (all_cells - active_cells)
  //  (6 is from two layers of three ghost zones)
  float a = (float)(all_cells-active_cells)/grids;
  if (g_dimension == 3) {
    float d = 2*a/9.0 - 12;
    if (d<=0) return 0.0;
    return 0.5*(-6.0 + sqrt(d));
  } else if (g_dimension == 2) {
    return (a - 36.0) / 12.0;
  }
}
