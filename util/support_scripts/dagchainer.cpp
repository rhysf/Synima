//  A. L. Delcher            modified by bhaas.
//
//  File:  dagchainer.cpp
//
//  Last Modified:  7 November 2003
//
//  Do DP on dag of matches to get chains of matches


#include  <stdio.h>
#include  <stdlib.h>
#include  <iostream>
#include  <iomanip>
#include  <fstream>
#include  <math.h>
#include  <string.h>
#include  <ctype.h>
#include  <limits.h>
#include  <float.h>
#include  <time.h>
#include  <assert.h>
#include  <errno.h>
#include  <vector>
#include  <algorithm>
#ifdef __sun__
#include <getopt.hh>
#else
#include <getopt.h>
#endif

using namespace std;

// gap extension penalty
float INDEL_SCORE = 1; // must be changed to <= 0

// gap open penalty.
float GAP_OPEN_PENALTY = 1; //must be changed <= 0

// length for a single gap in basepairs:
int BP_GAP_SIZE = -1;  // must be changed > 0


char filename[100] = {'\0'};
int MIN_ALIGNMENT_SCORE = -1;

int MAX_DIST_BETWEEN_MATCHES = 100000; // 100 kb default setting.


static bool  Reverse_Order = false;
  //  If set true by -r option, then use 2nd coordinates
  //  in reverse order

const int  MAX_ERROR_MSG_LEN = 1000;
  // Length of longest possible error message

static int  Max_Y = 0;
  //  The maximum column value of any match.  Used to adjust coords
  //  for reverse diagonals


char Clean_Exit_Msg_Line [MAX_ERROR_MSG_LEN];
  // String to write error messages to before exiting


struct  Cell_t {
  float  raw;
  int  score : 30;
  unsigned  from : 2;
};



struct  Score_t {
  int pairID;  // identifier of match pair
  int  x, y;  // x,y coordinates
  float  score;

  bool  operator < (const Score_t & node)  const {
    return  ( 
	     (x < node.x) 
	     || 
	     (x == node.x && y < node.y)
	     );
  }
};


struct  Path_t {
  float  score;
  int  rc;  // sum of row and column of last entry
  int  sub;
};



bool  Descending_Score (const Path_t & a, const Path_t & b) {
  return  ( 
	   (a.score > b.score)
	   || 
	   (a.score == b.score && a.rc > b.rc)
	   );
}


static void  Print_Chains (vector <Score_t> & score);

FILE* File_Open
    (const char * fname, const char * mode, const char * src_fname = NULL,
     size_t line_num = 0);

void  Clean_Exit
    (const char * msg, const char * src_fname = NULL, size_t line_num = 0);

const char* usage = "\n\n###########################\nusage: dagchainer -G single_gap_length_in_bp -O gap_open_penalty -E gap_extension_penalty -S min_alignment_score -D max_distance_between_matches -F filename [-r]\n\nSpecify -r to reverse 2nd coordinate list\n###########################\n";

int  Best_g = -1, Best_i, Best_j;


void process_arguments (int, char*[]);





int  main (int argc, char* argv[])  {

  FILE* fp;
  Score_t  s;
  vector<Score_t>  score; //list holds all inputted matches.
  vector<Score_t>  alignment;
  Cell_t *space, *p;
  Cell_t** A;
  int  align_ct = 0, max_i = 0, max_j = 0;
  int  f, i, j, k, n, pairID;
  
  // process parameter settings.
  process_arguments (argc, argv);

  // reading in coordinate file, add each match to a list.
  fp = File_Open (filename, "r");
  while  (fscanf (fp, "%d %d %d %f", &s.pairID, &s.x, &s.y, &s.score) == 4) {
    if  (s.y  > Max_Y)
      Max_Y = s.y;
    score.push_back(s); //copy match to score list.
  }
  
  n = score.size();

  if  (n == 0)  {
    fprintf (stderr, "No scores read.  Nothing to do\n");
    return (0);
  }
  

  if  (Reverse_Order)  {
    // reverse complement the second coordinate set.
    for  (i = 0;  i < n;  i++)
      score[i].y = Max_Y - score[i].y + 1;
  }

  sort (score.begin(), score.end());

  for  (i = 1, j = 0;  i < n;  i++) {
    if  (score[j].x == score[i].x && score[j].y == score[i].y)  {
      // consecutive score entries have the same identities.
      fprintf (stderr,
	       "Duplicate score entries:  (%d,%d,%.1f) and (%d,%d,%.1f)\n",
	       score[j].x, score[j].y, score[j].score, 
	       score[i].x, score[i].y, score[i].score);
      fprintf (stderr, "Discarding latter\n");
    } else {
      // not the same
      j++;
      if  (i != j) {
	score[j] = score[i];
      }
    }
  }
  
  j++;
  score.resize (j);
    
  Print_Chains (score);
       
  return  0;
}



static void  Print_Chains (vector<Score_t> & score) {

  //  Find and output highest scoring chains in  score  treating
  //  it as a DAG

  vector <float>  path_score;
  vector <int>  from, ans;
  vector <Path_t>  high;
  Path_t  p;
  bool  done;
  int  ali_ct = 0;
  int  i, j, k, m, n, s;
  
  do {
    done = true;
    n = score.size();
    path_score.resize(n);
    from.resize (n);
    for (i = 0;  i < n;  i ++)  {
      path_score[i] = score[i].score;
      from[i] = -1;
    }
    

    for (j = 1; j < n; j++) {
      for (i=j-1; i >= 0; i--) {
	
	int  del_x, del_y;
	
	del_x = score[j].x - score[i].x - 1;
	del_y = score[j].y - score [i].y - 1;
	
	if  (del_x >= 0 && del_y >= 0)  {

	  if (del_x > MAX_DIST_BETWEEN_MATCHES && del_y > MAX_DIST_BETWEEN_MATCHES) {
	    break;
	  }
	  if (del_x > MAX_DIST_BETWEEN_MATCHES || del_y > MAX_DIST_BETWEEN_MATCHES) {
	    continue;
	  }

	  	  
	  int num_gaps = (int) ( ((del_x + del_y)+ abs(del_x-del_y)) / (2 * BP_GAP_SIZE)  + 0.5);
	  
	  double  x = path_score[i] + score[j].score;
	  
	  if (num_gaps > 0) {
	    // Affine gap penalty:
	    // penalty = open + (num_gaps * extension_penalty)
	    x += GAP_OPEN_PENALTY + (num_gaps * INDEL_SCORE);
	  }
	  
	  
	  if  (x > path_score [j]) {
	    path_score [j] = x;
	    from [j] = i;
	  }
	}
      }
    }
    
    high.clear();
    for  (i = 0;  i < n;  i++) {
      if  (path_score[i] >= MIN_ALIGNMENT_SCORE)  {
	p.score = path_score[i];
	p.sub = i;
	p.rc = score[i].x + score[i].y;
	high.push_back(p);
      }
    }
    
    sort (high.begin(), high.end(), Descending_Score);
    
    m = high.size();
    for  (i = 0;  i < m;  i++) {
      if  (from[high[i].sub] != -2) {
	ans.clear();
	for  (j = high[i].sub;  from[j] >= 0;  j = from[j])  {
	  ans.push_back(j);
	}
	ans.push_back(j);
	if  (from[j] == -2)  {
	  done = false;
	  break;
	} else {
	  reverse(ans.begin(), ans.end ());
	  s = ans.size();
	  printf (">Alignment #%d  score = %.1f\n", ++ali_ct, path_score[high[i].sub]);
	  for  (j = 0;  j < s;  j++) {
	    int  printY;
	    from[ans[j]] = -2;
	    if  (Reverse_Order) {
	      printY = Max_Y - score[ans[j]].y + 1;
	    } else {
	      printY = score [ans[j]].y;
	    }
	    printf ("%3d: %d %6d %6d %7.1f %7.1f\n", 
		    j, score[ans[j]].pairID, score[ans[j]].x,
		    printY, score[ans[j]].score,
		    path_score[ans[j]]);
	  }
	}
      }
    }
    if  (! done)  {
      for  (i = j = 0;  i < n;  i++) {
	if  (from [i] != -2) {
	  if  (i != j)
	    score[j] = score[i];
	  j++;
	}
      }
      score.resize(j);
    }
  } while  (! done);
}



FILE *  File_Open
    (const char * fname, const char * mode, const char * src_fname,
     size_t line_num)

//  Open  fname  in  mode  and return a pointer to its control
//  block.  If fail, print a message and exit, assuming the call came from
//  source file  src_fname  at line  line_num .

  {
   FILE  *  fp;

   fp = fopen (fname, mode);
   if  (fp == NULL)
       {
        sprintf (Clean_Exit_Msg_Line,
                 "ERROR:  Could not open file  %s \n", fname);
        Clean_Exit (Clean_Exit_Msg_Line, src_fname, line_num);
       }

   return  fp;
  }





void  Clean_Exit
    (const char * msg, const char * src_fname, size_t line_num)

//  Write string  msg  to  stderr  and also a line indicating
//  the error happen in source file  src_fname  at line  line_num
//  if they are not  NULL  and  0  respectively.
//  Then exit with an error condition.

  {
   fprintf (stderr, "%s\n", msg);
   if  (src_fname != NULL)
       fprintf (stderr, "  in file  %s", src_fname);
   if  (line_num != 0)
       fprintf (stderr, "  at line  %lu", (long unsigned) (line_num));
   fprintf (stderr, "  errno = %d\n", errno);

   exit (EXIT_FAILURE);
  }




////////////////////////////////////////////////////


void process_arguments (int argc, char* argv[]) {

  int op;
  
  // Option parsing:
  while ((op = getopt(argc, argv, "O:E:F:S:rG:D:")) > 0) {
    switch (op)  
      {
      case 'O':
	GAP_OPEN_PENALTY = atof(optarg);
	break;
	
      case 'E' : 
	INDEL_SCORE = atof(optarg);
	break;
	
      case 'F':
	strcpy(filename, optarg);
	break;
	
      case 'S':
	MIN_ALIGNMENT_SCORE = atoi(optarg);
	break;
	
      case 'r':
	Reverse_Order = true;
	break;
	
      case 'G':
	BP_GAP_SIZE = atoi(optarg);
	break;
	
      case 'D':
	MAX_DIST_BETWEEN_MATCHES = atoi(optarg);
	break;
	
      default:
	fprintf (stderr, "Option %c is not recognized. %s", optopt, usage);
	exit(1); //unrecognizable option.
      }
  }
  
  
  if (argc == 1) {
    fprintf (stderr, "\n%s\n\n", usage);
    exit(1);
  }
  

  if (GAP_OPEN_PENALTY > 0) {
    fprintf (stderr, "The GAP open penalty must be set <= 0\n\n", usage);
    exit(2);
  }
  

  if (INDEL_SCORE > 0) { 
    fprintf (stderr, "The INDEL penalty must be set to number <= 0\n\n\n%s", usage);
    exit(2); // didn't specify the indel penalty.
  }
  if (*filename == '\0') {
    fprintf (stderr, "Must specify filename. %s", usage);
    exit(3); //must specify filename.
  }

  if (MIN_ALIGNMENT_SCORE <= 0) {
    fprintf (stderr, "Must provide a minimum alignment score > 0.  %s", usage);
    exit(4); //must specify min alignment score.
  }

  if (BP_GAP_SIZE <= 0) {
    fprintf (stderr, "Must provide a value for the length of a single gap in basepairs >= 0  using opt -G .  %s", usage);
    exit(5);
  }

  if (MAX_DIST_BETWEEN_MATCHES < 0) {
    fprintf (stderr, "The maximum distance between matches must be >= 0 using opt -D %s", usage);
  }

  
}














