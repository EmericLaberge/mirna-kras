//
// gcc -lm -std=c99 -O3 db2cm.c -o db2cm
//
// paul dallaire: xlirpu@gmail.com, paul.dallaire@umontreal.ca
//

/*
 * 20150105: change db2bat.c to compute the contact map and output svg
 * 20131010: add: median, mean, variance
 * 20131001: add: compute more stats on base pairs. parameter -mode value 2 and 3
 */

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

const double pi = 3.14159265358979323846;

////////////////////////////////////////// a square 2d matrix data type
struct square_t {
  double *data;
  int side;
};

typedef struct square_t square;

//access a value in square 
double * square_idx(  square *this, int i, int j ){
  if(NULL == this){
    fprintf(stderr, "ERROR: this == NULL in %s at line %i.\n", __FILE__, __LINE__ );
    exit(1);
  }
  
  if( i > this->side || j > this->side ){
    fprintf(stderr, "ERROR: bad index in %s at line %i.\n", __FILE__, __LINE__ );
    exit(1);
  }

  if( NULL == this->data ){
    fprintf(stderr, "ERROR: this data not allocated in %s at line %i.\n", __FILE__, __LINE__ );
    exit(1);    
  }

  return &(this->data[ ( i * this->side ) + j ]); 
}

//allocate a square matrix of doubles
square * square_alloc( int side ){

  square * sq = (square*) malloc( sizeof( square ) );
  if(NULL == sq){
    fprintf(stderr,"ERROR: No more memory in %s at line %i.\n", __FILE__, __LINE__ );
    exit(1);
  }

  sq->data = (double*) malloc( sizeof(double) * side * side );
  if(NULL == sq->data){
    fprintf(stderr, "ERROR: No more memory in %s at line %i.\n", __FILE__, __LINE__ );
    exit(1);
  }

  sq->side = side;

  return sq;
}

//initialize square matrix
void square_zero( square *this ){  
  for( int i=0; i < this->side; ++i )
    for( int j=0; j < this->side; ++j )
      *(square_idx( this, i, j )) = 0.0;
}

square * square_free( square *this ){
  free( this->data );
  free( this );
  return NULL;
}


//////////////////////////////////////// end square data type


//////////////////////////////////////// code to get the median (and other quartiles) rapidly

void swap( double *A, int a, int b ){
  double t = A[a];
  A[a]=A[b];
  A[b]=t;
}

// from sedgewick
int partition( double *A, int l, int r ){

  int i = l-1, j=r;
  double v = A[r];
  
  for(;;){
    while(A[++i]<v);
    while(A[--j]>v & j>l); 
    if(i>=j) break;
    swap(A,i,j);
  }
  swap(A,i,r);
  return i;
}

//worst case is O(n^2), scrambles input
void findKth( double *A, int n, int k ){
  int l=0,r=n-1;
  int kp;
  while( k != (kp=partition( A,l,r )) ){
    if( k>kp ){
      l=kp+1;
    } else if( k<kp ) {
      r=kp-1;
    } else {
      return;
    }
    if( l>=(r+1) ){
      return;
    }
  }
}


double median( double *A, int n ){
  if( n<1 )
    //should return nan but this screws the normalization up
    //return nan("char-sequence");
    return 0.0;
  findKth(A,n,n>>1);
  return A[n>>1];
}

double q1( double *A, int n ){
  if( n<1 )
    return 0.0;
  findKth(A,n,n>>2);
  return A[n>>2];
}

//////////////////////////////////////// two pass algorithm to calculate the variance (take sqrt for s.d.) 
//////////////////////////////////////// numerical recipies in c 2nd edition section 14.1 (eq 14.1.8) code adapted from moment(...)
//////////////////////////////////////// this algo minimizes roundoff errors for variance when x-<x> is small

typedef struct stats_t {
  double sum,
    max,
    q1, //first quartile (k==n/4)
    median, 
    mean, //average: <x>
    var, //normalized using n-1 
    sd, //sqrt(var)
    adev, // sum( abs( x - <x> ) ) absolute deviation to mean
    mad; //aboslute deviation to median
} stats;

void zero_stats(stats *m){
  m->sum=0.0;
  m->max=0.0;
  m->q1=0.0;
  m->median=0.0;
  m->mean=0.0;
  m->var=0.0;
  m->sd=0.0;
  m->adev=0.0;
  m->mad=0.0;
}


//tested on vectors of length 10000 of random data [0..1] using code testStats.c
//the result was compared to similar calculations using built-in functions in R
//There is a numerical imprecision in this process due to the format of the data 
//outputed which was fixed for vector data. Overall behaviour is perfect.
stats get_stats( double *A, int n ){
  stats m;
  zero_stats(&m);
  //should return nan but this screws the normalization up
  if( n<1 ) return m;
  //if( n<1 ) return nan("char-sequence");
  if( n<2 ) return m;

  //1- compute the mean
  double mean = 0.0, var=0.0, sum=0.0, max=0.0;

  for(int i=0;i<n;++i){
    sum += A[i];
    max = A[i]>max?A[i]:max;
  }
  
  mean = sum / n;

  
  double mid = median(A,n);
  double mad = 0.0; //absolute deviation to median
    
  //2- compute eq 14.1.8
  double s1=0.0, s2=0.0, adev = 0.0;
  for(int i=0;i<n;++i){
    mad += fabs( A[i] - mid );
    double s = A[i] - mean;
    adev += fabs( s );
    s1 += s*s;
    s2 += s;
  }
  mad /= n;
  adev /= n;

  double variance = (s1 - (s2*s2/n))/(n-1);
  
  m.sum = sum;
  m.max = max;
  m.q1 = q1(A,n);
  m.mean = mean;
  m.var = variance;
  m.sd = sqrt(variance);
  m.adev = adev;
  m.median = mid;
  m.mad = mad;
  return m;
}


//////////////////////////////////////// state machine used in db2bat
//////////////////////////////////////// to convert dot brackets to sequence of 
//////////////////////////////////////// integers

enum states { s1=0, s2, s3, s4, s5, s6, sSTART, sSTOP, 
	      s35, s467, sERROR, sEND,  };

const char trans[10][4] = {
  { s1,     s35,    s2,      sERROR  }, /* s1 */
  { s1,     s467,   s2,      sSTOP   }, /* s2 */
  { s1,     s3,     sERROR,  sERROR  }, /* s3 */
  { sERROR, s4,     s2,      sERROR  }, /* s4 */
  { sERROR, s5,     s2,      sERROR  }, /* s5 */
  { s1,     s6,     sERROR,  sERROR  }, /* s6 */
  { s1,     sSTART, sERROR,  sSTOP   }, /* sSTART */
  { sERROR, sSTOP,  sERROR,  sSTOP   }, /* sSTOP */
  { s1,     s35,    s2,      sERROR  }, /* s35 */
  { s1,     s467,   s2,      sSTOP   } /* s467 */
};



//////////////////////////////////////// accessory functions

double max( double a, double b ){ return a<b?b:a; }

void db_error( char *db, int i ){
  fprintf( stderr, "ERROR parsing dot bracket structure at position %i.\n%s\n", i, db );
  exit(1);
}


//////////////////////////////////////// stuff to parse a dot bracket

int convert_db_char( char c ){
  switch( c ){
  case '(' : return 0; break;
  case '.' : return 1; break;
  case ')' : return 2; break;
  default: return 3;
  }
}
/*
 * returns the (newly malloced) vector of indices in the batman offset 
 * 0:forward pairing, 1:reverse pairing, 2:forward bulge, 3: reverse bulge, 
 * 4: terminal hairpin loop, 5:unpaired in a multibranch, 6: dangling free ends
 * 7: gap, 8: other <- 7 and 8 are not processed.
 * 
 * db contains nothing but ()., no score etc
 */
int * db2bat( char *db, square *basepairs, double weight ){

  int len = (int)strlen(db);
  int *dbo = (int*)malloc( len * sizeof(int) );
  bzero( dbo, len * sizeof(int) );


  //use a stack to determine base pairs for the max version of the program
  //collect the base pairs data in basepairs and main will decide wether to use the info
  //or else use the sum
  
  int stack[len];
  int sp = -1;
  for( int i=0;i<len;++i ){
    if( db[i] == '(' ){
      //push
      stack[++sp]=i;
    } else if( db[i] == ')' ){
      //pop
      *(square_idx( basepairs, stack[sp--], i )) += weight;
    } else { //db[i] == ',' 
      ;
      *(square_idx( basepairs, i, i )) += weight;
    }
  }

  //run the state machine forward on the db.
  // maybe this is overkilled a bit...
  int state = sSTART;
  for( int i=0; i<len; ++i){
    dbo[i] = state = trans[state][convert_db_char(db[i]) ];
    if( sERROR == state )
      db_error( db, i );
    if( sSTOP == state )
      break;
  }

  //run the db backwards to spread knowledge about type of unpaired nucleotides
  for( int i=len-1; i>=0; --i ){
    if( dbo[i] == s467 ){
      if(i==len-1){
	dbo[i]=sSTART;//sSTOP;
      } else {
	switch(dbo[i+1]){
	case s1: dbo[i]=s6;break;
	case s2: dbo[i]=s4;break;
	default: dbo[i]=dbo[i+1]; 
	}
      }
    } else if( dbo[i] == s35 ){
      switch(dbo[i+1]){
      case s1: dbo[i]=s3;break;
      case s2: dbo[i]=s5;break;
      default: dbo[i]=dbo[i+1]; 	
      }  
    } 
  } 

  return dbo;
}  

int isCanonical( char a, char b ){
  a=toupper(a);
  b=toupper(b);
  b=(b=='T'?'U':b);
  return 
    (a=='A' & b=='U') | (a=='C' & b=='G' ) | (a=='G' & b=='C') | (a=='U' & b=='A');
}

//////////////////////////////////////// Generate .svg file showing the contact map

//provide: imgWidth_cm, imgHeight_cm, boxLeft, boxTop, boxWidth, boxHeight, imageDescription

//"<?xml version=\"1.O\" standalone=\"no\" ?>		\
//<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"	\
//  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">	\
//

const char* SVGHEADER = 
  "<svg width=\"%fcm\" height=\"%fcm\" viewBox=\"%f %f %f %f\"\
     xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\
     <desc>run of db2cm on %s</desc>\n";



//dbnum is number of 2D structures parsed 
void printSVG( square * basepairs, char *seq, double normalisation, 
	       double sensitivity, //typical: 4 3 2 1 0 -1 -2 -3 -4
	       double rotation, //degrees
	       double zoom,     //ratio of enlargement of plot within the right panel
	       double size,     //size of a nucleotide in cm (default should probably be about 0.2)
	       double threshold, //value below which no basepair is plotted (typical: 0.0001)
	       char * svgFileName, //name of file to generate
	       int isBoltzmann,    //set to 1 if the basepairs matrix was computed using Boltzmann weights 
	                         //or to 0 if it was computed by incrementing values with fixed values
	       int dbnum,       //number of dot bracket predictions (used solely for reporting to user in left panel's legend)
	       char *basename,  //name to report on  left panel and in .svg file description
	       char *extraUserdata, //data added to legend as text (few words only, no bounds checking)
	       int plot_arcs,           //do we plot the arc plot
	       int plot_dotPlot,        //do we plot the dot plot
	       int plot_seqOnDiagonal,  //do we draw the sequence on the diagonal
	       int plot_seqOnDotPlot,   //do we draw the seqeunces on the sides of the dot plot
	       int plot_legend          //do we plot the left panel
	       ){

	
  rotation = fmod(rotation,360.0); //probably not strictly necessary

  threshold = threshold > 1.0 ? 1.0 : threshold; //clamp threshold
  threshold = threshold < 0.0 ? 0.0 : threshold;

  zoom = zoom < 0.0 ? 0.0 : zoom; //clamp zoom 
	
  FILE *outf = stdout;
  if( svgFileName != NULL ){
    outf = fopen( svgFileName, "w+"); // if exists: wipe!
  }
  if( outf == NULL ) //fallback
    outf = stdout;

  int seqLen = basepairs->side; //length of sequence could also be take from stringlenght of seq
  const int pixelSz = 3; //size of a base pair dot in the image of the contact map matrix 
  const int borderWidth=1;  //thickness of black border line drawn around the image right and left widths
  const int borderHeight=1; //thickness of black border line drawn around the image top and bottom widths
	
  // Legend drawing switch. viewbox pixel rapidly becomes ridiculously too large for inteligent plot drawing at large sizes
  int leftPanel = plot_legend;
	
	
  //the viewbox has logical (0,0) at the top left of the dot plot Matrix (rotated and scaled afterwards)
  double bw =  (pixelSz * seqLen + 2.0*borderWidth );  //width of the viewbox
  double bh =  (pixelSz * seqLen + 2.0*borderHeight ); //height of the viewbox
	
  double legend_width_in_cm = leftPanel ? 5.0 : 0.0;
  double legend_width = leftPanel ? legend_width_in_cm * bw / (seqLen * size) : 0.0;
	
	
	
	
  const double pseudoCMperCM = 0.995;  //adjustement factor so that cm are cm in illustrator as determined by its ruler
  fprintf( outf, SVGHEADER,
	   (legend_width_in_cm + (double)seqLen * size ) / pseudoCMperCM,// /10.0, //imgWidth_cm, 
	   ( (double)seqLen * size ) / pseudoCMperCM, //10.0, //imgHeight_cm, 
	   -1.0 * legend_width, //boxLeft, 
	   0.0, //boxTop, 
	   bw + legend_width, //boxWidth, 
	   bh, //boxHeight, 
	   basename //commentary
	   );


  //the image is draws as a square with contact map upper right and mountain plot lower left
  // but the image gets rotated (180deg) to present the contact map in lower left of image
  // and the mountain lines in the top right
  // but It is possible to generate a mushroom image by rotating 135 degrees 
  // in this case, we must scale by 1/hypothenuse so that the square fits in the viewbox
  double sc = cos(pi/4.0)/cos(pi*(45.0-fabs(fmod(rotation,90.0)))/180.0); //scaling factor
  sc *= zoom; //user specified zoom
	
	
  fprintf( outf, 
	   "<g transform = \"rotate( %f, %f, %f),translate(%f,%f),scale(%f)\" >\n",
	   rotation,
	   bw/2.0, bh/2.0, // rotate :horizontal center, vertical center
	   (1-sc)*bw/2.0, (1-sc)*bh/2.0, //translate for scaling
	   sc //scaling factor
	   );

	
  //draw right panel
  fprintf( outf, "<g transform = \"translate( 1, 1 )\">\n" );

  //draw dot plot
  if( plot_dotPlot ) {
    for( int r = 0; r < seqLen; ++r ){
      for( int c = r; c < seqLen; ++c ){
	double freq = *(square_idx( basepairs, r, c ));
	freq /= normalisation;
	
		  if(freq > threshold ){ // (1.0-threshold)){
		freq = pow( freq, pow(2.0,-1.0*sensitivity) );
		freq = 100.0 * (1.0-freq); //now, freq is a color intensity actualy: lower values print darker

	
		  
	  //dots in the dotPlot (upper right part of the matrix)
	  // these are hollow squares if the base pair is canonical and 
	  // filled (sligthly smaller) circles if the base pair is not cannonical.
	  if( isCanonical(seq[r],seq[c]) ) {
	    //output a gray square (pixel) whose color density is given by freq
	    fprintf( outf,  
		     "<g><rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" fill=\"rgb(%f%%, %f%%, %f%%)\" stroke=\"rgb(%f%%, %f%%, %f%%)\" stroke-width=\"%i\" />\n",
		     c*3, //x, 
		     r*3, //y, 
		     pixelSz, //width, 
		     pixelSz, //height, 
		     freq, freq, freq, // fill ie: rgb(33%, 33%, 33%) 
		     freq, freq, freq, // stroke ie: rgb(33%, 33%, 33%) 
		     //"black", //stroke, 
		     0 //strokeWidth 
		     );
	  
	    fprintf( outf,  
		     "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" fill=\"white\" stroke=\"none\" /></g>\n",
		     c*3+1, //x, 
		     r*3+1, //y, 
		     pixelSz-2, //width, 
		     pixelSz-2 //height, 
		     );
	  } else { 
	    fprintf( outf,  
		     "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"rgb(%f%%, %f%%, %f%%)\" stroke=\"rgb(%f%%, %f%%, %f%%)\" stroke-width=\"%i\" />\n",
		     c*3+1.5, //cx, 
		     r*3+1.5, //cy, 
		     pixelSz/2.0, //r
		     freq, freq, freq, // fill ie: rgb(33%, 33%, 33%) 
		     freq, freq, freq, // stroke ie: rgb(33%, 33%, 33%) 
		     //"black", //stroke, 
		     0 //strokeWidth 
		     );
	  }
	}
      }
    }
    
    //draw sequences on the side of the dot plot
    if( plot_seqOnDotPlot ) {
      fprintf(outf,"<g>\n");
      fprintf( outf, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"black\" />\n",
	       0.0, -3.0, 
	       seqLen*3.0, 3.0 );
	
      for( int r = 0; r < seqLen; ++r ){
	fprintf(outf, 
		"<g transform=\"translate(%f,%f),rotate(%f,-0.5,0.5)\"><text x=\"-1.5\" y=\"1.5\" font-size=\"3\" fill=\"%s\">%c</text></g>\n", 
		(double)r*3+2,//1.5, //translate x
		-2.0,//(double)(r+1)*3-1.5, //translate y
		-1.0*rotation, 
		"white",//fill color
		seq[r] //char to print
		);
      }
      fprintf(outf,"</g>\n");

      fprintf(outf,"<g>\n");
      fprintf( outf, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"black\" />\n",
	       (seqLen)*3.0, 0.0, 
	       3.0, seqLen*3.0 );
      for( int r = 0; r < seqLen; ++r ){
	fprintf(outf, 
		"<g transform=\"translate(%f,%f),rotate(%f,-0.5,0.5)\"><text x=\"-1.5\" y=\"1.5\" font-size=\"3\" fill=\"%s\">%c</text></g>\n", 
		3.0*seqLen+2.0,//(double)r*3+2,//1.5, //translate x
		(double)(r)*3.0+1.0, //translate y
		-1.0*rotation, 
		"white",//fill color
		seq[r] //char to print
		);
		
      }
      fprintf(outf,"</g>\n");
    }
  }


  //draw the arcs plot
  if( plot_arcs ) {
    fprintf(outf,"<g>\n");
    for( int r = 0; r < seqLen; ++r ){
      for( int c = r; c < seqLen; ++c ){
	double freq = *(square_idx( basepairs, r, c ));
	freq /= normalisation;
		  
		  if(freq > threshold) {// (1.0-threshold)){
		
		freq = pow( freq, pow(2.0,-1.0*sensitivity) );
		freq = 100.0 * (1.0-freq); //now, freq is a color intensity actualy: lower values print darker

		
	  // plot of the corresponding arc in the lower left part of the matrix
	  fprintf( outf,  
		   "<path d=\"M %f %f Q %f,%f %f,%f\" fill=\"none\" stroke=\"black\" opacity=\"%f\" stroke-width=\"1\" />\n",
		   //stroke=\"rgb(%f\%, %f\%, %f\%)\" stroke-width=\"1\" />\n",
		   (c-1)*3+1.5, //from x 
		   (c+1)*3+1.5, //from y
		   (r-1)*3+1.5 + (r-c)*3/4, //quadratic x
		   (c+1)*3+1.5 - (r-c)*3/4, //quadratic y
		   (r-1)*3 +1.5, //to x
		   (r+1)*3 +1.5, //to y
		   //freq, freq, freq // stroke ie: rgb(33%, 33%, 33%) 
		   (100.0-freq)/100.0
		   );
	}
      }
    }
    fprintf(outf,"</g>\n");
  }
	
    //draw a rectangle filled with the character of the corresponding nucleotide character on the main diagonal of the matrix
  // the rectangle is a square of the matrix and thus rotates with the matrix
  // the text in each square is meant to be read and so it does rotate with the matrix
  // in a way that it is always possible to read the nucleotide when viewing the image.
  // the square is black and the text is white.
  if( plot_seqOnDiagonal ){
    fprintf(outf,"<g>\n");
    for( int r = 0; r < seqLen; ++r ){

      fprintf(outf, 
	      "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" fill=\"BLACK\" /> <g transform=\"translate(%f,%f),rotate(%f,-0.5,0.5)\"><text x=\"-1.5\" y=\"1.55\" font-size=\"3\" fill=\"%s\">%c</text></g>\n", 
	      (r-1)*3.0-0.5,(r+1)*3.0-0.5,4.0,4.0, //rect data
	      (double)(r-1)*3+2,//1.5, //translate x
	      (double)(r+2)*3-2, //translate y
	      -1.0*rotation, 
	      "white",//fill color
	      seq[r] //char to print
	      );
    }
    fprintf(outf,"</g>\n");
  }
	
	
  fprintf(outf,"</g></g>\n"); //that's it for the right panel

  //draw the left panel
  if( leftPanel ) {
		
    //overwrite band at left (leftPanel) of image where the sensitivity plot and other data will be drawn
    // * note: clipping does not work properly on groups (both in illustrator and in chrome)
    //		fprintf( outf, "<rect x=\"%f\" y=\"0\" width=\"%f\" height=\"%f\" fill=\"yellow\" />\n", -1.0*legend_width, legend_width, bh );
    fprintf( outf, "<rect x=\"%f\" y=\"0\" width=\"%f\" height=\"%f\" fill=\"white\" />\n", -1.0*legend_width, legend_width, bh );
			
		
    //sensitivity curve legend position
    double slside = legend_width * 0.6;
    double slleftspace = legend_width * 0.2;
    double slbottomlocation = slleftspace + slside + (bh*(1.0-(zoom>1.0?1.0:(zoom<0.5?0.5:zoom)))/2.0);//bh - (2.0*slleftspace);
    double fontsize = 0.3 * legend_width/legend_width_in_cm;
    double lineWidth = slside * 0.02;
		
    //move to the leftPanel
    fprintf(outf, "<g transform=\"translate(%f,0)\">\n", 
	    -1.0*legend_width //move to leftPanel location 
	    );
		

    //move to sensitivity graph x-axis decoration location
    fprintf(outf, "<g transform=\"translate(%f,%f)\">\n", 
	    slleftspace, slbottomlocation
	    );
		
    //x axis decoration (gray color bar at bottom)
    //draw rectangles side by side with proper color
    //draw two lines of these rectangles with an offset in their positions
    // so as to fill any gap between rectangles (drawing is not perfect)
    double graybarHeight = slside * 0.2;
    double barwidth = 0.05/size;
    double overlap = barwidth;
    for( double x=1.5*barwidth;x<slside;x+=barwidth ){
      double col=100.0 - 100.0 * ((double)x / (double)slside);
      fprintf( outf, "<rect x=\"%f\" width=\"%f\" y=\"%f\" height=\"%f\" fill=\"rgb(%f%%, %f%%, %f%%)\"  stroke=\"none\" />\n", 
	       x-barwidth, barwidth*0.8, 
	       2*lineWidth, graybarHeight,
	       col, col, col //opacity
	       );
    }			
    for( double x=barwidth;x<slside;x+=barwidth ){
      double col=100.0 - 100.0 * ((double)x / (double)slside);
      fprintf( outf, "<rect x=\"%f\" width=\"%f\" y=\"%f\" height=\"%f\" fill=\"rgb(%f%%, %f%%, %f%%)\"  stroke=\"none\" />\n", 
	       x-barwidth, barwidth, 
	       2*lineWidth, graybarHeight,
	       col, col, col //opacity
	       );
    }
		
		
    fprintf(outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >levels of gray</text>\n",
	    0.0,//x
	    2*lineWidth+graybarHeight+fontsize,//y
	    fontsize//font-size
	    );
		
		

    fprintf(outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" transform=\"rotate(-90)\" >basepair probability</text>\n",
	    0.0,//x (rotated becomes y)
	    -2.0 * fontsize,//y (rotated becomes x)
	    fontsize//font-size
	    );
		
		
    fprintf( outf, "</g>\n" ); //close move to sensitivity graph x-axis decoration location
		
		
    fprintf(outf, "<g transform=\"translate(%f,%f),scale(%f,%f)\">\n", 
	    slleftspace, slbottomlocation, //move to graph location 
	    1.0, -1.0 //invert y axis for ease (not particularly usefull as it turns out)
	    );

		
    //write a 0 on the y-axis
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >0</text>\n", 
	     -1.0*fontsize,//x
	     0.0+0.5*fontsize,//y
	     fontsize//fontsize
	     );
		

    fprintf( outf, "<rect x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\" stroke=\"lightgray\" stroke-width=\"%f\" fill=\"lightgray\" transform=\"scale(1,-1)\" />youp\n", 
	     0.0,//x
	     -1.0*threshold*slside,//y
	     slside*pow(threshold,pow(2,-1.0*sensitivity)),//width
	     threshold*slside,//height
	     0.5*lineWidth
	     );
		
    //put a dotted line at threshold on the y-axis
    const int segments=21;
    for( double i=0;i<segments;i+=2){
      fprintf( outf, "<line x1=\"%f\" x2=\"%f\" y1=\"%f\" y2=\"%f\" stroke=\"black\" stroke-width=\"%f\" />\n", 
	       i/segments*slside, (i+1)/segments*slside, slside*threshold, slside*threshold, lineWidth);
    }
    //write a T on the right y-axis
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" transform=\"scale(1,-1)\">T</text>\n", 
	     0.5*fontsize + slside,//x
	     -1.0*(-0.2*fontsize + threshold * slside),//y
	     fontsize//fontsize
	     );

		
		
		
    //write a 1 on the y-axis
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" transform=\"scale(1,-1)\">1</text>\n", 
	     -1.0*fontsize,//x
	     -1.0*(0.0-0.5*fontsize+slside),//y
	     fontsize //fontsize
	     );
		
		
    //y axis of the sensitivity plot
    fprintf(outf, "<line x1=\"0\" x2=\"0\" y1=\"0\" y2=\"%f\" stroke=\"black\" stroke-width=\"%f\"/>\n", slside, lineWidth);
    //x axis of the sensitivity plot
    fprintf(outf, "<line x1=\"0\" x2=\"%f\" y1=\"0\" y2=\"0\" stroke=\"black\" stroke-width=\"%f\"/>\n", slside, lineWidth );
		
    //this polyline is the data curve in the sensitivity plot
    double y=0.0;
    const double step=100;
    fprintf( outf, "<polyline stroke=\"black\" stroke-width=\"%f\" fill=\"none\" points=\"\n", lineWidth);
    fprintf(outf,"%f,%f ", slside * pow(y,pow(2,-1.0*sensitivity)), y*slside );
    y += 1.0/step;
    for( int i=1;i<=step;++i ){
      double x = pow(y,pow(2,-1.0*sensitivity));
      double y1 = y+1.0/step;
      double x1 =	pow(y1,pow(2,-1.0*sensitivity));
      fprintf( outf, "%f,%f ", x*slside, y*slside);
      y = y1;
    }
    fprintf(outf,"\" />\n");
		
    fprintf( outf, "</g>\n" ); //close move ot graph location
	
    double y0 = slbottomlocation + 2.0* slleftspace;
    double leftMargin = slleftspace / 2.0;
    const double VSPACE = 1.4; 
    //write name of input file (basename)
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >input: %s</text>\n", leftMargin, y0, fontsize, basename );
		
    y0 += VSPACE * fontsize;
		
    //write wether Boltzmann was used to weight this
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >dot brackets are %s weighted</text>\n", leftMargin, y0, fontsize, isBoltzmann?"Boltzmann":"linear" );
		
    y0 += VSPACE * fontsize;
    //write number of dot brackets in the run ?
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >%i structures</text>\n", leftMargin, y0, fontsize, dbnum );

    y0 += VSPACE * fontsize;
    //write sensitivity
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >sensitivity=%.1f</text>\n", leftMargin, y0, fontsize, sensitivity );
		
    y0 += VSPACE * fontsize;
    //write sensitivity
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >threshold=%.2e</text>\n", leftMargin, y0, fontsize, threshold );
		
    //write user extra annotation data
    y0 += VSPACE * fontsize;
    if( extraUserdata != NULL ){
      fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >%s</text>\n", leftMargin, y0, fontsize, extraUserdata );
    }
	  
	  
    //write trademark
    y0 = bh - (1.5 * fontsize);
    fprintf( outf, "<text x=\"%f\" y=\"%f\" font-size=\"%f\" >Generated using db2cm</text>\n", leftMargin, y0, fontsize );
	  
    fprintf( outf, "</g>\n" ); //close move to leftPanel location
	
  }
  fprintf( outf, "</svg>" );

  fflush( outf );
  if( outf != stdout )
    fclose( outf );
}



//these should likely be made dynamic (maybe latter)
const int LONGFILENAME = 10000;
const int LONGSEQTITLE = 1000;
const int LONGEST_SEQ = 10000;

//batman signature constant (nothing to do with dot plots)
const int ROWS = 9;

//gaz constant times temperature for Boltzmann
const double RT = 0.61597;



///////////////////////////////////////// Parameter type (command line parameters)
enum datatype { NUMBER, TEXT };

struct Ptype  {
  const char *name;
  const char *explanation;
  const enum datatype type;
  double *value;
  char ** text;//yeah, I know...
};

//////////////////////////////////////// Variables that will be filled at command line parsing
double xdb = -1;
double explore = 1e100;
double use_boltzmann = 0;
double weight = 1000.0;

//these are exclusively for batman signature
enum mode_Value_t { SUM=0, MAX, Q1, MEDIAN, MEAN, STDDEV, ADEV, MAD  };
double mode = 0;
double mode2 = -1;
double output_strNum = 0;
double silent_mode = 1;
double stats_mode = 0;
double explicit_w = -1; 

//these are exclusively for dot and arcs plots
double dotPlot_sensitivity = 0;
double dotPlot_rotation = 0.0;
char * dotPlot_fnameExt = "";
char * dotPlot_annotation = "";
double dotPlot_zoom = 0.9;
double dotPlot_size = 0.2;
double dotPlot_threshold = 0.0001;
double dotPlot_arcs_OFF = 1.0;

//turn off features of dot and arcs plots
double dotPlot_plot_arcs = 1.0;
double dotPlot_plot_dotPlot = 1.0;
double dotPlot_plot_seqOnDiagonal = 1.0;
double dotPlot_plot_seqOnDotPlot = 1.0;
double dotPlot_plot_legend = 1.0;

char * inputfilename = ""; //can be used for one file instead of reading filenames from stdin

 
///////////////////////////////////////// command line parameters as Ptype values
///////////////////////////////////////// (must end with {0,0,0})
struct Ptype parameters[] = {
  { "-xdb", "Maximum number of dot brackets to accept.", NUMBER, &xdb, 0 },
  { "-e", "Max fraction of MFE to accept (ppt:0..1000).", NUMBER, &explore, 0 },
  { "-b", "Should score be boltzmann based (don't use if energies are too big) (0:no, 1:yes)(default: 0)", NUMBER, &use_boltzmann, 0 },
  { "-w", "Weight of probability to apply (ppt:1..1000).", NUMBER, &weight, 0 },
  { "-mode", "Set the base pair stats collection mode (0:sum, 1:max, 2:Q1, 3:median, 4:mean, 5:sqrt(var), 6:ADEV, 7:ADEV2, 100+: dump)(default is sum).", NUMBER, &mode, 0 },
  { "-mode2", "Secondary base pairs stats (-1:zero, 0:sum, 1:max, 2:Q1, 3:median, 4:mean, 5:sqrt(var), 6:ADEV, 7:ADEV2, 100+: dump)(default is NONE).", NUMBER, &mode2, 0 },
  { "-o", "Output number of structures considered (0:no (default), 1:stdout, 2:stderr)", NUMBER, &output_strNum, 0 },
  { "-s", "Silent mode. Do not output signatures. (0:no (default), 1: yes, shut up)", NUMBER, &silent_mode, 0 },
  { "-h", "Help message", 0, 0 },
  { "-sensitivity", "Dot plot sensitivity 0=normal, higher=more base pairs, lower=less base pairs", NUMBER, &dotPlot_sensitivity, 0 },
  {"-rotation", "Dot plot image rotation in degrees", NUMBER, &dotPlot_rotation, 0 },
  {"-extra", "Dot plot filename extra annotation", TEXT, 0, &dotPlot_fnameExt },
  {"-annot", "Dot plot. User annotation to be printed in figure legend.", TEXT, 0, &dotPlot_annotation }, 
  {"-zoom", "Dot plot image zoom factor (default is 0.9)", NUMBER, &dotPlot_zoom, 0 },
  {"-size", "Dot plot image size scaling factor (1 is about 1cm/nt. Default is 0.2)", NUMBER, &dotPlot_size, 0 },
  {"-threshold", "Dot plot. Do not plot basepairs whose probability is this much or lower (default is 0.0001)", NUMBER, &dotPlot_threshold, 0 },
  {"-arcs", "Control drawing of arcs plot (Default = 1, set to 0 to block.)", NUMBER, &dotPlot_plot_arcs, 0},
  {"-dotPlot", "Control drawing of dot plot (Default = 1, set to 0 to block.)", NUMBER, &dotPlot_plot_dotPlot, 0},
  {"-seqd", "Control drawing of sequence on main diagonal (Default = 1, set to 0 to block.)", NUMBER, &dotPlot_plot_seqOnDiagonal, 0},
  {"-seqs", "Control drawing of sequences on side of dot plot (Default = 1, set to 0 to block.)", NUMBER, &dotPlot_plot_seqOnDotPlot, 0},
  {"-legend", "Control drawinf of legend elements at left of plots (Default = 1, set to 0 to block.)", NUMBER, &dotPlot_plot_legend, 0},
  {"-f", "If you have only one input file, you can name it using this parameter", TEXT, 0, &inputfilename},
	{ 0,0,0 }
};


void printHelpAndExit(int pnum){
  for( int h=0;h<pnum;++h ){
    printf("%s: %s\n", parameters[h].name, parameters[h].explanation );
  }
  exit(1);
}

void parseArgs( int argc, char ** argv ){
  
  int pnum = 0;
  while(parameters[pnum].name != 0) ++pnum;

  for( int i=1;i<argc;++i ){
    int pmatch = -1;
    for( int x = 0; x<pnum; ++x ){
      if( strcmp( parameters[x].name, argv[i] ) == 0 ){
	pmatch = x;
	if( strcmp( parameters[x].name, "-h" ) == 0 ){
	  printHelpAndExit(pnum);
	} 
	if( i+1>=argc ){
	  fprintf(stderr,"ERROR: No value provided for parameter %s.\n",argv[i]);
	  exit(1);
	}
	if( parameters[x].type == NUMBER )
	  *(parameters[x].value) = strtod(argv[++i],NULL);
	else
	  *(parameters[x].text) = strdup(argv[++i]);
      }
    }
    if(pmatch < 0){
      printHelpAndExit(pnum);
    }
  }
}

/*
 * this expression would take an mcfold html output and generate a .db file devoid of coaxial energies suitable for input to db2cm
 * sed -n '1,/Explored/{d};/BP/{q};/^>/,${p}' | cut -f 1,2 -d ' '
 */
int main( int argc, char **argv ){

  parseArgs( argc, argv );

  char dbfilename[LONGFILENAME];
  char seq_title[LONGSEQTITLE];
  char sequence[LONGEST_SEQ + 1];


  if( silent_mode == 0 ){
    printf("Magic Batman cache\n");
  }

	int cmd_line_fname_given = strcmp( inputfilename, "" );
  while( cmd_line_fname_given || 1 == scanf("%s",dbfilename) ){    

	if( cmd_line_fname_given )
	  strcpy( dbfilename, inputfilename );
    char db[LONGEST_SEQ];
    float score;
    
    FILE *DBF = fopen( dbfilename, "r" );
    if( DBF == NULL ){
      sleep(1);
      DBF = fopen( dbfilename, "r" );
      if( DBF == NULL ){
	fprintf(stderr,"ERROR: Unable to open file <%s>\n", dbfilename );
	continue;//return 1;
      }
    }

    fscanf(DBF, "%s", seq_title); 
    if( strcmp(seq_title,"Explored") == 0 )
      fscanf(DBF, "%s", seq_title);
    fscanf(DBF, "%s", sequence );




    //we read the dot-brackets from the file twice !
    // 1- in the first pass, we don't care about the structures themselves, just it's energy. 
    //    we determine MFE and partition function (sum of exponentials of energies)
    // 2- In the second pass, we parse each dot-bracket and 
    //    (a) add their base pairs to the basepairs matrix with correct weight (boltzmann, 1, other fixed)
    //    (b) add their base pairs to the batman signature with correct weight (boltzmann, 1, other fixed)
    
    long int START_OF_STRUCTURES = ftell( DBF );
    if(-1==START_OF_STRUCTURES){
      fprintf(stderr,"ERROR:bad file (%s).\n", dbfilename);
      fclose(DBF);
      continue;
    }

    double partition = 0.0;
    double mfe = 0.0;
    int is_mfe_set = 0;
    double maxBoltzmannProb = 0;

    int dbnum=0;
    while( 2 == fscanf( DBF, "%s %f", db, &score ) ){
	
      char c; 
      while( ' ' <= (c=fgetc(DBF) ) ) 
	;
      fputc(c,DBF); //eat any junk after energy value from mcff output

      if( ! is_mfe_set ){
	is_mfe_set = 1;
	mfe = score;
      }

      ++dbnum;

      //enforce explore parameter
      if( ((mfe - score)/mfe) > explore/1000.0 ){
	--dbnum;
	break;
      }
      //enforce xdb parameter
      if( xdb > 0 && dbnum > xdb ){
	--dbnum;
	break;
      }

      if( use_boltzmann ){
	double boltzmannProb = exp( (-1.0*score)/ RT );
	if( !isnormal( boltzmannProb) ){
	  perror( "Energy values too large for calculation of Boltzmann partition. Use -b 0." );
	  exit(0);
	}
	maxBoltzmannProb = max( maxBoltzmannProb, boltzmannProb );
	partition += boltzmannProb;
      } else {
	partition += 1;
      }
    }


    switch( (int)output_strNum ){
    case 0: break;
    case 1:
      fprintf( stdout, "%d\n", dbnum );
      break;
    case 2:
      fprintf( stderr, "%d\n", dbnum );
      break;
    default:
      fprintf(stderr, "WARNING: BAD VALUE (%f) PASSED IN FOR PARAMETER -o.\nShould be in { 0 1 2 }\n", output_strNum );
    }


    //initialise the batman signature
    double bat[ROWS][strlen(sequence)];
    for( int i=0;i<ROWS;++i)
      for( int j=0;j<strlen(sequence);++j)
	bat[i][j]=0.0;



    //initialise the basepairs matrix
    square *basepairs = square_alloc( strlen(sequence) );
    square_zero( basepairs );

    //second pass of reading file: (a) reset file, (b) parse
    if( -1 == fseek( DBF, START_OF_STRUCTURES, SEEK_SET ) ){
      fprintf(stderr,"ERROR: There was an error reading from file %s.\n", dbfilename );
      continue;
    }
    clearerr(DBF);

    dbnum=0;
    while( 2 == fscanf( DBF, "%s %f", db, &score ) ){
     
      char c; 
      while( ' ' <= (c=fgetc(DBF) ) ) 
	; 
      fputc(c,DBF); //eat any junk after energy value from mcff output

      ++dbnum; //one more dot bracket read.

      //enforce explore command line parameter
      if( ((mfe - score)/mfe) > explore/1000.0 ){
	--dbnum;
	break;
      }

      //enforce xdb command line parameter
      if( xdb>0 && dbnum > xdb ){
	--dbnum;
	break;
      }


      //the weight of this sequence in the contact map
      double db_weight = 1.0;
      if( use_boltzmann ){
	double boltzmannProb = exp( (-1.0*score)/ RT) / partition;
	if( !isnormal( boltzmannProb) ){
	  perror( "Energy values too large for calculation of Boltzmann partition. Use -b 0." );
	  exit(0);
	}
	db_weight = ((maxBoltzmannProb/partition) * (1.0 - (weight/1000.0))) + (boltzmannProb * (weight/1000.0));
      } 
      //force the use of an fixed value for w
      if( explicit_w >= 0 ){
	db_weight = explicit_w;
      }

      int *dbo = db2bat(db, basepairs, db_weight);

      for( int i=0; i<strlen(db); ++i ){
	bat[dbo[i]][i] += db_weight;
      }
      free(dbo);
      dbo=NULL;
    }
    fclose(DBF);


    //////////////////////////////////////// OUTPUT A .svg FILE WITH CONTACT MAP (basepairs)

	   
    char dotPlot_svgFileName[strlen(dbfilename)+strlen(dotPlot_fnameExt)+4];
    strcpy(dotPlot_svgFileName, dbfilename);
    strcat(dotPlot_svgFileName, dotPlot_fnameExt);
    strcat(dotPlot_svgFileName, ".svg");
	  
    printSVG( basepairs, //basepairs matrix to draw 
	      sequence,  //input RNA sequence 
	      use_boltzmann!=0?1.0:dbnum, //scaling factor for values in basepairs matrix (if boltzmann then already scaled)
	      dotPlot_sensitivity, //non-linearity of gray-ness plotting function (exponential)
	      dotPlot_rotation, //rotation of the image in degrees (it's a tiny bit complex because scalling must by applied to box the matrix diagonal)
	      //could be made better by accounting for arcs that do not fit the box. but it is not clear how to do this because the arcs
	      //are bezier curves.
	      dotPlot_zoom, //adjust image zoom
	      dotPlot_size, //adjust image size
	      dotPlot_threshold, //do not plot basepairs whose probabilities are this much or lower
	      dotPlot_svgFileName, //where to store the svg file
	      use_boltzmann,//isboltzmann
	      dbnum,//xdb
	      dbfilename, //basename
	      strlen(dotPlot_annotation) >0 ? dotPlot_annotation : NULL,
	      dotPlot_plot_arcs,
	      dotPlot_plot_dotPlot,
	      dotPlot_plot_seqOnDiagonal,
	      dotPlot_plot_seqOnDotPlot,
	      dotPlot_plot_legend
	      );
    

    //////////////////////////////////////// CONVERT THE RAW BATMAN SIGNATURE TO A VERSION NORMALIZED 
    //////////////////////////////////////// USING THE STATISTICS SPECIFIED ON ENTRY


    //if user chooses -max -sum -median -mean -stddev then set values of bat[c(0,1)][] to reflect this
    if( mode < 100 ){

      int len = basepairs->side;
      //forward pairings


      if( mode < 100 ){
	double *data = (double*)malloc(sizeof(double)*len);
	for( int i = 0; i < len; ++i ){
	  for( int j = i; j < len; ++j ){
	    data[j] = (double)*(square_idx(basepairs,i,j));
	  }
	  int n = len-i;
	  stats m = get_stats(data+i,n);
	  switch((int)mode){
	  case MAX : bat[0][i] = m.max; break;
	  case SUM : bat[0][i] = m.sum; break;
	  case Q1 : bat[0][i] = m.sum; break;	    
	  case MEDIAN : bat[0][i] = m.median; break;
	  case MEAN :   bat[0][i] = m.mean; break;
	  case STDDEV:  bat[0][i] = m.sd; break;
	  case ADEV:  bat[0][i] = m.adev; break;
	  case MAD:  bat[0][i] = m.mad; break;
	  default:
	    fprintf(stderr,"ERROR: indeterminate mode in %s at line %d\n", __FILE__, __LINE__);
	    exit(1);
	  }

	  if( mode2 >= 0 )
	    switch((int)mode2){
	    case MAX : bat[7][i] = m.max; break;
	    case SUM : bat[7][i] = m.sum; break;
	    case Q1 : bat[7][i] = m.sum; break;	    
	    case MEDIAN : bat[7][i] = m.median; break;
	    case MEAN :   bat[7][i] = m.mean; break;
	    case STDDEV:  bat[7][i] = m.sd; break;
	    case ADEV:  bat[7][i] = m.adev; break;
	    case MAD:  bat[7][i] = m.mad; break;
	    default:
	      fprintf(stderr,"ERROR: indeterminate mode in %s at line %d\n", __FILE__, __LINE__);
	      exit(1);
	    }
	}
	
	for( int i = 0; i < len; ++i ){
	  double reverse = 0;
	  for( int h = 0; h < i; ++h ){
	    data[h] = *(square_idx(basepairs,h,i));
	  }
	  int n = i;
	  stats m = get_stats(data,n);
	  switch((int)mode){
	  case MAX : bat[1][i] = m.max; break;
	  case SUM : bat[1][i] = m.sum; break;
	  case Q1 : bat[1][i] = m.sum; break;	    
	  case MEDIAN : bat[1][i] = m.median; break;
	  case MEAN :   bat[1][i] = m.mean; break;
	  case STDDEV:  bat[1][i] = m.sd; break;
	  case ADEV:  bat[1][i] = m.adev; break;
	  case MAD:  bat[1][i] = m.mad; break;
	  default:
	    fprintf(stderr,"ERROR: indeterminate mode in %s at line %d\n", __FILE__, __LINE__);
	    exit(1);
	  }
	  if( mode2 >= 0 ) 
	    switch((int)mode2){
	    case MAX : bat[8][i] = m.max; break;
	    case SUM : bat[8][i] = m.sum; break;
	    case Q1 : bat[8][i] = m.sum; break;	    
	    case MEDIAN : bat[8][i] = m.median; break;
	    case MEAN :   bat[8][i] = m.mean; break;
	    case STDDEV:  bat[8][i] = m.sd; break;
	    case ADEV:  bat[8][i] = m.adev; break;
	    case MAD:  bat[8][i] = m.mad; break;
	    default:
	      fprintf(stderr,"ERROR: indeterminate mode in %s at line %d\n", __FILE__, __LINE__);
	      exit(1);
	    }
	}	  	
	free(data);
      } else {
	fprintf(stderr,"ERROR: indeterminate mode in %s at line %d\n", __FILE__, __LINE__);
	exit(1);
      }
    } else { //if user chooses mode in 100.. , then dump base pairs to stdout
      int len = basepairs->side;
      //forward pairings

      const double min1 = -1;
      for( int i = 0; i < len; ++i ){
	//padd lines with negative values
	int j;
	for( j=0; j<i; ++j ){
	  fprintf(stdout,"%f ", min1);
	}
	for(j = i; j < len; ++j ){
	  fprintf(stdout,"%f ", *(square_idx(basepairs,i,j)) );
	}
	fprintf(stdout,"\n");
      }

      for( int i = 0; i < len; ++i ){
	int h;
	for( h = 0; h < i; ++h ){
	  fprintf(stdout,"%f ", *(square_idx(basepairs,h,i)) );
	}
	for( h = i; h<len; ++h ){
	  fprintf(stdout,"%f ", min1 );
	}
	fprintf(stdout,"\n");
      }	    
    }
   

    //normalize column wise
    for( int i=0;i<strlen(db);++i ){
      double sum=0.0;
      for( int j=0; j<ROWS; ++j ){
	sum+=bat[j][i];
      }
      if(sum>0.0) 
	for( int j=0; j<ROWS; ++j ){
	  bat[j][i] /= sum;
	}      
    }
    
    if( silent_mode == 0 ){
      printf("%s\n",dbfilename);
      printf("%zi", strlen(sequence));
      for( int i=0; i<strlen(sequence);++i ){
	printf("(%f,%f,%f,%f,%f,%f,%f,%f,%f)",
	       bat[0][i],bat[1][i],bat[2][i],bat[3][i],bat[4][i],bat[5][i],bat[6][i],bat[7][i],bat[8][i] );	
      }
      printf("\n");
    }
    //clean up ram
    basepairs = square_free( basepairs );

	  if(cmd_line_fname_given)
		  break;
  }//while more filenames
}
