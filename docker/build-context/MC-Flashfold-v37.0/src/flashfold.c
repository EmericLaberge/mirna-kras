/*
 * flashfold.c
 * Paul Dallaire
 *
 * f37d: (2018-03-23) corrected stupid bug that affected -zd -zzd and -sd
 *                    The problem was caused in main where the reactivity mask lenght was allocated prior to the duplex modes definition (command line application) 
 *                    resulting in reactivity mask being defined on the lenght of the sequence passed via -s only not taking in to account that the actual length of the sequence
 *                    once duplex sequence is combined with -s (plus the intervening AA)
 *                    Also corrected a minor memory leak.
 *
 * f37: (2018-01-15) add option -nri (reactivity mask as int)
 *
 * f36: (2018-01-09) add reactivity mask
 *                   adapted code for windows compatibility (some of it is sorta hacked and could be made better)
 *                   removed variable lenght arrays because unable to turn that on in visual studio
 *
 * f35: (2016-05-25) rare bug Mathieu causing repeated search with identical threshold when structure pool empty.
 *
 * f34: (2015-05-14) debug rare empty output on some sequences by introducing MAX_EMPTY_THRESHOLD.
 *                   elevated by 1 the verbosity level necessary to print the location of the tables used (was 1 is now 2)
 *                   Gabriel C-Parent: changed file namimg scheme to flashfold.c
 *
 * f33: (2013-08-16) add parameter -e
 *
 * f32: (2013-07-11) changed lookup for directory tables (1-parameter,2-environment variable "MCFTABLES",3-./tables, 4-~/MC-Flashfold/tables, 5-/usr/local/MC-Flashfold/tables)
 *
 * f31: (2013-07-11) bugfix -> Corrected improper logic in masks allowed free nucleotides in terminal free loops and at MB boundaries causing some masks to fail.
 *
 * f30: bugfix -> Correct verifications of unpaired nucleotides in loops for calculations of RO[] values.
 *
 * f29:
 *    : bugfix -> detect situations where -ft X is used but the number of possible solutions is < X, print a warning, replace X with returned solutions number
 *    : bugfix -> some masks mislead to a wrong (slightly lower) threshold: fix: when -t=0 and -ft is unspecified, used -ft 1
 *
 * f28: thiner output for -gs as -gst
 *    : output parse in unsorted outputs using -show-parse
 *
 * f27: started to add bulges or arbitrary lengths in graph connectivity calculation
 *      duplex modes cleanup
 *
 * f26: augment the rule set of reachable dot brackets
 *
 * f25: add parameter -et (explicit rules to match two dot brackets as reachable in the graph)
 *      output masks that match any member of a component is marked in the graph summary.
 *      BUG FIX: -om masks did not match all occurences since f21 when -alt used
 *      full balanced masks can contain non-balanced information making them more friendly to the user.
 *
 * f22: add canonical base pairs mask
 *
 * f21: Bug fixes in backtrack (no change on output and no real speed gain)
 *      time measures for benchmarcking (added -times)
 *      new command line parameter to specify the location of the tables directory (-tables dir)
 *      Bug fix in forward() causing possible memory corruption. (no observed change in behaviour)
 *      Added command line parameter -alt to output dot brackets that are different when base pair is not GC, UA
 *      2013-03-19/20
 *
 * f20: Iterative refinement of threshold estimate when heap is too small (changed -max to -ft)
 *      Removed -ltl replaced with -nltl making long terminal loops default to TRUE.
 *      Some bugfixes
 *
 * f19: Unbroken. New stable.
 *      MIND can now be set to 1. (but all 1_1_3 NCM are +Inf).
 *
 * f18: Stable but broken.
 *      fixed bug in graph calculation where radius was not updated properly upon bucket fusion.
 *      fixed masks bug mentionned in f17. -um and -m now give structures for calculated mfe.
 *      The error was due to lack of checking masks for permission of unpaired nucleotides
 *      in NCMs in forward().
 *      other bug fixes.
 *      augmented explanatory comments.
 *      2013-03-13
 *
 * f17: add balanced full masks (-m) and unbalanced masks (-um)
 *      there is a difference between reported MFE and best solution found when using some masks on some sequences...
 *      Quick fix: use -max with larger threshold. The difference apears due to an overshoot of MB[i,j] and not to a backtrack problem.
 *      2013-03-11
 *
 * f16: add long terminal loops (-ltl)
 *
 * f15: add parameter -STEM_OF_ONE
 * f15: use SIMD for comparing db as short ints in graph components extraction (no apparent benefit, shit.)
 *
 * f14: use bucket list metric data structure for graph components extraction. (some speed gain in graph calculation. unclear)
 *
 * f13: proper parsing of command parameters
 *
 * f12: started working on comments. More cleanup to be done.
 * f12: bugfix in backtrack producing twice same structure once through the suite BRANCH->R->L and once through BRANCH->L
 * f12: add a simple duplex mode (-sd, -zd, -zzd) for use with microRNAs such as in the stable mariages project 2013-03-05
 *
 * f11: add a simple output balanced mask facility (-om)
 *
 * f10: add arguments -g and -gs. ie: calculate the graph of reachable components
 *
 * f9: add argument -max (a heap to count the proper treshold that will output the requested number of structures)
 * f9: memory reduction: use only half size matrices and arrays whenever possible (this is painfull)
 *
 * f8: work on memory reduction.
 * 1- StackSize cut from N^2 to N
 * 2- doubles replaced by floats (this does introduce numerical rounding error,
 *    so the value of iota is changed and threshold is modified accordingly in the backtrack routine.
 *
 * f7: correct output stable clean (but see f12, backtrack not corrected)
 *
 * boring compile line:
 *     gcc  -o mcff flashfold.c -lm -std=c99 -O3
 * cool compile line:
 *     gcc -o mcff flashfold.c -lm -std=c99 -O3 -march=native -mfpmath=sse -msse2 -msse4.1 -D___SIMD_COMPILATION__=1
 *     strip -s mcff (stripping shrinks a bit the size of executable)
 * on majsrv4:
 *     gcc -o mcff-majsrv4 flashfold.c -lm -std=c99 -O3 -mfpmath=sse -msse2 -DNOTIMESTAMP=1
 *
 * if you want square brackets for shapes define SQUARE_BRACKET_SHAPES to some value: -DSQUARE_BRACKET_SHAPES=1
 */


#include <assert.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>

#include <string.h>
 //#include <sys/times.h>
#include <time.h>
//#include <sys/resource.h>
#include <ctype.h>

//#define WindowsCompilation 1
#ifndef WindowsCompilation

#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>

#else


typedef size_t ssize_t;

#include <io.h>
//#include <Shlwapi.h> //for PathIsDirectory

struct rusage {
	int a;
};

int getrusage(void* a, void* b) {
	return 42;
}


#define S_IRGRP 0
#define S_IROTH 0 

#endif // !WindowsCompilation




#ifdef ___SIMD_COMPILATION__
#include <emmintrin.h> // - X86 SSE2, posix_memalign
#endif

//return code for calls to exit()
#define NORMAL_EXIT 1
#define ABORT_EXIT 0

// to keep track of the version
#ifndef SOFTWARE_VERSION
#define SOFTWARE_VERSION "37"
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  BENCHMARKING

int KEEP_TIME = 0; //if set then output times and other statistics of various passes. (set in main via the command line option -time)

long long int BACKTRACK_NODES = 0; //Nodes in the backtracking stack (number of times around the while IN THE LAST BACKTRACK INVOCATION)
int SOLUTIONS = 0; //Number of solutions found during the last backtrack invocation
int COMPONENTS = 0; //Number of graph components identified
int INPUT_LENGTH = 0; //N

//number of time points
#define TIME_POINTS 10

//the following spots are measured
enum {
	AT_BEGIN = 0, BEFORE_FORWARD, AFTER_FORWARD, BEFORE_BACKTRACK, AFTER_BACKTRACK, BEFORE_LAST_BACKTRACK,
	AFTER_LAST_BACKTRACK, BEFORE_OUTPUT, AFTER_OUTPUT, AT_END
};

struct rusage ticks[TIME_POINTS];

//get the number of seconds between the two time points
//a happens before b and both values are <TIME_POINTS
double ds(int b, int a) {
#ifndef WindowsCompilation
	struct timeval ta = ticks[a].ru_utime;
	struct timeval tb = ticks[b].ru_utime;
	return ((double)(tb.tv_sec - ta.tv_sec))
		+ ((double)tb.tv_usec - ta.tv_usec) / 1e6;
#else
	return 0;
#endif
}

//print collected stats
//if you CHANGE this function, then you must UPDATE report_times_header() (please)
void report_times(void) {
	double f = ds(AFTER_FORWARD, BEFORE_FORWARD);
	double b = ds(AFTER_BACKTRACK, BEFORE_BACKTRACK);
	double b2 = ds(AFTER_LAST_BACKTRACK, BEFORE_LAST_BACKTRACK);
	double out = ds(AFTER_OUTPUT, BEFORE_OUTPUT);
	double all = ds(AT_END, AT_BEGIN);

	fprintf(stderr, "%i %i %i %lli %5.3f %5.3f %5.3f %5.3f %5.3f\n",
		INPUT_LENGTH, SOLUTIONS, COMPONENTS, BACKTRACK_NODES,
		f, b, b2, out, all);
}

//output header line for stats reported by report_times()
void report_times_header(void) {
	fprintf(stderr, "N solutions components nodes forward backtrack backtrack2 output global\n");
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Varia

//all heap memory alligned on 16 bytes boundaries if SIMD else, malloc
void *alloc_aligned(ssize_t sz) {
#ifdef ___SIMD_COMPILATION__
	void*ptr;
	if (0 != posix_memalign(&ptr, 16, sz)) {
		fprintf(stderr, "ERROR: Out of memory.\n");
		exit(ABORT_EXIT);
	}
	return ptr;
#else
	void * ptr = malloc(sz);
	if (NULL == ptr) {
		fprintf(stderr, "ERROR: Out of memory.\n");
		exit(ABORT_EXIT);
	}
	return ptr;
#endif
}


//DEPRECATED FROM standard unix because of security risks
char *strdup(const char *str)
{
	ssize_t n = strlen(str) + 1;
	char *dup = alloc_aligned(n);
	if (dup)
	{
		strcpy(dup, str);
	}
	return dup;
}

//analoguous to strdup for any type
void *memdup(void* m, size_t N) {
	void *dup = alloc_aligned(N);
	memcpy(dup, m, N);
	return dup;
}

//verifies that a pointer is not null and exits otherwise.
void check_out_of_memory(void*ptr, const char *file, const int line_num) {
	if (NULL == ptr) {
		fprintf(stderr, "ERROR: OUT OF MEMORY in file %s at line %i.\n", file, line_num);
		exit(ABORT_EXIT);
	}
}

char *mergeToNewString(const char *a, const char *b) {
	char *r = (char*)alloc_aligned((1 + strlen(a) + strlen(b)) * sizeof(char));
	check_out_of_memory(r, __FILE__, __LINE__);
	strcpy(r, a);
	strcpy(r + strlen(a), b);
	return(r);
}


/*
	utility to split input parameter in the form "1.2 3.2 -0.2 0 0 0",
	convert the values to double and dump them in a preallocated vector
	Used in main to parse an argument
*/
int splitToFloat(char * values, float * reactivity, const int expectedEntries)
{
	int nt = 0;
	int vlen = strlen(values);
	for (int i = 0; i < vlen; ++i)
	{
		while (values[i] <= ' ' && i < vlen)
			++i;
		if (i < vlen) {
			if (nt < expectedEntries) {
				reactivity[nt] = atof(values + i);
			}
			++nt;
		}
	}
	return nt;
}
/////////////////////////////////////////////////////////////////////////////////////////////// HEAP STUFF

//heap[0] is max number of items in heap
//heap[1] is number of current items in heap
//heap[2..(2+(n-1))] heap data (where n is heap[0]: capacity)
//the heap is inverted so as to easily remove a smallest value

float *allocate_heap(size_t capacity) {
	float *heap = (float *)alloc_aligned((2 + capacity) * sizeof(float));
	if (heap != NULL) {
		heap[0] = capacity;
		heap[1] = 0;
	}
	return heap;
}


void sift_down_heap(float *heap, size_t i) { //sift node i down
	const size_t n = heap[1];
	heap += 1;//from now on heap is 1 based and can be indexed from 1..n inclusive (don't write in heap[0] though)
	size_t k = i, j;
	do {
		j = k;
		if (2 * j <= n && heap[2 * j] > heap[k]) k = 2 * j;
		if (2 * j < n && heap[(2 * j) + 1] > heap[k]) k = (2 * j) + 1;
		float t = heap[j]; heap[j] = heap[k]; heap[k] = t;
	} while (j != k);
}

void percolate_heap(float *heap, size_t i) {
	heap += 1;//from now on heap is 1 based and can be indexed from 1..n inclusive (don't write in heap[0] though)
	size_t k = i, j;
	do {
		j = k;
		if (j > 1 && heap[j / 2] < heap[k]) k = j / 2;
		float t = heap[j]; heap[j] = heap[k]; heap[k] = t;
	} while (j != k);
}


//add v to the heap if either the heap is not full or v is smaller than the least value
void add_to_heap(float *heap, float v) {
	assert(heap[0]);
	if (heap[1] < heap[0]) { //there is room, insert
		//not full
		size_t idx = ((size_t)++heap[1]);
		heap[idx + 1] = v;
		percolate_heap(heap, idx);
	}
	else if (heap[2] > v) {
		//v smaller than least value in heap
		heap[2] = v;
		sift_down_heap(heap, 1);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////// ARRAYS FACILITIES
/*
 * X dimensions matrix class wrapped arround R multidim arrays
 * X <= 5
 */
struct array_t {
	float *data;
	int DIM[5]; //ROWS, COLS, OTHERS.
	int dimensions;
};


//general array accessor functions (use ONLY with normal multidimensional arrays allocated with aXalloc (see lower)
#define a1d( M, c) (M.data[(c)])
#define a2d( M, r, c ) (M.data[((c)*M.DIM[0])+(r)])
#define a3d( M, a, b, c ) (M.data[ ((c)*(M.DIM[0]*M.DIM[1])) +((b)*M.DIM[0]) + (a)])
#define a4d( M, a, b, c, d ) (M.data[((d)*(M.DIM[2]*M.DIM[0]*M.DIM[1]))+((c)*(M.DIM[0]*M.DIM[1]))+((b)*M.DIM[0])+(a)])
#define a5d( M, a, b, c, d, e ) (M.data[ ((e)*(M.DIM[3]*M.DIM[2]*M.DIM[0]*M.DIM[1])) + ((d)*(M.DIM[2]*M.DIM[0]*M.DIM[1])) +((c)*(M.DIM[0]*M.DIM[1]))+((b)*M.DIM[0])+(a)])


//specialized accessor used for square 2D matrices allocated using square2Halloc
#define a2dh( M, r, c ) ( M.data[ ( (M.DIM[0]) * ((r)) )-(((r)*((r)+(1)))>>(1))+(c)] )

//specialized accessor used for the array F (allocated using square3Halloc)
#define a3dh( M, a, b, c ) ( M.data[ ( M.DIM[2] * (( (M.DIM[0]) * ((a)) )-(((a)*((a)+(1)))>>(1))+(b) ) ) + c])


//use with normal multidimensional arrays only: crash hard if something went wrong
void aXalloc(struct array_t *A, int dimensions, int N1, int N2, int N3, int N4, int N5) {
	A->dimensions = dimensions;
	A->DIM[0] = N1;
	A->DIM[1] = N2;
	A->DIM[2] = N3;
	A->DIM[3] = N4;
	A->DIM[4] = N5;
	int howBig = 1;
	for (int i = 0; i < dimensions; ++i) howBig *= A->DIM[i];
	A->data = (float*)alloc_aligned(sizeof(float) * howBig);
	if (A->data == NULL) {
		printf("OUT OF MEMORY at line %i in %s\n", __LINE__, __FILE__);
		exit(ABORT_EXIT);
	}
}

//main diag + upper right square matrix
void square2Halloc(struct array_t *A, int N) {
	A->dimensions = 2;
	A->DIM[0] = N;
	int howBig = ((N + 1)*N) / 2;
	A->data = (float*)alloc_aligned(sizeof(float) * howBig);
	if (A->data == NULL) {
		printf("OUT OF MEMORY at line %i in %s\n", __LINE__, __FILE__);
		exit(ABORT_EXIT);
	}
}

//used only for F (a square matrix whose elements are vectors of NCMs
void square3Halloc(struct array_t *A, int N, int K) {
	A->dimensions = 3;
	A->DIM[0] = N;
	A->DIM[1] = N;
	A->DIM[2] = K;
	int howBig = K * (((N + 1)*N) / 2);//(N+1)*N is always even
	A->data = (float*)alloc_aligned(sizeof(float) * howBig);
	if (A->data == NULL) {
		printf("OUT OF MEMORY at line %i in %s\n", __LINE__, __FILE__);
		exit(ABORT_EXIT);
	}
}

//use with normal multidimensional arrays only
void initArray(struct array_t M, float val) {
	size_t len = 1;
	for (int d = 0; d < M.dimensions; ++d)
		len *= M.DIM[d];
	for (size_t i = 0; i < len; ++i) {
		*(M.data + i) = val;
	}
}

//main diag + upper right square matrix
void initHArray(struct array_t M, float val) {
	int N = M.DIM[0];
	size_t len = ((N + 1)*N) / 2;
	for (size_t i = 0; i < len; ++i) {
		*(M.data + i) = val;
	}
}

//used only for F (a square matrix whose elements are vectors of NCMs but whose only main diag + upper right are allocated.)
void init3HArray(struct array_t M, float val) {
	int N = M.DIM[0];
	size_t len = M.DIM[2] * (((N + 1)*N) / 2);
	for (size_t i = 0; i < len; ++i) {
		*(M.data + i) = val;
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////// FLASHFOLD DEFINITIONS

static const float LARGEFLOAT = FLT_MAX - 1000; //leaway for overflow
static const float iota = (float)1e-9;//1e-10 for double

const int MIND = 2; //Minimum unpaired nucleotides in a terminal loop. (can set MIND=1 if you want but energies in MCFold original database for NCMs of type (.) (id=1 or id=0 with zero based offset) are all set to +Inf.
const int MIN_MB = 3;//8; //(2*MIND) + 4; (min space required to start a multi-branch), changed to 3 to allow for terminal long loops
//const int MIN_RO = 10;//MIND + 4 + 4; (min space required to start a right open (long) internal loop)
const int MAX_NCM_TERMINAL_LOOP = 6; //max number of nucleotides in a terminal loop described by an NCM (including closing bp)
const int NCMs = 21; //max NCM identifier (or NCMID)
const int FIRST_TWO_STRANDS_NCM = 4; //Minimum NCMID of 2 strands. Smaller values are terminal loop NCMs (or one strand NCMs).

int ALLOW_STEMS_OF_ONE = 0; //hairpin loops consisting of only one terminal 1_1_X NCM are normally dissallowed ( 0 ).
int ALLOW_LONG_TERMINAL_LOOP = 0; //allow long unpaired nucleotides loop to cap stems


/*
 * DUPLEX MODE
 * Introduced in the code (of version 12) to satisfy the simple hybridization modes of one oligonucleotide (or micoRNA)
 * to anoter strand of RNA (or mRNA in the vicinity of a MRE).
 * CONSTRAINTS:
 *  The first nucleotide of the second strand (or microRNA) pairs to the last nucleotide of the first strand (mRNA).
 *  There are no terminal loops in either strand.
 * MODEL:
 *  The user supplies two strands. One (the mRNA) via the normal parameter (-s) and the other via a second parameter (-sd, -zd or -zzd).
 *  If -sd is used then the mode is SIMPLE. If -zd is used then the mode is NCM (or zip).
 *  If -zzd is used then the mode is NCM4 (or ZIPZIP).
 *  The strings are concatenated with a AA in between in main().
 *  The value of DUPLEX_HAIRPIN_START is set to the length of the first sequence.
 *  The value of DUPLEX_HAIRPIN_LENGTH is set to the length of AA. (ie: 2)
 *  The value of DUPLEX_MODE is set to either SIMPLE, NCM or NCM4.
 *  In the forward pass (inside), only one loop (1_1_2 XAAY; where X is the last nt of the first strand (mRNA) and Y is the first
 *  nt of the second strand (microRNA) is initialized to the value 0.0. The other loops are inited to LARGE_FLOAT.
 *  Updates to the tables are limited or modified in the various cases according to the value of DUPLEX_MODE.
 *  Backtracking is modified accordingly.
 *  Output is only the score by default. If verbose is set to non-zero then the structure is output also.
 * MODES:
 *  SIMPLE -> allow long loops in either arm anywhere normally allowed by the folding model.
 *  NCM    -> allow any NCM and possibly a single long loop at the start of the 5' strand (mRNA) but nowhere else.
 *  NCM4   -> allow the use of only 2_2_2 and nothing else; the lengths of both sequences must match exactly.
 * NOTES:
 *  Finding the proper treshold is troublesome for backtracking.
 */

int DUPLEX_MODE = 0; //this flag set if -sd -zd -zzd used at command line
int DUPLEX_HAIRPIN_START = -1; //initialized in main
int DUPLEX_HAIRPIN_LENGTH = -1; //initialized in main

enum DUPLEX_MODE_VALUE {
	NONE = 0,
	SIMPLE = 1,
	NCM = 2,
	NCM4 = 3
};

/////////////////////////////////////////////////////////////////////////////////////////////// UTILITIES
/*
 * Generate an index value usable to access the NCMs energy values from the table energies.
 * seqAsInt is the sequence of nucleotides coded as integers (A=1, C=2, G=3, U=4)
 * i..k is the interval for the 5' strand of the NCM
 * j..l is the interval for the 3' strand of the NCM
 * if there is only one strand in the NCM, then set j and l to -1
 * The value computed is subtracted by 20 because the energies table does not contain values for the first 20 entries (dinucleotides)
 * The value is further subtracted by 1 because the formula computes indexes usable in R where tables and vectors start with value 1.
 */
int seq2idx(const int *seqAsInt, const int i, int k, const int j, int l) {
	int val = 0;
	int posInWord = 1;
	if (l > 0) while (l >= j) { val += seqAsInt[l] * (float)trunc((float)pow(4, posInWord - 1)); ++posInWord; --l; }
	while (k >= i) { val += seqAsInt[k] * (float)trunc((float)pow(4, posInWord - 1)); ++posInWord; --k; }
	if (val < 21) {
		printf("ERROR in seq2idx ikjl(%i,%i,%i,%i)\n", i, k, j, l);
		volatile int *N = 0;
		*N = 1;
	}
	return val - 20 - 1;
}

/*
 * convert a RNA as character string to integer string where A->1, C->2, G->3, U->4, T->4
 * crash program if another character found in stream.
 */
int *toSeqAsInt(const char *seq) {
	if (!seq) {
		printf("ERROR: no sequence buffer passed at line %i in %s\n.", __LINE__, __FILE__);
		exit(ABORT_EXIT);
	}
	size_t seqLen = strlen(seq);
	if (seqLen == 0) {
		return 0;
	}
	int *seqAsInt = (int*)alloc_aligned(sizeof(int) * (seqLen + 1));
	if (seqAsInt == NULL) {
		printf("ERROR: OUT OF MEMORY at line %i in %s\n", __LINE__, __FILE__);
		exit(ABORT_EXIT);
	}

	for (int i = 0; i < seqLen; ++i) {
		switch (seq[i]) {
		case 'A':
		case 'a':
			seqAsInt[i] = 1;
			break;
		case 'C':
		case 'c':
			seqAsInt[i] = 2;
			break;
		case 'G':
		case 'g':
			seqAsInt[i] = 3;
			break;
		case 'U':
		case 'u':
		case 'T':
		case 't':
			seqAsInt[i] = 4;
			break;
		default:
			printf("ERROR: bad char %c in sequence\n", seq[i]);
			exit(ABORT_EXIT);
		} //switch
	} //for i
	return seqAsInt;
}//toSeqAsInt()


#ifdef SQUARE_BRACKET_SHAPES
#define OPEN_SHAPE_CHAR '['
#define CLOSE_SHAPE_CHAR ']'
#else
#define OPEN_SHAPE_CHAR '('
#define CLOSE_SHAPE_CHAR ')'
#endif

char *db2shape(char *dbwithdots) {
	//char db[strlen(dbwithdots)];
	char *db = alloc_aligned(sizeof(char) * strlen(dbwithdots));
	int dblen = 0;
	for (int i = 0; i < strlen(dbwithdots); ++i) {
		if (dbwithdots[i] != '.')
		{
			char c = dbwithdots[i];
			c = c == '>' ? ')' : (c == '<' ? '(' : c);
#ifdef SQUARE_BRACKET_SHAPES
			c = c == ')' ? ']' : (c == '(' ? '[' : c);
#endif 
			db[dblen++] = c;
		}
	}

	//compute the idb (integer dot bracket) of the db
	//int buddies[dblen];
	//int stack[dblen];
	int * buddies = alloc_aligned(sizeof(int) * dblen);
	int *stack = alloc_aligned(sizeof(int) * dblen);
	int sp = -1;
	int maxsp = -1; //needed to know if a db is only filled with '.'
	for (int i = 0; i < dblen; ++i) {
		if (db[i] == OPEN_SHAPE_CHAR) {
			stack[++sp] = i;
		}
		else if (db[i] == CLOSE_SHAPE_CHAR) {
			buddies[i] = stack[sp];
			buddies[stack[sp]] = i;
			--sp;
		}
		maxsp = maxsp > sp ? maxsp : sp;
	}

	int slen = 1; //1 and not 0 because if db[0] != '.' then there is no room for the first (
	//1find length of shape
	for (int i = 1; i < dblen; ++i) {
		if ((buddies[i - 1] - buddies[i] != 1) | (db[i] != db[i - 1])) ++slen;
	}
	//2alloc ram and fill shape
	char *shape = (char*)alloc_aligned(sizeof(char)*(slen + 1));
	if (shape == NULL) {
		fprintf(stderr, "ERROR: OUT OF MEMORY in %s at line %i\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}

	int next = 0;
	if (maxsp > -1)shape[next++] = OPEN_SHAPE_CHAR;
	for (int i = 1; i < dblen; ++i) {
		if (((buddies[i - 1] - buddies[i] != 1) | (db[i] != db[i - 1])))
			shape[next++] = db[i];
	}
	shape[next] = 0;

	free(stack);
	free(buddies);
	free(db);

	return shape;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////// USER EXTRA REACTIVITY PER NUCLEOTIDE
/*
 * Energy values added to any pair involving specified nucleotides.
 * This vector is created as soon as sequence length is known so that it is always usable by any computation.
 * The values specified by the user are added to this vector at options parse.
*/
float * reactivity = NULL;


/*
 * output mask facility
 * outputs that respect the constraints specified in these masks are identified as such in some outputs.
 */
 /////////////////////////////////////// UMASKS DEFINES

 // global mask if set are in these variable
int * FULL_MASK = NULL;
int * UNBALANCED_MASK = NULL;


//set to 1 in main() if user specified -alt flag
//this allows proper output masks matching
//because there is a semantic change in the value of a () base pair if -alt is on
int USING_ALTERNATE_OUTPUT_DB = 0;

#define UMASK_FORWARD_PAIRED     1 //(((int)1))  // (
#define UMASK_BACK_PAIRED        (1<<1) //(((int)1)<<1)  // )
#define UMASK_PAIRED             (1<<2) //(((int)1)<<2)  // |
#define UMASK_UNPAIRED           (1<<3) //(((int)1)<<3)  // .

#define UMASK_NOT_C_PAIRED       (1<<4) //(((int)1)<<4) // - is not involved in a canonical pair
#define UMASK_FORWARD_C_PAIRED   (1<<5) //(((int)1)<<5)  // [ is forward canonical paired
#define UMASK_BACK_C_PAIRED      (1<<6) //(((int)1)<<6)  // ] is reverse canonical paired
#define UMASK_C_PAIRED           (1<<7) //(((int)1)<<7)  // + is involved in canonical pairing

#define UMASK_NOT_NC_PAIRED      (1<<8) //(((int)1)<<8)  // _ is not involved in a non-canonical pair
#define UMASK_FORWARD_NC_PAIRED  (1<<9) //(((int)1)<<9)  // < is forward non-canonicaly paired
#define UMASK_BACK_NC_PAIRED     (1<<10) //(((int)1)<<10)  // > is reverse non-canonicaly paired
#define UMASK_NC_PAIRED          (1<<11) //(((int)1)<<11)  // ! is involved in a non-canonical pair

#define UMASK_DONTCARE           (1<<12) // (((int)1)<<12)  // x

#define UMASK_NOT_REVERSE_PAIRED (1<<13) //q
#define UMASK_NOT_FORWARD_PAIRED (1<<14) //p

//following defines used in __upair() and __ufree macros
//O:Open bracket, C: Closing bracket, C: canonical, NC: Non-canonical

#define UMASK_NT_CAN_BE_FREE  ( UMASK_UNPAIRED | UMASK_NOT_C_PAIRED | UMASK_NOT_NC_PAIRED | UMASK_DONTCARE | UMASK_NOT_REVERSE_PAIRED | UMASK_NOT_FORWARD_PAIRED )

#define UMASK_O_C_CAN_PAIR   ( UMASK_DONTCARE | UMASK_PAIRED | UMASK_FORWARD_PAIRED | UMASK_NOT_NC_PAIRED | UMASK_C_PAIRED | UMASK_FORWARD_C_PAIRED | UMASK_NOT_REVERSE_PAIRED )

#define UMASK_C_C_CAN_PAIR   ( UMASK_DONTCARE | UMASK_PAIRED | UMASK_BACK_PAIRED | UMASK_NOT_NC_PAIRED | UMASK_C_PAIRED | UMASK_BACK_C_PAIRED | UMASK_NOT_FORWARD_PAIRED )

#define UMASK_O_NC_CAN_PAIR  ( UMASK_DONTCARE | UMASK_PAIRED | UMASK_FORWARD_PAIRED | UMASK_NOT_C_PAIRED | UMASK_NC_PAIRED | UMASK_FORWARD_NC_PAIRED | UMASK_NOT_REVERSE_PAIRED )

#define UMASK_C_NC_CAN_PAIR  ( UMASK_DONTCARE | UMASK_PAIRED | UMASK_BACK_PAIRED | UMASK_NOT_C_PAIRED | UMASK_BACK_NC_PAIRED | UMASK_NC_PAIRED | UMASK_NOT_FORWARD_PAIRED )

//following defines used for matching output masks only

#define UMASK_O_ANY_CAN_PAIR      ( UMASK_O_C_CAN_PAIR | UMASK_O_NC_CAN_PAIR )

#define UMASK_C_ANY_CAN_PAIR      ( UMASK_C_C_CAN_PAIR | UMASK_C_NC_CAN_PAIR )




#define OUTPUT_MASK_BALANCED     1
#define OUTPUT_MASK_UNBALANCED   2

struct output_mask_t {
	int * balanced_idb;
	int * mask_idb;
	char *identifier;
	int type;
};

struct output_mask_t *omasks = NULL;
int omasksSize = 0;
int omasksCapacity = 0;

//true if the passed dot bracket as char matches the dot bracket passed as vector of integers
int doesDBmatchOutputMask(const char *db, const int mi) {// mi is index of mask in omasks const int *mask_idb, const int maskType ){

	//starting with f26, output balanced masks are also unbalanced ones

	int dblen = (int)strlen(db);
	if (omasks[mi].type == OUTPUT_MASK_BALANCED) {
		//compute the idb (integer dot bracket) of the db
		//int buddies[dblen];
		//int stack[dblen];
		int *buddies = alloc_aligned(sizeof(int) * dblen);
		int *stack = alloc_aligned(sizeof(int) * dblen);

		int sp = -1;
		for (int i = 0; i < dblen; ++i) {
			if (db[i] == '(' || db[i] == '<') {
				stack[++sp] = i;
			}
			else if (db[i] == ')' || db[i] == '>') {
				buddies[i] = stack[sp];
				buddies[stack[sp]] = i;
				--sp;
			}
			else {//a .
				buddies[i] = i;
			}
		}

		//int*balanced_idb=(int*)mask_idb;

		for (int i = 0; i < dblen; ++i) {
			if (omasks[mi].balanced_idb[i] > -1 && buddies[i] != omasks[mi].balanced_idb[i])
				//missmatch
				free(stack);
				free(buddies);
			return 0;
		}
		free(stack);
		free(buddies);
		//return 1;
	} //else if( maskType == OUTPUT_MASK_UNBALANCED ){ //maskType is unbalanced
	int* mask = omasks[mi].mask_idb;
	for (int i = 0; i < dblen; ++i) {
		if (!USING_ALTERNATE_OUTPUT_DB) {
			if (
				!((db[i] == '.' && (mask[i] & UMASK_NT_CAN_BE_FREE))
					||
					(db[i] == '(' && (mask[i] & UMASK_O_ANY_CAN_PAIR))
					||
					(db[i] == ')' && (mask[i] & UMASK_C_ANY_CAN_PAIR))
					)) {
				return 0;
			}
		}
		else {
			//db include <> symbols for non-canonical base pairs.
			int OK = 0;
			switch (db[i]) {
			case '.':
				OK = (mask[i] & UMASK_NT_CAN_BE_FREE);
				break;
			case '<':
				OK = (mask[i] & UMASK_O_NC_CAN_PAIR);
				break;
			case '>':
				OK = (mask[i] & UMASK_C_NC_CAN_PAIR);
				break;
			case '(':
				OK = (mask[i] & UMASK_O_C_CAN_PAIR);
				break;
			case ')':
				OK = (mask[i] & UMASK_C_C_CAN_PAIR);
				break;
			default: {
				fprintf(stderr, "ERROR: Output dot bracket contains unexpected character (%c) in file %s at line %i.\n", db[i], __FILE__, __LINE__);
				exit(ABORT_EXIT);
			}
			}
			if (0 == OK) {
				return 0;
			}
		}
	}
	return 1;
}



/*
 * Given a dot bracket as char 2D structure passed in,
 * this code returns a alloc_aligneded char* representing all output masks (from the global omasks set in main() ) that match the provided dot bracket as char passed in
 * The returned char* comprises one label ( identifier from output_mask_t ) for each matching mask.
 * NOTE: caller frees the return value
 */
char * maskLabelsForOmasks(const char *db) {
	int hitlen;
	hitlen = 0;

	for (int i = 0; i < omasksSize; ++i) {
		if (doesDBmatchOutputMask(db, i))//omasks[i].mask_idb,omasks[i].type ) )
			hitlen += 1 + strlen(omasks[i].identifier);
	}
	if (hitlen == 0) return NULL;
	char *labels = (char*)alloc_aligned(sizeof(char)*(1 + hitlen));
	if (labels == NULL) {
		fprintf(stderr, "ERROR: OUT OF MEMORY in file %s at line %i\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
	labels[0] = 0;
	for (int i = 0; i < omasksSize; ++i) {
		if (doesDBmatchOutputMask(db, i)) {//omasks[i].mask_idb, omasks[i].type) ){
			strcat(labels, omasks[i].identifier);
			strcat(labels, " ");
		}
	}
	return labels;
}

/*
 * Given an dot bracket as char 2D structure passed in and a vector of int representing matching output masks,
 * update the vector by considering all output masks (+1 if omasks[i] matches db).
 */
void update_matches_vector(const char *db, int *matches) {
	for (int i = 0; i < omasksSize; ++i) {
		if (doesDBmatchOutputMask(db, i))
			++matches[i];
	}
}
/*
 * companion to update_matches_vector()
 * create a string filled with labels of mathing output masks.
 * NOTE: caller frees return string.
 */
char *labels_for_matches(const int * matches) {
	int hitlen = 0;

	for (int i = 0; i < omasksSize; ++i) {
		if (matches[i] > 0)
			hitlen += 1 + (int)strlen(omasks[i].identifier);
	}
	if (hitlen == 0) {
		char* empty = (char*)alloc_aligned(sizeof(char));
		empty[0] = (char)0;
		return empty;
	}
	char *labels = (char*)alloc_aligned(sizeof(char)*(1 + hitlen));
	if (labels == NULL) {
		fprintf(stderr, "ERROR: OUT OF MEMORY in file %s at line %i\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
	labels[0] = (char)0;
	for (int i = 0; i < omasksSize; ++i) {
		if (matches[i] > 0) {
			strcat(labels, omasks[i].identifier);
			strcat(labels, " ");
		}
	}
	return labels;

}

/*
 * user supplied balanced mask is converted to an integer dot bracket representation and checked for integrity at the same time.
 * masks are composed of chars from ().x where ( must match ), x is -1 and . is self index (ie: idb[i]=i);
 */
int * string2balanded_idb(const char *user_mask) {
	int len = (int)strlen(user_mask);
	//int stack[len];
	int *stack = alloc_aligned(sizeof(int) * len);
	int sp = -1;
	int *idb = (int*)alloc_aligned(sizeof(int)*len);
	if (idb == NULL) {
		fprintf(stderr, "ERROR: OUT OF MEMORY in file %s at line %i.\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
	for (int i = 0; i < len; ++i) {
		if (user_mask[i] == '(') {
			stack[++sp] = i;
		}
		else if (user_mask[i] == ')') {
			idb[i] = stack[sp];
			idb[stack[sp]] = i;
			--sp;
		}
		else if (user_mask[i] == '.') {
			idb[i] = i;
		}
		else { //if( user_mask[i]=='x' ) {
			idb[i] = -1;
		} //else {
		  //  fprintf(stderr, "ERROR: mask comprises an unknown character (%c).\n",user_mask[i]);
		  //  exit( ABORT_EXIT );
		  // }
	}
	if (sp != -1) {
		fprintf(stderr, "ERROR: unbalanced mask %s.\n", user_mask);
		exit(ABORT_EXIT);
	}
	free(stack);
	return idb;
}



/*
 * user supplied mask is converted to an integer representation of unbalanced masks
 *
 * BALANCED MASKS are composed of x .()
 * UNBALANCED MASKS are composed of x .()| -[]+ _<>! whose meanings are as follows:
 *
 * x is dont care
 *
 * . is not paired
 * ( is forward paired
 * ) is reverse paired
 * | is paired
 *
 * - is not involved in a canonical pair
 * [ is forward canonical paired
 * ] is reverse canonical paired
 * + is involved in canonical pairing
 *
 * _ is not involved in a non-canonical pair
 * < is forward non-canonicaly paired
 * > is reverse non-canonicaly paired
 * ! is involved in a non-canonical pair
 *
 * p is NOT involved in a reverse pair
 * q is NOT involved in a forward pair
 */
int * string2UMASK_idb(const char *user_mask) {
	int len = (int)strlen(user_mask);
	int *idb = (int*)alloc_aligned(sizeof(int)*len);
	for (int i = 0; i < len; ++i) {
		if (user_mask[i] == '(') {
			idb[i] = UMASK_FORWARD_PAIRED;
		}
		else if (user_mask[i] == ')') {
			idb[i] = UMASK_BACK_PAIRED;
		}
		else if (user_mask[i] == '.') {
			idb[i] = UMASK_UNPAIRED;
		}
		else if (user_mask[i] == 'x') {
			idb[i] = UMASK_DONTCARE;
		}
		else if (user_mask[i] == '|') {
			idb[i] = UMASK_PAIRED;
		}
		else if (user_mask[i] == '-') {
			idb[i] = UMASK_NOT_C_PAIRED;
		}
		else if (user_mask[i] == '[') {
			idb[i] = UMASK_FORWARD_C_PAIRED;
		}
		else if (user_mask[i] == ']') {
			idb[i] = UMASK_BACK_C_PAIRED;
		}
		else if (user_mask[i] == '+') {
			idb[i] = UMASK_C_PAIRED;
		}
		else if (user_mask[i] == '_') {
			idb[i] = UMASK_NOT_NC_PAIRED;
		}
		else if (user_mask[i] == '<') {
			idb[i] = UMASK_FORWARD_NC_PAIRED;
		}
		else if (user_mask[i] == '>') {
			idb[i] = UMASK_BACK_NC_PAIRED;
		}
		else if (user_mask[i] == '!') {
			idb[i] = UMASK_NC_PAIRED;
		}
		else if (user_mask[i] == 'p') {
			idb[i] = UMASK_NOT_REVERSE_PAIRED;
		}
		else if (user_mask[i] == 'q') {
			idb[i] = UMASK_NOT_FORWARD_PAIRED;
		}
		else {
			fprintf(stderr, "ERROR: mask comprises an unknown character (%c).\n", user_mask[i]);
			exit(ABORT_EXIT);
		}
	}
	return idb;
}


/*
 * given that full masks have been read (both balanced and unbalanced), scan the balanced one and replace all
 * characters that are not in .()x with x, updating the unbalanced mask at the same time.
 */

const int NUMBER_OF_ALLOWED_MASK_SYMBOLS = 15;
const char allowedMaskSymbols[] = ")(.x|-[]+_<>!pq";
int symbol2idx(char s) { for (int i = 0; i < NUMBER_OF_ALLOWED_MASK_SYMBOLS; ++i) if (allowedMaskSymbols[i] == s) return i; return -1; }

const int CONCILIATION_MATRIX[] = {
	/*             )  (  .  x  |  -  [  ]  +  _  <  >  !  p  q */
	/*             0  1  2  3  4  5  6  7  8  9 10 11 12 13 14*/
	/* )  0 */     0,-1,-1, 0, 0,11,-1, 7, 7, 7,-1,11,-1,-1, 0,
	/* (  1 */    -1, 1,-1, 1, 1,10, 6,-1, 6, 6,10,-1,-1, 1,-1,
	/* .  2 */    -1,-1, 2, 2,-1, 5,-1,-1,-1, 9,-1,-1,-1, 2, 2,
	/* x  3 */     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,
	/* |  4 */     0, 1,-1, 4, 4,12, 6, 7, 8, 8,10,11,12,13,14,
	/* -  5 */    11,10, 5, 5,12, 5,-1,-1,-1, 2,10,11,12,-1,-1,
	/* [  6 */    -1, 6,-1, 6, 6,-1, 6,-1, 6, 6,-1,-1,-1, 6,-1,
	/* ]  7 */     7,-1,-1, 7, 7,-1,-1, 7, 7, 7,-1,-1,-1,-1, 7,
	/* +  8 */     7, 6,-1, 8, 8,-1, 6, 7, 8, 8,-1,-1,-1,-1,-1,
	/* _  9 */     7, 6, 9, 9, 4, 2, 6, 7, 8, 9,-1,-1,-1,-1,-1,
	/* < 10 */    -1,10,-1,10,10,10,-1,-1,-1,-1,10,-1,10,10,-1,
	/* > 11 */    11,-1,-1,11,11,11,-1,-1,-1,-1,-1,11,11,-1,11,
	/* ! 12 */    11,10,-1,12,12,12,-1,-1,-1,-1,10,11,12,-1,-1,
	/* p 13 */    -1, 1, 2,13,13,-1, 6,-1,-1,-1,10,-1,-1,13,-1,
	/* q 14 */     0,-1, 2,14,14,-1,-1, 7,-1,-1,-1,11,-1,-1,14
};

void conciliate_full_masks(int seqLen, char *bmask, char *umask) { //b:balanced, u:unbalanced
	if (NULL == bmask) {
		return;
	}

	int errors = 0;
	//char bad_positions[ seqLen +1 ];
	char *bad_positions = alloc_aligned(sizeof(char) * (seqLen + 1));
	bad_positions[seqLen] = 0;

	char * bmask_orig = strdup(bmask);
	char * umask_orig = strdup(umask);

	for (int i = 0; i < seqLen; ++i) {

		bad_positions[i] = ' ';

		int f = symbol2idx(bmask[i]);
		int ub = symbol2idx(umask[i]);


		if (f < 0) {
			fprintf(stderr, "ERROR: Unknown symbol (%c) in full mask.\n", bmask[i]);
			exit(ABORT_EXIT);
		}
		if (ub < 0) {
			fprintf(stderr, "ERROR: Unknown symbol (%c) in full unbalanced mask.\n", umask[i]);
			exit(ABORT_EXIT);
		}

		int concil = CONCILIATION_MATRIX[f + (NUMBER_OF_ALLOWED_MASK_SYMBOLS * ub)];
		if (concil < 0) {
			++errors;
			bad_positions[i] = '^';
		}
		else {
			umask[i] = allowedMaskSymbols[concil];
			if (f > symbol2idx('x')) {
				bmask[i] = 'x';
			}
		}
	}
	if (errors > 0) {
		fprintf(stderr, "ERROR: Can not reconcile %i constraint(s) from -m with those from -um.\n%s\n%s\n%s\n",
			errors,
			bmask_orig,
			umask_orig,
			bad_positions);
		exit(ABORT_EXIT);
	}
	free(umask_orig);
	free(bmask_orig);
	free(bad_positions);
}


/*
 * main calls this function for each user supplied omask.
 * each mask has a name and a mask.
 * the mask is appended to the global omasks
 */
void add_output_user_mask(const char *identifier, const char *user_mask, int maskType) {
	if (omasks == NULL) {
		omasksCapacity = 8;
		omasksSize = 0;
		omasks = (struct output_mask_t*)alloc_aligned(sizeof(struct output_mask_t) * omasksCapacity);
	}
	if (omasksCapacity <= omasksSize) {
		omasksCapacity *= 2;
		omasks = (struct output_mask_t*)realloc(omasks, sizeof(struct output_mask_t) * omasksCapacity);
	}
	if (omasks == NULL) {
		fprintf(stderr, "ERROR: OUT OF MEMORY in file %s at line %i\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
	omasks[omasksSize].identifier = strdup(identifier);
	omasks[omasksSize].balanced_idb = NULL;
	omasks[omasksSize].type = maskType;
	if (maskType == OUTPUT_MASK_BALANCED) {
		//starting with f26, a balanced mask is also an unbalanced one.
		omasks[omasksSize].balanced_idb = (int*)string2balanded_idb(user_mask);
		omasks[omasksSize].mask_idb = (int*)string2UMASK_idb(user_mask);
	}
	else if (maskType == OUTPUT_MASK_UNBALANCED) {
		omasks[omasksSize].mask_idb = (int*)string2UMASK_idb(user_mask);
	}
	else {
		fprintf(stderr, "ERROR: Unknown type of user mask in file %s at line %i.\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}

	++omasksSize;
}


/*
 * prints a 2D array to screen.
 * --> only works on arrays alloc'd using aXalloc
 */
void show2DArray(struct array_t A) {
	for (int i = 0; i < A.DIM[0]; ++i) {
		for (int j = 0; j < A.DIM[1]; ++j) {
			if (a2d(A, i, j) >= LARGEFLOAT)
				printf("  --  ");
			else
				printf("%+5.2f ", a2d(A, i, j));
		}
		printf("\n");
	}
}

/*
 * prints a 3D array to screen.
 * --> only works on arrays alloc'd using aXalloc
 */
void show3DArray(struct array_t A) {
	for (int k = 0; k < A.DIM[2]; ++k) {
		printf("\n%i\n", k);
		for (int i = 0; i < A.DIM[0]; ++i) {
			for (int j = 0; j < A.DIM[1]; ++j) {
				if (a3d(A, i, j, k) >= LARGEFLOAT)
					printf("  --  ");
				else
					printf("%+5.2f ", a3d(A, i, j, k));
			}
			printf("\n");
		}
	}
}

void showUpper2DArray(struct array_t A, int N) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i > j) {
				printf("  LL  ");
			}
			else if (i == j) {
				printf("  MM  ");
			}
			else {
				if (a2dh(A, i, j) >= LARGEFLOAT)
					printf("  --  ");
				else
					printf("%+5.2f ", a2dh(A, i, j));
			}
		}
		printf("\n");
	}
	printf("\n");
}

void showUpper3DArray(struct array_t A, int N, int K) {
	for (int k = 0; k < K; ++k) {
		printf("\nk=%i\n", k);
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				if (i > j) {
					printf("  LL  ");
				}
				else if (i == j) {
					printf("  MM  ");
				}
				else {
					if (a3dh(A, i, j, k) >= LARGEFLOAT)
						printf("  --  ");
					else
						printf("%+5.2f ", a3dh(A, i, j, k));
				}
			}
			printf("\n");
		}
		printf("\n\n");
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////// FILES FACILITIES
//////////////////////////////// ENERGIES DATA STRUCTURES

char * TABLES_SOURCE_PATH = NULL;

struct array_t strands, //the length of each strand of each NCM is described. NCMID x [5'strand = 0, 3'strand = 1]
	energies, //energies for each NCM (sequence x NCMID) ((otherwise known as cycles)
	transitions, //transitions energies 4 dimensions NCMID x NCMID x 5'nt x 3'nt
	hinges, //energies of base pairs used in hinges calculations as well as in capping stems 5'nt x 3'nt
	junctions; //junctions energies NCMID x NCMID



	/*
	 * read floats from text data in fname to store, return number of reads
	 * don't mess with data files!
	 */
size_t readFloatTextDataFileDebug(const char* path, char *infname, float *store, int debug) {

	//fixme: protect against too small read buffer (store)

	char *fname = mergeToNewString(path, infname);

	FILE * fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("Could not find file %s\n", fname);
		exit(ABORT_EXIT);
	}
	size_t i = 0;
	while (fscanf(fp, "%f", store) != EOF) {
		if (0 && debug) {//ha ha
			printf("%3.4f ", *(store));
		}
		++store;
		++i;
	}

	fclose(fp);

	free(fname);

	return i;
}

//wrapper with debug set to 0
size_t readFloatTextDataFile(const char *path, char *fname, float *store) {
	return readFloatTextDataFileDebug(path, fname, store, 0);
}

void readFloatBinaryDataFile(const char * path, char *infname, float *store, size_t sz) {

	char *fname = mergeToNewString(path, infname);
#ifndef WindowsCompilation
	int fd = open(fname, O_RDONLY);
#else
	int fd = open(fname, O_RDONLY | O_BINARY);
#endif
	if (fd == -1) {
		printf("ERROR: Could not open file %s for reading.\n", fname);
		exit(ABORT_EXIT);
	}
	sz *= sizeof(float);
	while (sz > 0) {
		size_t nread = read(fd, store, sz);
		if (nread == 0) {
			printf("ERROR: Could not read file %s or file is corrupted.\n", fname);
			exit(ABORT_EXIT);
		}
		store += nread;
		sz -= nread;
	}
	close(fd);
	free(fname);
}

/* dump data to binary file. Crash on error */
void dumpFloatAsBinaryToFile(const char *path, char *infname, float *data, size_t sz) {

	char *fname = mergeToNewString(path, infname);

#ifndef WindowsCompilation 
	int fd = open(fname, O_CREAT | O_WRONLY | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
#else 
	int fd = open(fname, O_CREAT | O_BINARY | O_NOINHERIT | O_RDWR | O_TRUNC);// | O_WRONLY | O_TRUNC);
#endif
	if (fd == -1) {
		printf("ERROR: Can not create new file %s.\n", fname);
		exit(ABORT_EXIT);
	}
	sz *= sizeof(float);
	size_t writen = write(fd, data, sz);
	if (writen != sz) { //refine this by placing the write in a loop please
		printf("ERROR: Can not write to file %s.\n", fname);
		exit(ABORT_EXIT);
	}
	close(fd);
	free(fname);
}



/*
 * These tables are digested by R and dumped as csv files
 * The R code is in a file called ~dallaire/code/mcscratch/source/mcscratch_load_tables-V20.R
 * The original MCfold tables from Parisien version are found under ~dallaire/code/mcscratch/MCFOLD-DB-orig/
 * They are digested manually to those found under ~dallaire/code/mcscratch/tables/ *.scratch
 * The R processed files are dumped to ~dallaire/code/mcscratch/tables/ *.f2.csv
 * Finally the binary versions are dumped to ./tables/ *.f2
 */
int LOADEDENERGYTABLES = 0;
void loadEnergyTables() {

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ALLOCATE MEMORY FOR ARRAYS
	aXalloc(&strands, 2, NCMs, 2, 0, 0, 0);
	aXalloc(&energies, 2, 87360, NCMs, 0, 0, 0);
	aXalloc(&hinges, 2, 4, 4, 0, 0, 0);
	aXalloc(&junctions, 2, 21, 21, 0, 0, 0);
	aXalloc(&transitions, 4, 21, 21, 4, 4, 0);


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! load binary tables if available
	char *testFname = mergeToNewString(TABLES_SOURCE_PATH, "strands.f8");
	int fd = open(testFname, O_RDONLY);
	free(testFname);
	if (fd != -1) {
		close(fd);

		readFloatBinaryDataFile(TABLES_SOURCE_PATH, "strands.f8", strands.data, strands.DIM[0] * strands.DIM[1]);
		readFloatBinaryDataFile(TABLES_SOURCE_PATH, "transitions.f8", transitions.data, transitions.DIM[0] * transitions.DIM[1] * transitions.DIM[2] * transitions.DIM[3]);
		readFloatBinaryDataFile(TABLES_SOURCE_PATH, "hinges.f8", hinges.data, hinges.DIM[0] * hinges.DIM[1]);
		readFloatBinaryDataFile(TABLES_SOURCE_PATH, "junctions.f8", junctions.data, junctions.DIM[0] * junctions.DIM[1]);
		readFloatBinaryDataFile(TABLES_SOURCE_PATH, "energies.f8", energies.data, energies.DIM[0] * energies.DIM[1]);

	}
	else {
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! otherwise, load text tables and generate binary tables

		size_t reads = readFloatTextDataFile(TABLES_SOURCE_PATH, "strands.f2.csv", strands.data);
		if (reads != strands.DIM[0] * strands.DIM[1]) {
			fprintf(stderr, "ERROR: BAD FILE strands.f2.csv\n");
			exit(ABORT_EXIT);
		}
		dumpFloatAsBinaryToFile(TABLES_SOURCE_PATH, "strands.f8", strands.data, strands.DIM[0] * strands.DIM[1]);

		reads = readFloatTextDataFile(TABLES_SOURCE_PATH, "energies.f2.csv", energies.data);
		if (reads != energies.DIM[0] * energies.DIM[1]) {
			fprintf(stderr, "ERROR: BAD FILE energies.f2.csv\n");
			exit(ABORT_EXIT);
		}
		dumpFloatAsBinaryToFile(TABLES_SOURCE_PATH, "energies.f8", energies.data, energies.DIM[0] * energies.DIM[1]);

		reads = readFloatTextDataFile(TABLES_SOURCE_PATH, "transitions.f2.csv", transitions.data);
		printf("INFO: done with reading\n");
		if (reads != transitions.DIM[0] * transitions.DIM[1] * transitions.DIM[2] * transitions.DIM[3]) {
			fprintf(stderr, "ERROR: BAD FILE transitions.f2.csv\n");
			exit(ABORT_EXIT);
		}
		dumpFloatAsBinaryToFile(TABLES_SOURCE_PATH, "transitions.f8", transitions.data, transitions.DIM[0] * transitions.DIM[1] * transitions.DIM[2] * transitions.DIM[3]);

		reads = readFloatTextDataFileDebug(TABLES_SOURCE_PATH, "hinges.f2.csv", hinges.data, 1);
		if (reads != hinges.DIM[0] * hinges.DIM[1]) {
			fprintf(stderr, "ERROR: BAD FILE hinges.f2.csv\n");
			exit(ABORT_EXIT);
		}
		dumpFloatAsBinaryToFile(TABLES_SOURCE_PATH, "hinges.f8", hinges.data, hinges.DIM[0] * hinges.DIM[1]);

		reads = readFloatTextDataFile(TABLES_SOURCE_PATH, "junctions.f2.csv", junctions.data);
		if (reads != junctions.DIM[0] * junctions.DIM[1]) {
			fprintf(stderr, "ERROR: BAD FILE junctions.f2.csv\n");
			exit(ABORT_EXIT);
		}
		dumpFloatAsBinaryToFile(TABLES_SOURCE_PATH, "junctions.f8", junctions.data, junctions.DIM[0] * junctions.DIM[1]);
	}

	LOADEDENERGYTABLES = 1; //done loading

}



/////////////////////////////////////////////////////////////////////////////////////////////// forward: INSIDE ALGORITHM


//ATTENTION: UNBALANCED_MASK contains values defined earlier with #define UMASK_xxx whereas FULL_MASK has values as follows: -1 if don't care, FULL_MASK[i]==i if unpaired only and [i] == j if i and j can pair.
//and these macros are used to verify that a nuleotide pair can pair or not


#define __upair( i, j, seqAsInt ) ( (NULL==UNBALANCED_MASK) || ( \
( (5==((seqAsInt)[(i)] + (seqAsInt)[(j)])) \
? ( ( UNBALANCED_MASK[(i)]&(UMASK_O_C_CAN_PAIR) ) \
     && \
     ( UNBALANCED_MASK[(j)]&(UMASK_C_C_CAN_PAIR) ) \
   ) \
: \
  ( ( UNBALANCED_MASK[(i)]&(UMASK_O_NC_CAN_PAIR) ) \
     && \
    ( UNBALANCED_MASK[(j)]&(UMASK_C_NC_CAN_PAIR) ) \
  ) \
) \
) )

#define mpair( i, j, seqAsInt ) ( ((NULL!=FULL_MASK)?( (FULL_MASK[(i)]==(j))?(1):( (FULL_MASK[(i)]==-1) && (FULL_MASK[(j)]==-1) ) ):(1)) && ( ((i)!=(j)) && ( __upair(((i)<(j)?(i):(j)),((i)<(j)?(j):(i)), (seqAsInt) ) ) ) )

#define __ufree( i ) ( (NULL==UNBALANCED_MASK) || (UNBALANCED_MASK[(i)] & UMASK_NT_CAN_BE_FREE) )
#define mfree( i ) ( ( (NULL==FULL_MASK) || ((FULL_MASK[(i)] == i) || (FULL_MASK[(i)]==-1))   ) && ( __ufree( (i) ) ) )


//inline
int mfree_loop(int i, int j) {
	for (int a = (i); a <= (j); ++a) {
		if (!(mfree(a)))
			return 0;;
	}
	return 1;

	//    int ok = 0;
	//    for(int a=(i);a<=(j);++a){
	//        ok+=(mfree(a)>0);
	//    }
	//    return ok == (j-i+1);
}


void print_mask_as_table(int N, int *seqAsInt) {

	printf(" ");
	for (int j = 0; j < N; ++j) {
		printf("%c", mfree(j) ? 'f' : ' ');
	}
	printf("\n");
	for (int i = 0; i < N; ++i) {
		printf("%c", mfree(i) ? 'f' : ' ');
		for (int j = 0; j < N; ++j) {
			printf("%c", mpair(i, j, seqAsInt) ? 'p' : ' ');
		}
		printf("\n");
	}
}

/* ////////////////////////////// STRUCTURES FOR THE FORWARD PASS
 * indexes on the sequence : 0--------------------i--p----------q--j----------(n-1)
 *
 * in the tables (F, L, R, RO, MB and MBO):
 *
 *                                    -- j --
 *
 *            0--------------------i--p----------q--j----------(n-1)
 *            |                                                |
 *            |                                                |
 *            |                           (i,j)                |
 *            |                                                |
 *      |     i                             i<j                i
 *      |     |                                                |
 *            p                                                p
 *      i     |                                                |
 *            |                                                |
 *      |     |                                                |
 *      |     |   lower left: nothing                          |
 *            q                                                q
 *            |                                                |
 *            j         i>j                                    j
 *            |                                                |
 *            0--------------------i--p----------q--j----------(n-1)
 *
 *   All tables use the upper right part.
 *
 *forward part of the algorithm
 * seqAsInt is the sequence of nucleotides coded as integers (A=1, C=2, G=3, U=4)
 * N is the length of the sequence
 * W is the window size (MIND<W<=N) (W=N for full sequence folding)
 * F is a table that will receive the forward scores (of size NxNxNCMs)
 * L and R are the Best Closing tables used to receive the value for left/right structures (size NxN)
 * MB is the MultiBranch table used to receive the values for right structures (size NxN)
 * To understand MBO and RO you will probably need to understand flashfolds gramar.
 * this routine is O( N^3 x NCMs )
 */
inline float minf(float a, float b) { return a < b ? a : b; }
void forward(const int* seqAsInt, const int N, const int W, struct array_t F, struct array_t L, struct array_t R, struct array_t RO, struct array_t MB, struct array_t MBO) {

	assert(N >= W);

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX all values in F, BC, MB, ...
	init3HArray(F, LARGEFLOAT);
	initHArray(L, LARGEFLOAT);
	initHArray(R, LARGEFLOAT);
	initHArray(RO, LARGEFLOAT);
	initHArray(MB, LARGEFLOAT);
	initHArray(MBO, LARGEFLOAT);

	if (DUPLEX_MODE != NONE) {
		//this scores a terminal loop with unpaired nts AA at DUPLEX_HAIRPIN_START-1
		// we get (..) where ( is the last nt from the 5' arm and ) is the first nt from the 3' arm and .. are AA
		// this is an NCM of type 2 (1 based) so we index using 1 (0 based)
		a3dh(F, DUPLEX_HAIRPIN_START - 1, DUPLEX_HAIRPIN_START + DUPLEX_HAIRPIN_LENGTH, 1) = (float)0.0;
	}

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! recurse on (i,j) fill by diagonal starting from the
	//                               one nearest from the main diagonal but still addressable by an NCM of type 1_1_X
	for (int d = MIND; d < W; ++d) {
		for (int i = 0; i < N - d; ++i) {
			int j = (i + d);// + 1);

			//---------- update F
			if (mpair(i, j, seqAsInt)) {
				for (int k = 0; k < NCMs; ++k) {
					if (DUPLEX_MODE == NCM4 && k != 4) continue;
					if (k < FIRST_TWO_STRANDS_NCM && (1 + j - i) <= MAX_NCM_TERMINAL_LOOP) { //1 strand NCM (terminal NCM)
						if (DUPLEX_MODE == NONE
							&& ((i + 1 < j - 1 && mfree_loop(i + 1, j - 1)) || (i + 1 == j - 1 && mfree(i + 1)))
							) {
							int idx = seq2idx(seqAsInt, i, j, -1, -1);
							float kNCMe = a2d(energies, idx, k);
							a3dh(F, i, j, k) = kNCMe;
						}
					}
					else { //2 strands NCM
						int p = i + (a2d(strands, k, 0) - 1); //the 5' arm of NCM k.
						int q = j - (a2d(strands, k, 1) - 1); //the 3' arm of NCM k.
						if ((q - p > 1 && p > i && q < j)
							&& mpair(p, q, seqAsInt)
							&& (p - i < 2 || mfree_loop(i + 1, p - 1))
							&& (j - q < 2 || mfree_loop(q + 1, j - 1)))
						{ //can iterate on NCMs (kk) at (p,q)

							float kNCMe = a2d(energies, seq2idx(seqAsInt, i, p, q, j), k);

							//----- F from kk
							for (int kk = 0; kk < NCMs; ++kk) {
								float forward = a3dh(F, p, q, kk);
								float hinge = a2d(hinges, seqAsInt[p] - 1, seqAsInt[q] - 1); //base pair at (p,q); -1 because there are no 0 values in seqAsInt
								float nucleotideTax = reactivity[p] + reactivity[q];
								float transition = a4d(transitions, kk, k, seqAsInt[p] - 1, seqAsInt[q] - 1);
								float junction = a2d(junctions, k, kk);
								float join = junction + ((transition + hinge > (float)1.0) ? (float)1.0 : (transition + hinge));
								//float score = forward + join + kNCMe;
								float score = forward + join + kNCMe + nucleotideTax;
								if (DUPLEX_MODE != NONE && 1 == kk) {
									score = hinge + kNCMe + forward + nucleotideTax;
								}
								a3dh(F, i, j, k) = score < a3dh(F, i, j, k) ? score : a3dh(F, i, j, k);
							}//for all kk

							//----- F from MB

							if (DUPLEX_MODE < NCM)
								if ((q - p - 1) > MIN_MB && p + 1 > i && q - 1 < j) {
									float hinge = a2d(hinges, seqAsInt[p] - 1, seqAsInt[q] - 1);
									float nucleotideTax = reactivity[p] + reactivity[q];
									hinge += nucleotideTax;
									float branch_score = a2dh(MBO, p + 1, q - 1);
									if (DUPLEX_MODE == SIMPLE) branch_score = a2dh(RO, p + 1, q - 1);
									float score = kNCMe + hinge + branch_score;
									a3dh(F, i, j, k) = (score < a3dh(F, i, j, k)) ? (score) : (a3dh(F, i, j, k));
								}
						} //q-p>=MIND
					} //end 2 strands NCMs
				} //for all k (update F)
			} //if(mpair(i,j)

			//---------- update L
			float Lscore = mfree(i) ? a2dh(L, i + 1, j) : LARGEFLOAT;//a2dh(L,i,j) ; //carry from L[i+1,j] ( or .() )

			float ijhinge = a2d(hinges, seqAsInt[i] - 1, seqAsInt[j] - 1);
			float nucleotideTax = reactivity[i] + reactivity[j];
			ijhinge += nucleotideTax;
			if (mpair(i, j, seqAsInt)) {
				for (int k = ALLOW_STEMS_OF_ONE ? 0 : FIRST_TWO_STRANDS_NCM; k < NCMs; ++k) { //find best score from F[i,j]
					float score = a3dh(F, i, j, k) + ijhinge;
					Lscore = (score < Lscore ? score : Lscore);
				}
			}
			a2dh(L, i, j) = Lscore;

			//---------- update R
			float Rscore = mfree(j) ? a2dh(R, i, j - 1) : LARGEFLOAT;//a2dh(R,i,j); //().
			float RfromL = a2dh(L, i, j);   //.() and ()
			a2dh(R, i, j) = Rscore < RfromL ? Rscore : RfromL;

			//---------- update RO
			float NotNCM = LARGEFLOAT;
			//      if( i+4<j-0-MIND & mfree_loop(i,i+4) ) NotNCM = a2dh(R,i+4,j-0);
			//if( i+3<j-2-MIND & mfree_loop(i,i+3) & mfree_loop(j-2,j) ) NotNCM = NotNCM<a2dh(R,i+3,j-2)?NotNCM:a2dh(R,i+3,j-2);
			//if( i+2<j-3-MIND & mfree_loop(i,i+2) & mfree_loop(j-3,j) ) NotNCM = NotNCM<a2dh(R,i+2,j-3)?NotNCM:a2dh(R,i+2,j-3);
			//if( i+0<j-4-MIND & mfree_loop(j-4,j)) NotNCM = NotNCM<a2dh(R,i+0,j-4)?NotNCM:a2dh(R,i+0,j-4);//i+1,j-4 overlaps this one

			if ((i + 4 < j - 1) & mfree_loop(i, i + 3)) NotNCM = a2dh(R, i + 4, j - 0);
			if ((i + 3 < j - 2) & mfree_loop(i, i + 2) & mfree_loop(j - 1, j)) NotNCM = NotNCM < a2dh(R, i + 3, j - 2) ? NotNCM : a2dh(R, i + 3, j - 2);
			if ((i + 2 < j - 3) & mfree_loop(i, i + 1) & mfree_loop(j - 2, j)) NotNCM = NotNCM < a2dh(R, i + 2, j - 3) ? NotNCM : a2dh(R, i + 2, j - 3);
			//if( (i+0<j-4) & mfree_loop(j-4,j)) NotNCM = NotNCM<a2dh(R,i+0,j-4)?NotNCM:a2dh(R,i+0,j-4);//i+1,j-4 overlaps this one
			if ((i + 0 < j - 4) & mfree_loop(j - 3, j)) NotNCM = NotNCM < a2dh(R, i + 0, j - 4) ? NotNCM : a2dh(R, i + 0, j - 4);//i+1,j-4 overlaps this one
			a2dh(RO, i, j) = NotNCM;

			//----------- Long Terminal Loops contributions
			if (ALLOW_LONG_TERMINAL_LOOP & (j - i + 1 >= 5)) {
				if (mfree_loop(i, j)) {
					a2dh(RO, i, j) = a2dh(RO, i, j) < 0 ? a2dh(RO, i, j) : 0;
					a2dh(R, i, j) = a2dh(R, i, j) < 0 ? a2dh(R, i, j) : 0;
				}
			}


			//---------- update MB and MBO
			float score = LARGEFLOAT;

			for (int s = i + MIND; (s + 1) < (j - MIND); ++s) {
				float left = a2dh(L, i, s);
				float rightMB = a2dh(MB, s + 1, j); // L(i,s) + MB(s+1,j)
				float s_score = left + rightMB; //min for (i,s)(s+1,j)
				score = score < s_score ? score : s_score; //min for all s
			} //for all s (update MB)

			a2dh(MB, i, j) = score < a2dh(R, i, j) ? score : a2dh(R, i, j);
			a2dh(MBO, i, j) = score < a2dh(RO, i, j) ? score : a2dh(RO, i, j);;

		} //for all i
	} //for all d


} //forward()




/////////////////////////////////////////////////////////////////////////////////////////////// BACKTRACKING STACK

enum node_t { RIGHT = 0, LEFT = 1, K = 2, KK = 3, BRANCH = 4, TLRIGHT = 5 }; //types of the nodes on the stack (TLRIGHT is terminal long loop))

const int BLOCKED_STATE = INT_MAX; //a special value of status_t::subtype that indicates that an item on the stack should not be revisited

enum stackBehaviour_t { NORMALBEHAVIOUR = 0, RIGHTBEHAVIOUR = 1, LEFTBEHAVIOUR = 2 }; //stack pop() is defined in terms of the top element's behaviour

struct status_t {
	int i, j,
		oi, oj,   //init these to i, j
		k, kk, s; //init these to -1
	enum node_t type;
	int subtype; //init to -1
	enum stackBehaviour_t behaviour;
	float reverse, promise, forward, score; //ATTENTION: forward and score used for debugging only (do not rely on them in calculation)
};

/*
 * new states are initialized using this function
 * see the code in backtrack()
 */
void setStatus(struct status_t *S,
	int i, int j, int oi, int oj, int k, int kk, int s,
	enum node_t type, int subtype, enum stackBehaviour_t behaviour,
	float reverse, float promise) {

	S->i = i; S->j = j; S->oi = oi; S->oj = oj; S->k = k, S->kk = kk; S->s = s;
	S->type = type; S->subtype = subtype; S->behaviour = behaviour;
	S->reverse = reverse; S->promise = promise;
	S->forward = S->score = 0;  //<<<<<<<<<<<<<<<<<<<<<<<<<<<< watch out for this. Do not rely on score and forward values being set in all states
}

/*
 * Use with great care and don't use this outside of pop() and terminal states KK and TLRIGHT !
 * that is because behaviour and subtype in particular should not be duplicated liberaly on the stack
 */
void copyStatus(struct status_t *stack, int from, int to) {
	memcpy(stack + to, stack + from, sizeof(struct status_t));
}

/*
 * we push on the bottom stack using ++bsp (bottom stack pointer)
 * we push on the top stack using --tsp (top stack pointer)
 * but we ALWAYS POP USING the function pop()
 * pop() takes care of poping both stacks according to the behaviour or their top elements
 */
void pop(struct status_t *S, int *bsp, int *tsp) {
	switch (S[*bsp].behaviour)
	{
	case NORMALBEHAVIOUR:
		--(*bsp);
		break;

	case LEFTBEHAVIOUR:
		--(*bsp);
		++(*tsp);
		break;

	case RIGHTBEHAVIOUR:
		copyStatus(S, (*bsp)--, --(*tsp)); //move right node on the bottom stack to top of (or bottom since stack is reversed) top stack
		S[*tsp].subtype = -1; //need to reset the state of the right node after poping or it will prevent further backtraking
		S[*tsp].i = S[*tsp].oi; //reset i and j to their initial values to allow backtracking in RIGHT and LEFT nodes.
		S[*tsp].j = S[*tsp].oj;
		S[*tsp].s = -1; //reset s so BRANCH nodes can restart backtrack
		break;

	default:
		fprintf(stderr, "ERROR: Found a bug in file %s at line %i \n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
}

/*plots one stack node for debugging purposes.*/
static const char *IDs[8] = { "R ", "L ", "K ", "KK", "MB", "TLR" }; //table of nodes short identifiers for printouts
static const char *STEPNAME[8] = { "RIGHT", "LEFT", "INIT-STEM", "EXTENSION", "MULTI-BRANCH", "TERMINAL LONG LOOP" }; //Long names of derivations synchronized to state names in article
void showStatus(struct status_t *stack, const int bsp, const int tsp) { //print node to screen
	printf("\ti(%i),j(%i),oi(%i),oj(%i),s(%i),k(%i),kk(%i),type(%s),subtype(%i),behaviour(%i),bsp(%i),tsp(%i),reverse(%+5.3f), promise(%+5.3f), forward(%+5.3f), sum(%+5.3f),score(%+5.10f)\n",
		(stack + bsp)->i, (stack + bsp)->j, (stack + bsp)->oi, (stack + bsp)->oj, (stack + bsp)->s,
		(stack + bsp)->k, (stack + bsp)->kk, IDs[(stack + bsp)->type], (stack + bsp)->subtype, (stack + bsp)->behaviour, bsp, tsp,
		(stack + bsp)->reverse, (stack + bsp)->promise, (stack + bsp)->forward, (stack + bsp)->reverse + (stack + bsp)->promise + (stack + bsp)->forward,
		(stack + bsp)->score);
}


/////////////////////////////////////////////////////////////////////////////////////////////// UNSORTED SIMPLE DUPLEX SOLUTIONS OUTPUTER

int iscanon(const int i, const int j, const int *vi) { return(5 != ((vi)[(i)] + (vi)[(j)])); }
int allIsCanon(const int i, const int j, const int *vi) { return 0; }
int(*specialBracketsFor)(const int, const int, const int *) = &allIsCanon;


/*
 * Used only for the duplex mode.
 * You can customize it to draw a hybrid if you have time to loose.
 * This outputer cannot be intermixed with the other outputers like heap, graph and such.
 */
float SIMPLEDUPLEX_SolutionOutputer(struct status_t *stack, const int bsp, const float score, int I, int J, float threshold, const int *seqAsInt) {
	//Get length of dotb to output and start offset (I). I > 0 if scanning.
	int solLen = J - I + 1;
	//char dotb[ solLen +1 ]; //build the dot bracket solution in here
	char *dotb = alloc_aligned(sizeof(char) * (solLen + 1));
	for (int i = 0; i < solLen; ++i)
		dotb[i] = '.';
	dotb[solLen] = 0;


	//go down the stack, setting the dot brackets to dotb
	for (int sp = bsp; sp >= 0; --sp) {
		//showStatus( stack, sp, 0 );
		switch (stack[sp].type) {
		case RIGHT: //don't care
			break;
		case LEFT:  //don't care
			break;
		case K:     //don't care
			break;
		case KK:    //pair (i,j)

		{
			int c = specialBracketsFor(stack[sp].i - I, stack[sp].j - I, seqAsInt);

			dotb[stack[sp].i - I] = (c ? '<' : '(');
			dotb[stack[sp].j - I] = c ? '>' : ')';
			//if k internal then (p,q) also paired
			// k is NCMs+1 if terminal (you know, to prevent backtracking... see case KK in backtrack
			if (stack[sp].k < NCMs && stack[sp].k >= FIRST_TWO_STRANDS_NCM) {
				int p = stack[sp].i + (a2d(strands, stack[sp].k, 0) - 1); //the 5' arm of NCM k,
				int q = stack[sp].j - (a2d(strands, stack[sp].k, 1) - 1); //the 3' arm of NCM k.
				int c = specialBracketsFor(p - I, q - I, seqAsInt);
				dotb[p - I] = c ? '<' : '(';
				dotb[q - I] = c ? '>' : ')';
			}
		}
		break;
		case BRANCH: //don't care
			break;
		case TLRIGHT: //don't care
			break;
		default:    //report weird error
			fprintf(stderr, "ERROR: Found a BUG ( inexistant node type in SIMPLEDUPLEX_SolutionOutputer() ) at line %i in file %s.\n", __LINE__, __FILE__);
			exit(ABORT_EXIT);
			break;
		}
	}

	//insert a spacer in stdout dot bracket
	for (int i = DUPLEX_HAIRPIN_START; i < DUPLEX_HAIRPIN_START + DUPLEX_HAIRPIN_LENGTH; ++i) {
		dotb[i] = ' ';
	}
	printf("%s %+5.3f\n", dotb, score);
	free(dotb);
	return threshold;
}

/////////////////////////////////////////////////////////////////////////////////////////////// PARSE SOLUTIONS OUTPUTER
/*
 * Show the parse of a solution
 */
float parseSolutionOutputer(struct status_t *stack, const int bsp, const float score, int I, int J, float threshold, const int * seqAsInt) {
	//Get length of dotb to output and start offset (I). I > 0 if scanning.
	int solLen = J - I + 1;
	//char dotb[ solLen +1 ]; //build the dot bracket solution in here
	char *dotb = alloc_aligned(sizeof(char) * (solLen + 1));
	for (int i = 0; i < solLen; ++i)
		dotb[i] = '-';
	dotb[solLen] = 0;

	//go down the stack, setting the dot brackets to dotb
	for (int sp = 0; sp <= bsp; ++sp) {
		//showStatus( stack, sp, 0 );

		switch (stack[sp].type) {
		case RIGHT:
			for (int i = stack[sp].i - I; i <= stack[sp].j - I; ++i)dotb[i] = 'R';
			for (int i = stack[sp].j - I + 1; i <= stack[sp - 1].j - I; ++i) dotb[i] = '.';//'r';
			break;
		case LEFT:  for (int i = stack[sp].i - I; i <= stack[sp].j - I; ++i)dotb[i] = 'L';
			for (int i = stack[sp - 1].i - I; i < stack[sp].i - I; ++i) dotb[i] = '.';//'l'
			break;
		case K:     //don't care
			for (int i = stack[sp].i - I; i <= stack[sp].j - I; ++i)dotb[i] = 'I';
			dotb[stack[sp].i - I] = '(';
			dotb[stack[sp].j - I] = ')';
			break;
		case KK:    //pair (i,j)
		{ for (int i = stack[sp].i - I; i <= stack[sp].j - I; ++i)dotb[i] = 'X';
		dotb[stack[sp].i - I] = '(';
		dotb[stack[sp].j - I] = ')';
		//if k internal then (p,q) also paired
		// k is NCMs+1 if terminal (you know, to prevent backtracking... see case KK in backtrack
		if (stack[sp].k < NCMs && stack[sp].k >= FIRST_TWO_STRANDS_NCM) {
			int p = stack[sp].i + (a2d(strands, stack[sp].k, 0) - 1); //the 5' arm of NCM k,
			int q = stack[sp].j - (a2d(strands, stack[sp].k, 1) - 1); //the 3' arm of NCM k.

			dotb[p - I] = '(';
			dotb[q - I] = ')';

			for (int i = stack[sp].i - I + 1; i < p; ++i) dotb[i] = '.';//'5';
			for (int i = q + 1; i < stack[sp].j - I; ++i) dotb[i] = '.';//'3';

		}
		else {
			for (int i = stack[sp].i - I + 1; i < stack[sp].j - I; ++i) {
				dotb[i] = '.';//'c';
			}
		}
		}
		break;
		case BRANCH: for (int i = stack[sp].i - I; i <= stack[sp].j - I; ++i)dotb[i] = 'M';
			break;
		case TLRIGHT: for (int i = stack[sp].i - I; i <= stack[sp].j - I; ++i)dotb[i] = '.';//'t';
			break;
		default:    //report weird error
			fprintf(stderr, "ERROR: Found a BUG ( inexistant node type in dotBracketSolutionOutputer() ) at line %i in file %s.\n", __LINE__, __FILE__);
			exit(ABORT_EXIT);
			break;
		}
		printf("%s %s\n", dotb, STEPNAME[stack[sp].type]);
	}
	free(dotb);
	return threshold;
}



/////////////////////////////////////////////////////////////////////////////////////////////// UNSORTED SOLUTIONS OUTPUTER

int PRINT_PARSE_BEFORE_DOT_BRACKET_SOLUTIONS = 0;

/*
 * Generate dot bracket
 * Output it to stdout followed by its free folding energy
 * I and J are respectively the first and last nucleotide values (normally this is I=0, J=N-1
 */
float dotBracketSolutionOutputer(struct status_t *stack, const int bsp, const float score, int I, int J, float threshold, const int * seqAsInt) {
	if (PRINT_PARSE_BEFORE_DOT_BRACKET_SOLUTIONS)
		parseSolutionOutputer(stack, bsp, score, I, J, threshold, seqAsInt);

	//Get length of dotb to output and start offset (I). I > 0 if scanning.
	int solLen = J - I + 1;
	//char dotb[ solLen +1 ]; //build the dot bracket solution in here
	char *dotb = alloc_aligned(sizeof(char) * (solLen + 1));
	for (int i = 0; i < solLen; ++i)
		dotb[i] = '.';
	dotb[solLen] = 0;

	//go down the stack, setting the dot brackets to dotb
	for (int sp = bsp; sp >= 0; --sp) {
		//showStatus( stack, sp, 0 );
		switch (stack[sp].type) {
		case RIGHT: //don't care
			break;
		case LEFT:  //don't care
			break;
		case K:     //don't care
			break;
		case KK:    //pair (i,j)

		{
			int c = specialBracketsFor(stack[sp].i - I, stack[sp].j - I, seqAsInt);

			dotb[stack[sp].i - I] = (c ? '<' : '(');
			dotb[stack[sp].j - I] = c ? '>' : ')';
			//if k internal then (p,q) also paired
			// k is NCMs+1 if terminal (you know, to prevent backtracking... see case KK in backtrack
			if (stack[sp].k < NCMs && stack[sp].k >= FIRST_TWO_STRANDS_NCM) {
				int p = stack[sp].i + (a2d(strands, stack[sp].k, 0) - 1); //the 5' arm of NCM k,
				int q = stack[sp].j - (a2d(strands, stack[sp].k, 1) - 1); //the 3' arm of NCM k.
				int c = specialBracketsFor(p - I, q - I, seqAsInt);
				dotb[p - I] = c ? '<' : '(';
				dotb[q - I] = c ? '>' : ')';
			}
		}
		break;
		case BRANCH: //don't care
			break;
		case TLRIGHT: //don't care
			break;
		default:    //report weird error
			fprintf(stderr, "ERROR: Found a BUG ( inexistant node type in dotBracketSolutionOutputer() ) at line %i in file %s.\n", __LINE__, __FILE__);
			exit(ABORT_EXIT);
			break;
		}
	}

	char * maskLabels = maskLabelsForOmasks(dotb);
	printf("%s %+5.3f %s\n", dotb, score, maskLabels ? maskLabels : "");
	free(maskLabels);
	free(dotb);
	return threshold;
}


/////////////////////////////////////////////////////////////////////////////////////////////// SORTED SOLUTIONS OUTPUTER
/*
 * cummulate the solutions so that they can be sorted before output
 * in the global table solArray of scored_db_t items
 * current size of solArray is SolCapacity
 * current number of items in solArray is solNum
 * scored_db::db is a character* containing the solutions as dot brackets
 */
struct scored_db_t {
	float score;
	char *db;
};

size_t solNum = 0;

size_t solCapacity = 0;

struct scored_db_t *solArray = NULL;

void free_sorted() {
	for (size_t i = 0; i < solNum; ++i) {
		free(solArray[i].db);
	}
	if (solCapacity > 0) {
		free(solArray);
	}
	solArray = NULL;
	solCapacity = 0;
	solNum = 0;
}

/*
 * Compare two entries of type scored_db_t to determine which comes first (has lower score)
 * Utility usable by qsort -1 if scored_db1 < scored_db2, ...
 */
int cmpScoredDB(const void *one, const void *two) {
	struct scored_db_t *sdb1 = (struct scored_db_t *)one;
	struct scored_db_t *sdb2 = (struct scored_db_t *)two;
	if (sdb1->score < sdb2->score) return -1;
	if (sdb1->score > sdb2->score) return 1;
	return 0;
}
/*
 * Once all solutions bellow or equal to threshold have been collected
 * sort them by score (best first) and output them to stdout.
 * output is followed by a shape level 5 representation
 * and by any identifiers associated to output masks currently defined
 */
void outputSortedDBSolutions() {
	qsort(solArray, solNum, sizeof(struct scored_db_t), &cmpScoredDB);
	for (size_t i = 0; i < solNum; ++i) {
		char *shape = db2shape(solArray[i].db);
		char * maskLabels = maskLabelsForOmasks(solArray[i].db);
		printf("%s %+5.3f %s %s\n", solArray[i].db, solArray[i].score, shape, maskLabels ? maskLabels : "");
		free(shape);
		free(maskLabels);
	}
}
/*
 * Invoqued each time a new solution is found during backtrack
 * -) compute dot bracket from stack
 * -) insert in solArray table (increase size if storage too small)
 */
float sortedDotBracketSolutionOutputer(struct status_t *stack, const int bsp,
	const float score, int I, int J, float threshold, const int * seqAsInt) {
	static const size_t initial_capacity = 64;
	if (solArray == NULL) {
		solCapacity = initial_capacity;
		solArray = (struct scored_db_t *)alloc_aligned(sizeof(struct scored_db_t) * solCapacity);
		if (solArray == NULL) {
			fprintf(stderr, "ERROR: OUT OF MEMORY (in file %s, at line %i).\n", __FILE__, __LINE__);
			exit(ABORT_EXIT);
		}
		solNum = 0;
	}

	if (solCapacity <= (solNum + 1)) {
		solCapacity = solCapacity * 2;

		solArray = (struct scored_db_t*) realloc((void*)solArray, sizeof(struct scored_db_t) * solCapacity);

		if (solArray == NULL) {
			fprintf(stderr, "ERROR: OUT OF MEMORY (in %s, at line %i).\n", __FILE__, __LINE__);
			exit(ABORT_EXIT);
		}
	}

	//Get length of dotb to output and start offset (I). I > 0 if scanning.
	int solLen = J - I + 1;
	char *dotb = (char*)alloc_aligned(sizeof(char*) * (solLen + 1));
	if (NULL == dotb) {
		fprintf(stderr, "ERROR: OUT OF MEMORY (in %s, at line %i).\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
	for (int i = 0; i < solLen; ++i)
		dotb[i] = '.';
	dotb[solLen] = 0;

	solArray[solNum].score = score;
	solArray[solNum].db = dotb;
	++solNum;

	//go down the stack, setting the dot brackets to dotb
	for (int sp = bsp; sp >= 0; --sp) {
		switch (stack[sp].type) {
		case RIGHT: //don't care
			break;
		case LEFT:  //don't care
			break;
		case K:     //don't care
			break;
		case KK:    //pair (i,j)

		{
			int c = specialBracketsFor(stack[sp].i - I, stack[sp].j - I, seqAsInt);

			dotb[stack[sp].i - I] = (c ? '<' : '(');
			dotb[stack[sp].j - I] = c ? '>' : ')';
			//if k internal then (p,q) also paired
			// k is NCMs+1 if terminal (you know, to prevent backtracking... see case KK in backtrack
			if (stack[sp].k < NCMs && stack[sp].k >= FIRST_TWO_STRANDS_NCM) {
				int p = stack[sp].i + (a2d(strands, stack[sp].k, 0) - 1); //the 5' arm of NCM k,
				int q = stack[sp].j - (a2d(strands, stack[sp].k, 1) - 1); //the 3' arm of NCM k.
				int c = specialBracketsFor(p - I, q - I, seqAsInt);
				dotb[p - I] = c ? '<' : '(';
				dotb[q - I] = c ? '>' : ')';
			}
		}
		break;
		case BRANCH: //don't care
			break;
		case TLRIGHT: //don't care
			break;

		default:    //report weird error
			fprintf(stderr, "ERROR: Found a BUG ( inexistant node type in sortedDotBracketSolutionOutputer() ) at lines %i in file %s.\n", __LINE__, __FILE__);
			exit(ABORT_EXIT);
			break;
		}
	}

	//printf( "%s %+5.3f\n", dotb, score );
	return threshold;
}


/////////////////////////////////////////////////////////////////////////////////////////////// GRAPH SOLUTIONS OUTPUTER
/*
 * cummulate the solutions in here (homologous to solArray under sorted solutions outputer
 *
 * here, db is a 2D structure loosely analogous to dot brackets but represented as vectors of short ints.
 * db[i] is the nucleotide with which i pairs.
 * if db[i] == i then i is unpaired.
 *
 * sometimes in this code, idb is used to represent these structures (integer db)
 */
struct graph_t {
	float score;
	int component; //init to self index in array
	unsigned short *db;
	int solLen; //length of db array (this is a waste because it is always N)
	int size; // size of component

	//data for metric tree
	int radius; //max differences between this and any entry in metric sub-tree: init to 0
	int next_bucket; //init to i+1
	int sibling; //init to -1
	int last_sibling; //init to -1
};


size_t graphSize = 0;

size_t graphCapacity = 0;

struct graph_t *graphArray = NULL;


/*
 * cleanup function
 */
void free_graph() {
	for (size_t i = 0; i < graphSize; ++i) {
		free(graphArray[i].db);
	}
	if (graphCapacity > 0)
		free(graphArray);

	graphArray = NULL;
	graphSize = 0;
	graphCapacity = 0;
}

/*
 * Compute the number of indices where a and b differ (a and b are idb) (hamming a,b)
 * the code can return max if the number of differences is >= max
 */
int db_differences(short *a, short *b, int N, int max) {
	int diff = 0;
	for (int i = 0; i < N; ++i) {
		diff += a[i] != b[i];
		if (diff >= max)
			return diff;
	}
	return diff;
}

#ifdef ___SIMD_COMPILATION__

typedef  short v8hi __attribute__((vector_size(16)));
typedef union { v8hi v;  short i[8]; } v8hi_u;

/*
 * same as db_differences but using SIMD
 */
int db_differences_SIMD(short *a, short *b, int N, int max) {

	v8hi *av = (v8hi*)a;
	v8hi *bv = (v8hi*)b;

	v8hi_u diff;
	diff.v = av[0] - av[0];

	static const v8hi ones = { 1,1,1,1, 1,1,1,1 };

	int i = 0;
	int res = 0;

	if (max > N) {

		for (; i < N; ++i) {
			v8hi tmp = (av[i] - bv[i]);
			diff.v += __builtin_ia32_pminuw128(tmp, ones);//must use  version !
		}

		res = diff.i[0] + diff.i[1] + diff.i[2] + diff.i[3]
			+ diff.i[4] + diff.i[5] + diff.i[6] + diff.i[7];

	}
	else {
		//regularly check if res > max

		div_t bl = div(N, 8);
		int Blocks = bl.quot + (bl.rem ? 1 : 0);

		for (; i < (max + 1) / 8; ++i) {
			v8hi tmp = (av[i] - bv[i]);
			diff.v += __builtin_ia32_pminuw128(tmp, ones);//must use  version !
		}

		res = diff.i[0] + diff.i[1] + diff.i[2] + diff.i[3]
			+ diff.i[4] + diff.i[5] + diff.i[6] + diff.i[7];

		for (; i < Blocks / 4; ++i) {
			v8hi tmp = (av[i] - bv[i]);
			diff.v += __builtin_ia32_pminuw128(tmp, ones);//must use  version !
		}

		res = diff.i[0] + diff.i[1] + diff.i[2] + diff.i[3]
			+ diff.i[4] + diff.i[5] + diff.i[6] + diff.i[7];

		if (res > max) return res;

		for (; i < Blocks / 2; ++i) {
			v8hi tmp = (av[i] - bv[i]);
			diff.v += __builtin_ia32_pminuw128(tmp, ones);//must use  version !
		}

		res = diff.i[0] + diff.i[1] + diff.i[2] + diff.i[3]
			+ diff.i[4] + diff.i[5] + diff.i[6] + diff.i[7];
		if (res > max) return res;

		for (; i < (3 * Blocks) / 4; ++i) {
			v8hi tmp = (av[i] - bv[i]);
			diff.v += __builtin_ia32_pminuw128(tmp, ones);//must use  version !
		}

		res = diff.i[0] + diff.i[1] + diff.i[2] + diff.i[3]
			+ diff.i[4] + diff.i[5] + diff.i[6] + diff.i[7];
		if (res > max) return res;

		for (; i < Blocks; ++i) {
			v8hi tmp = (av[i] - bv[i]);
			diff.v += __builtin_ia32_pminuw128(tmp, ones);//must use  version !
		}

		res = diff.i[0] + diff.i[1] + diff.i[2] + diff.i[3]
			+ diff.i[4] + diff.i[5] + diff.i[6] + diff.i[7];

	}
	return res;
}
#else
/*SIMD compilation flag not set. Fallback.*/
int db_differences_SIMD(short *a, short *b, int N, int max) {
	return db_differences(a, b, N, max);
}
#endif



int(*db_differences_fptr)(short *a, short *b, int N, int max) = &db_differences; //reset to ..SIMD by main if user option set so

/*
 * We use a metric data structure (list of buckets) to find graph connected component
 * where each component has its own bucket.
 * add_to_bucket fuses the buckets that contain node i and node j
 * The resulting bucket is that whose lowest energy structure is lowest. (sic).
 * Because the list of bucket is incrementaly built on the array graphArray and because this
 * array is sorted by energy of structures before the graph is built, the lower the index of
 * the first component of the bucket, the lower its energy.
 * after the fuse, the list of buckets (via nextBucket) is not up to date (because updating it would be O(n))
 * The next access to the previous bucket must splice the list.
 * this is done by checking x->nextBucket->component != x in disjointComponents() (see lower)
 * for explanations of radius, see disjointComponents()
 */
void add_to_bucket(int i, int j) {

	//update component number for i and for j
	i = graphArray[i].component;

	j = graphArray[j].component;

	//b has lower minimum free energy because it's index is lower than that of n
	//so b's bucket is augmented with n's bucket
	int b = i < j ? i : j;
	int n = i < j ? j : i;

	/*  assert(n>b);
	 assert( graphArray[n].component != graphArray[b].component);
	 assert( graphArray[n].sibling != n );
	 assert( graphArray[b].sibling != b );
	 */

	 //new size for resulting bucket
	graphArray[b].size += graphArray[n].size;

	//n in inserted into b
	int newBucket = graphArray[b].component;

	//find last sibling of b
	int l = b;
	if (graphArray[b].last_sibling != -1) //bucket of b contains only b?
		l = graphArray[b].last_sibling;

	//assert( l != n );

	//fusion would be O(1) if components numbers did not require updating
	graphArray[l].sibling = n;

	//update b->last_sibling
	l = n;
	if (graphArray[n].last_sibling != -1) //bucket n contains only n?
		l = graphArray[n].last_sibling;
	graphArray[b].last_sibling = l;

	//update component numbers and radius of bucket //O(n) ;-(
	int N = graphArray[b].solLen;
	while (n != -1) {
		graphArray[n].component = newBucket;

		float r = graphArray[b].radius;
		//mask to short because SIMD type is short (no consequences because this routine checks for equallity)
		float diff = (*db_differences)((short*)graphArray[b].db, (short*)graphArray[n].db, N, N + 1);
		r = r < diff ? diff : r;
		graphArray[b].radius = r;

		n = graphArray[n].sibling;
	}

}

/*
 * for type 1 graph connections we have the luxury of optimizing the algo and get better than n^2 behaviour
 * types refers to the type of connection linking two 2D.
 * type 0 allows a connection to be called when the number of differences between two intdb are <= 3.
 * type 1 forces migrating bulges to immediate neighbours.
 *
 * bucket list
 * uses a metric data structure to compute the disjoint components in graphArray.
 * added in f14
 * the idea is to skip the comparison of a new element with entries whose radiuses are too large for bucket i
 *
 * radius is a parameter set for each bucket.
 * radius is the maximum number of differences (db_differences) between a bucket 2D structure and any 2D structure found in its bucket.
 *
 */
void disjointComponents_type0() {
	//add new element one by one starting with 1 to graph already computed
	for (int i = 1; i < graphSize; ++i) {
		int bucket = 0;
		int nextBucket = (int)graphSize;
		while (bucket != (int)graphSize && bucket < i) {
			//Because we use a singly linked list, we need to remove buckets in the list
			// that have been subject to fusion as mentionned in the comment of add_to_bucket()
			nextBucket = graphArray[bucket].next_bucket;
			while (nextBucket < graphSize && graphArray[nextBucket].component != nextBucket) {
				nextBucket = graphArray[nextBucket].next_bucket;
			}
			graphArray[bucket].next_bucket = nextBucket;

			//mask to short because SIMD type is short (no consequences because this routine checks for equallity)
			int bucketdiff = (*db_differences_fptr)((short*)graphArray[i].db, (short*)graphArray[bucket].db, graphArray[i].solLen, graphArray[bucket].radius + 4);
			//enter the bucket ?
			if (bucketdiff <= graphArray[bucket].radius + 3) {
				if (bucketdiff <= 3) { //same component: add to bucket
					add_to_bucket(i, bucket); //could involve fusion
				}
				else { //scan bucket
					int elt = graphArray[bucket].sibling;
					while (elt != -1) {
						//mask to short because SIMD type is short (no consequences because this routine checks for equallity)
						int eltdiff = (*db_differences_fptr)((short*)graphArray[i].db, (short*)graphArray[elt].db, graphArray[i].solLen, 4);
						if (eltdiff <= 3) { //found a friend
							add_to_bucket(i, bucket);
							elt = -1;
						}
						else { //found no friend yet
							elt = graphArray[elt].sibling;
						}
					}
				}
			}
			//done with bucket, get next bucket
			bucket = nextBucket;
		}
	}
}

/*
 * L^2 algo and type 1 accepted connections between 2D structures.
 * added in f26
 */
int MAX_ERRORS_ALLOWED = 1; //controled by user parameter -MAX_ERRORS
enum Flash_DERR_TYPE { Flash_NO_ERROR = 0, BROKEN_BASE_PAIR = 1, BULGE_MIGRATION_1, BULGE_MIGRATION_2, LONG_BULGE_MIGRATION_1, LONG_BULGE_MIGRATION_2 };
int sameComponent_type1(int a, int b) {

	unsigned short * db1 = graphArray[a].db;
	unsigned short * db2 = graphArray[b].db;
	int n = graphArray[a].solLen;

	static char *oked = NULL;
	static int oked_capacity = 0;
	if (oked == NULL || oked_capacity < n) {
		oked = (char*)realloc(oked, sizeof(char)*n);
		oked_capacity = n;
	}

	//memset(oked, 0, sizeof(char)*n);
	memset(oked, 0, n);

	int max_errors_reached = 0;

	for (int i = 0; i < n; ++i) {
		if (!oked[i]) {
			if (db1[i] != db2[i]) {
				if (max_errors_reached >= MAX_ERRORS_ALLOWED) { return 0; }
				else { ++max_errors_reached; }

				int errType = Flash_NO_ERROR;
				//swap so that if dbx[i] is upaired then x==1
				//swap so that otherwise if dbx[i]<dby[i] then x==2
				if (db2[i] == i) {//|| db1[i]<db2[i]){
					unsigned short *swap = db1;
					db1 = db2;
					db2 = swap;
				}
				if (db1[i] == i && db1[db2[i]] == db2[i]) { //broken base pair?
					errType = BROKEN_BASE_PAIR;
					oked[db2[i]] = 1;
				}
				if (errType == Flash_NO_ERROR) {
					if (db1[i + 1] == db2[i] && db1[i] == i && db2[i + 1] == i + 1) { //this strand bulge migration
						errType = BULGE_MIGRATION_1;
						oked[i + 1] = 1;
						oked[db2[i]] = 1;
					}
					if (errType == Flash_NO_ERROR) {
						//could this be a long migration?
						int lm = db1[db2[i]] == db2[i] && db2[db1[i]] == db1[i];
						if (lm) for (int k = db2[i] + 1; k < db1[i]; ++k) {
							if (!(db1[k] == db2[k] == k)) {
								lm = 0;
								break;
							}
						}
						if (lm) {
							errType = LONG_BULGE_MIGRATION_1;
							oked[db1[i]] = 1;
							oked[db2[i]] = 1;
						}
					}
				}
				if (errType == Flash_NO_ERROR) {
					if (db1[i] < db2[i]) {
						unsigned short *swap = db1;
						db1 = db2;
						db2 = swap;
					}
					if (db1[i] == db2[i] + 1) {//db1[db1[i]-1] == db1[i]-1 && db2[db2[i]+1] == db2[i]+1 ) { //oposite strand bulge migration
						errType = BULGE_MIGRATION_2;
						oked[db2[i]] = 1;
						oked[db1[i]] = 1;
					}
					if (errType == Flash_NO_ERROR) {
						//could this be a long migration?
						int lm = db1[db2[i]] == db2[i] && db2[db1[i]] == db1[i];
						if (lm) for (int k = db2[i] + 1; k < db1[i]; ++k) {
							if (!(db1[k] == db2[k] == k)) {
								lm = 0;
								break;
							}
						}
						if (lm) {
							errType = LONG_BULGE_MIGRATION_2;
							oked[db1[i]] = 1;
							oked[db2[i]] = 1;
						}
					}
				}
				if (errType == Flash_NO_ERROR && ++max_errors_reached >= MAX_ERRORS_ALLOWED) {//other differences
					return 0;
				}
			}
		}
	}
	return 1;
}


void disjointComponents_type1() {
	//int tofuse[graphSize];
	int *tofuse = alloc_aligned(sizeof(int) * graphSize);
	for (int i = 1; i < graphSize; ++i) {
		int p = 0;
		tofuse[p++] = graphArray[i].component;
		for (int j = 0; j < i; ++j) {
			int jcomp = graphArray[j].component;
			int marked = 0;
			for (int k = 0; k < p; ++k) { //already in same component?
				if (jcomp == tofuse[k]) {
					marked = 1;
					break;
				}
			}
			if (!marked) {//i and j are not in the same component
				if (sameComponent_type1(i, j)) { // should they ?
					tofuse[p++] = jcomp; //Let's fuse them.
				}
			}
		}
		//go fuse components
		if (p > 1) {
			//qsort(tofuse,p,sizeof(int),&compare_int_smaller);
			int destination_component = tofuse[1];
			for (int j = 0; j <= i; ++j) {
				for (int k = 0; k < p; ++k) {
					if (graphArray[j].component == tofuse[k]) {
						graphArray[j].component = destination_component;
						break;
					}
				}
			}
		}
	}
	free(tofuse);
}

//1st sorting pass (no graph considerations)
//just sort according to energy score.
int cmpGraphByScoreOnly(const void *one, const void *two) {
	struct graph_t *sdb1 = (struct graph_t *)one;
	struct graph_t *sdb2 = (struct graph_t *)two;

	float score1 = sdb1->score;
	float score2 = sdb2->score;
	if (score1 < score2) return -1;
	if (score1 > score2) return +1;
	return 0;
}
//2nd sorting pass (1-component, 2-score)
int cmpGraphDB(const void *one, const void *two) { //usable by qsort -1 if scored_db1 < scored_db2, ...
	struct graph_t *sdb1 = (struct graph_t *)one;
	struct graph_t *sdb2 = (struct graph_t *)two;

	float score1 = sdb1->score;
	float score2 = sdb2->score;
	int component1 = sdb1->component;
	int component2 = sdb2->component;

	if (component1 < component2) return -1;
	if (component1 > component2) return +1;
	if (score1 < score2) return -1;
	if (score1 > score2) return +1;
	return 0;
}

//convert a structure as ints(internal representation of db) with char dot brackets for output
char *intdb2chardb(char *buf, unsigned short *idb, int len, int *seqAsInt) {
	if (buf == NULL) buf = (char*)alloc_aligned(sizeof(char)*(len + 1));
	buf[len] = 0;
	for (int i = 0; i < len; ++i) {
		//        buf[i] = idb[i] == i ? ('.') : ( idb[i]>i?'(':')');
		if (idb[i] == i) {
			buf[i] = '.';
		}
		else if (idb[i] > i) {
			int c = specialBracketsFor(i, idb[i], seqAsInt);
			if (c) {
				buf[i] = '<';
			}
			else {
				buf[i] = '(';
			}
		}
		else {
			int c = specialBracketsFor(idb[i], i, seqAsInt);
			if (c) {
				buf[i] = '>';
			}
			else {
				buf[i] = ')';
			}
		}
	}
	return buf;
}

/*
 * called after all solutions have been collected
 * sort, get components, sort, output
 * if fullOutput then output all solutions after outputing summary
 * if thin output, don't output shape and remove headers
 */
void outputGraphDBSolutions(int fullOutput, int thinOutput, int* seqAsInt, int graph_type) {

	//1- sort by mfe,
	qsort(graphArray, graphSize, sizeof(struct graph_t), &cmpGraphByScoreOnly);
	//float mfe=graphArray[0].score;
	//float threshold = graphArray[graphSize-1].score;

	//1a- reinitialize the values of next_bucket
	for (int i = 0; i < graphSize; ++i) {
		graphArray[i].next_bucket = 1 + i;
		graphArray[i].component = i;
		graphArray[i].sibling = -1;
		graphArray[i].size = 1;
		graphArray[i].last_sibling = -1;
	}

	//2- compute graph with high dimension metric data structure approach
	switch (graph_type) {
	case 0:
		disjointComponents_type0();
		break;
	case 1:
		disjointComponents_type1();
		break;
	default:
		fprintf(stderr, "ERROR: Found a BUG: Unknown type of graph type:%i. In file %s, at line %i.", graph_type, __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}


	//we have to renumber the components before anything else
	//so: 1) sort by score only
	//    2) starting from mfe, update to component number to a decreasing number starting with -1
	//    3) change the values back to a zero based positive numbering scheme
	// at that point, the components are numbered starting with the global mfe and growing
	// so that the mfe of component with number c is more stable than that of another component's mfe whose component number is c+something>0

	int newc = -1;
	for (int i = 0; i < graphSize; ++i) {
		int c = graphArray[i].component;
		if (c > -1) {
			for (int j = i; j < graphSize; ++j) {
				if (graphArray[j].component == c)
					graphArray[j].component = newc;
			}
			--newc;
		}
	}
	for (int i = 0; i < graphSize; ++i) {
		graphArray[i].component = -graphArray[i].component - 1;
	}

	//now re sort again the array but by placing the components first in the sort order
	qsort(graphArray, graphSize, sizeof(struct graph_t), &cmpGraphDB);


	//output graph summary
	if (!thinOutput)
		printf("components summary\n");
	//we output the array or just the mfe of each component labeled with the size of it's component
	char *cdb = NULL;
	int currentComponent = -1;
	int componentSize = 1;
	//int om_matches[ omasksSize ];
	int *om_matches = alloc_aligned(sizeof(int) * omasksSize);
	for (int omi = 0; omi < omasksSize; ++omi)
		om_matches[omi] = 0;

	//float worse_component_mfe = threshold-3>mfe?threshold-threshold-3:threshold;
	//size_t size_yet =0;
	//int done=0;
		//for( size_t i = 0; i < graphSize && graphArray[i].score < worse_component_mfe; ++i ){
	for (size_t i = 0; i < graphSize; ++i) {
		if (currentComponent != graphArray[i].component) {

			++COMPONENTS;//benchmarking

			currentComponent = graphArray[i].component;
			if (currentComponent > 0) {
				char *matching_labels = labels_for_matches(om_matches);
				if (thinOutput) {
					printf(" %s %i\n", matching_labels, componentSize);
				}
				else {
					printf(" %s %5i structure(s) in component\n", matching_labels, componentSize);
				}
				free(matching_labels);
				componentSize = 1;
			}

			// if((graphSize*0.9)<=(size_yet)){
			//     done=1;
			//     break;
			// }
			// size_yet+=graphArray[i].size;

			for (int omi = 0; omi < omasksSize; ++omi)
				om_matches[omi] = 0;

			cdb = intdb2chardb(cdb, graphArray[i].db, graphArray[i].solLen, seqAsInt);
			update_matches_vector(cdb, om_matches);

			if (thinOutput) {
				printf("%s %+5.3f ", cdb, graphArray[i].score);
			}
			else {
				char * shape = db2shape(cdb);
				printf("%5i %s %+5.3f %s ", graphArray[i].component, cdb, graphArray[i].score, shape);//,maskLabels?maskLabels:"" );
				free(shape);
			}
		}
		else {
			++componentSize;
			cdb = intdb2chardb(cdb, graphArray[i].db, graphArray[i].solLen, seqAsInt);
			update_matches_vector(cdb, om_matches);
		}
	}
	//if(!done && currentComponent>-1){ //end of last line not caught in the for
	if (currentComponent > -1) { //end of last line not caught in the for
		char *matching_labels = labels_for_matches(om_matches);
		if (thinOutput) {
			printf(" %s %i\n", matching_labels, componentSize);
		}
		else {
			printf(" %s %5i structure(s) in component\n", matching_labels, componentSize);
		}
		free(matching_labels);
	}

	//output all structures, sorted (1) by component and (2) by energy
	if (fullOutput) {
		if (thinOutput) {
			printf("\n");
		}
		else {
			printf("\n\nFull output sorted by component and then by energy.\nComponent, size pro-fusion, dot bracket, energy, shape level 5, matching om\n");
		}
		for (size_t i = 0; i < graphSize; ++i) {
			if (i && graphArray[i].component != graphArray[i - 1].component) {
				printf("\n");
			}
			cdb = intdb2chardb(cdb, graphArray[i].db, graphArray[i].solLen, seqAsInt);
			char * shape = db2shape(cdb);
			char * maskLabels = maskLabelsForOmasks(cdb);
			if (thinOutput) {
				printf("%5i %5i %s %+5.3f %s\n", graphArray[i].component, graphArray[i].size, cdb, graphArray[i].score, maskLabels ? maskLabels : "");
			}
			else {
				printf("%5i %5i %s %+5.3f %s %s\n", graphArray[i].component, graphArray[i].size, cdb, graphArray[i].score, shape, maskLabels ? maskLabels : "");
			}
			free(shape);
			free(maskLabels);
		}
	}

	if (cdb != NULL) free(cdb);
	free(om_matches);
}


/*
 * backtrack() calls this function whenever it finds a solution.
 * capture the solution by reading the contents of the bottom stack.
 * place that structure in the array graphArray for latter analysis.
 */
float graphSolutionOutputer(struct status_t *stack, const int bsp,
	const float score, int I, int J, float threshold, const int * seqAsInt) {
	static const size_t initial_capacity = 64;
	if (graphArray == NULL) {
		graphCapacity = initial_capacity;
		graphArray = (struct graph_t *)alloc_aligned(sizeof(struct graph_t) * graphCapacity);
		if (graphArray == NULL) {
			fprintf(stderr, "ERROR: OUT OF MEMORY (in %s, at line %i).\n", __FILE__, __LINE__);
			exit(ABORT_EXIT);
		}
		graphSize = 0;
	}

	if (graphCapacity <= (graphSize + 1)) {
		graphCapacity = graphCapacity * 2;

		graphArray = (struct graph_t*) realloc((void*)graphArray, sizeof(struct graph_t) * graphCapacity);

		if (graphArray == NULL) {
			fprintf(stderr, "ERROR: OUT OF MEMORY (in %s, at line %i).\n", __FILE__, __LINE__);
			exit(ABORT_EXIT);
		}
	}

	//Get length of dotb to output and start offset (I). I > 0 if scanning.
	int solLen = J - I + 1;
	int endpadd = div(solLen, 8).rem;
	unsigned short *dotb = (unsigned short*)alloc_aligned(sizeof(unsigned short*) * (solLen + endpadd));
	if (NULL == dotb) {
		fprintf(stderr, "ERROR: OUT OF MEMORY (in %s, at line %i).\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
	for (int i = 0; i < solLen; ++i)
		dotb[i] = i;


	graphArray[graphSize].score = score;
	graphArray[graphSize].db = dotb;
	graphArray[graphSize].solLen = solLen;
	graphArray[graphSize].component = (int)graphSize;
	graphArray[graphSize].size = 1;

	graphArray[graphSize].radius = 0;
	graphArray[graphSize].next_bucket = (int)graphSize + 1;
	graphArray[graphSize].sibling = -1;
	++graphSize;

	//go down the stack, setting the dot brackets to dotb
	for (int sp = bsp; sp >= 0; --sp) {
		switch (stack[sp].type) {
		case RIGHT: //don't care
			break;
		case LEFT:  //don't care
			break;
		case K:     //don't care
			break;
		case KK:    //pair (i,j)
			dotb[stack[sp].i - I] = stack[sp].j - I;
			dotb[stack[sp].j - I] = stack[sp].i - I;
			//if k internal then (p,q) also paired
			// k is NCMs+1 if terminal (you know, to prevent backtracking... see case KK in backtrack
			if (stack[sp].k < NCMs && stack[sp].k >= FIRST_TWO_STRANDS_NCM) {
				int p = stack[sp].i + (a2d(strands, stack[sp].k, 0) - 1); //the 5' arm of NCM k,
				int q = stack[sp].j - (a2d(strands, stack[sp].k, 1) - 1); //the 3' arm of NCM k.
				dotb[p - I] = q - I;
				dotb[q - I] = p - I;
			}
			break;
		case BRANCH: //don't care
			break;
		case TLRIGHT: //don't care
			break;
		default:    //report weird error
			fprintf(stderr, "ERROR: Found a BUG (inexistant node type in graphSolutionOutputer() ) at line %i in file %s.\n", __LINE__, __FILE__);
			exit(ABORT_EXIT);
			break;
		}
	}


	//insert the new structure in the graph
	//int newc = (int)graphSize-1;

	//this is the L^2 solution :

	/* for( int c = 0; c < (int)graphSize-1; ++c ){ */
	/*     if( 1 == sameComponent(c, newc) ){ */
	/*         ++graphArray[c].degree; */
	/*         ++graphArray[newc].degree; */
	/*         if( graphArray[c].component != graphArray[newc].component ) { //make sure that both entries are not already in same component */
	/*             int gc = graphArray[c].component; */
	/*             int nc = graphArray[newc].component; */
	/*             if( nc == graphSize-1 ){ */
	/*                 graphArray[newc].component = gc; */
	/*             } else { */
	/*                 //fuse both components to smallest of the two */
	/*                 int minc = gc<nc?gc:nc; */
	/*                 int maxc = gc>nc?gc:nc; */
	/*                 for( int i=0;i<graphSize;++i ){ //this could be tuned to use a pfor (sse or something) */
	/*                     if( graphArray[i].component == maxc ) */
	/*                         graphArray[i].component = minc; */
	/*                 } */
	/*             } */
	/*         } */
	/*     } */
	/* } */

	return threshold;
}



/////////////////////////////////////////////////////////////////////////////////////////////// COUNT SOLUTIONS OUTPUTER
int numberOfSolutions = 0;
/*
 * backtrack() calls this function whenever it finds a solution,
 * just keep count of solutions.
 */
float countSolutionOutputer(struct status_t *stack, const int bsp, const float score, int I, int J, float threshold, const int * seqAsInt) {
	++numberOfSolutions;
	return threshold;
}


/////////////////////////////////////////////////////////////////////////////////////////////// DEBUG SOLUTIONS OUTPUTER (NOT USEFUL AS IS)
/*
 * Debuging solution outputer
 * Output it followed by its free folding energy
 * I and J are respectively the first and last nucleotide values (normally this is I=0, J=N-1)
 * If a score higher in stack is smaller than one lower in stack, print the whole stack + dot bracket soln else print nothing.
 */
float debugSolutionOutputer(struct status_t *stack, const int bsp, const float score, int I, int J, float threshold, const int * seqAsInt) {

	//Get length of dotb to output and start offset (I). I > 0 if scanning.
	int solLen = J - I + 1;
	//char dotb[ solLen +1 ]; //build the dot bracket solution in here
	char *dotb = alloc_aligned(sizeof(char) * (solLen + 1)); //build the dot bracket solution in here
	for (int i = 0; i < solLen; ++i)
		dotb[i] = '.';
	dotb[solLen] = 0;


	float bigger = (float)DBL_MAX;
	int sp = bsp;

	int error = 0;
	for (; sp >= 0; --sp) {
		float smaller = stack[sp].score;
		if (smaller != 0) {
			if (smaller > (iota + bigger)) {
				//if(bigger-score>iota){
				++error;
				printf("%5.3f sp=>%i. score(%5.3f -> %5.3f) type(%s -> %s) subtype(%i -> %i), s(%i -> %i), ij(->(%i,%i)),  big-final(%5.3f)\t", bigger - smaller, sp, bigger, smaller, IDs[stack[sp].type], IDs[stack[sp + 1].type], stack[sp].subtype, stack[sp + 1].subtype, stack[sp].s, stack[sp + 1].s, stack[sp + 1].i, stack[sp + 1].j, bigger - score);

				if (1) {
					if ((sp + 1) < bsp) {
						showStatus(stack, sp + 2, 0);
					}
					if (sp < bsp) {
						showStatus(stack, sp + 1, 0);
					}
					showStatus(stack, sp, 0);
					if (sp > 0) {
						showStatus(stack, sp - 1, 0);
					}
					printf("\n");
				}
			}
		}
		bigger = smaller;
		//}
	}

	if (0)if (error) {
		//go down the stack, setting the dot brackets to dotb
		for (int sp = bsp; sp >= 0; --sp) {
			//showStatus(stack,sp,0);
			switch (stack[sp].type) {
			case RIGHT: //don't care
				break;
			case LEFT:  //don't care
				break;
			case K:     //don't care
				break;
			case KK:    //pair (i,j)
				dotb[stack[sp].i - I] = '(';
				dotb[stack[sp].j - I] = ')';
				//if k internal then (p,q) also paired
				if (stack[sp].k >= FIRST_TWO_STRANDS_NCM) {
					int p = stack[sp].i + (a2d(strands, stack[sp].k, 0) - 1); //the 5' arm of NCM k,
					int q = stack[sp].j - (a2d(strands, stack[sp].k, 1) - 1); //the 3' arm of NCM k.
					dotb[p - I] = '(';
					dotb[q - I] = ')';
				}
				break;
			case BRANCH: //don't care
				break;
			case TLRIGHT://don't care
				break;
			default:    //report weird error
				fprintf(stderr, "ERROR: Found a BUG ( inexistant node type in debugSolutionOutputer() ) at line %i in file %s.\n", __LINE__, __FILE__);
				exit(ABORT_EXIT);
				break;
			}
			//printf( "%s %6.4f\n", dotb, score );
		}

		printf("  %s %6.4f\n", dotb, score);
	}
	free(dotb);
	return threshold;
}


/////////////////////////////////////////////////////////////////////////////////////////////// HEAP SOLUTIONS OUTPUTER (LOCATE PROPER THRESHOLD)
/*
 * captures scores of solutions whose scores are better than threshold and modifies
 * threshold if a sufficient number of solutions of better score have already been encountered.
 * We just keep the energies of the solutions, not the solutions themselves.
 *
 * the capacity of the heap is found at heap[0]
 * the size of the heap is found at heap[1]
 * the max element of the heap is found at heap[2]
 */
float * heap = NULL;
void heapSetup(int numberOfSolutions) {
	heap = allocate_heap(numberOfSolutions);
}
void freeHeap() {
	if (heap != NULL) {
		free(heap);
	}
	heap = NULL;
}
float heapSolutionOutputer(struct status_t *stack, const int bsp,
	const float score, int I, int J, float threshold, const int * seqAsInt) {
	float newt = 0;
	add_to_heap(heap, score);
	if (heap[0] == heap[1]) { //heap full?
		newt = heap[2];
	}
	else {
		newt = threshold;
	}
	return newt;
}


float estimateLargerThreshold(float current_threshold, float mfe, const int N) {

	//added in f34
	//At this energy depth it is very unlikely that if the last pass did not add new structures
	//that increasing the THRESHOLD would help.
	const int MAX_EMPTY_THRESHOLD = 1000;

	float sz = heap[1];
	float cap = heap[0];

	static float previous_size = -1;
	if (previous_size == sz && current_threshold > MAX_EMPTY_THRESHOLD) {
		//no solutions between threshold steps
		heap[0] = heap[1]; //lie and make as if the heap had always been of the proper size to get all solutions
		int FOUND_SOLUTIONS = (int)heap[1];
		int REQUIRED_SOLUTIONS = (int)cap;
		fprintf(stderr, "WARNING: Exausted structure pool, found %i solutions and not %i.\n", FOUND_SOLUTIONS, REQUIRED_SOLUTIONS);
		return current_threshold;
	}
	previous_size = sz;

	if (NULL == heap) {
		fprintf(stderr, "ERROR: Found a BUG. Call to estimateLargerThreshold with NULL heap in file %s, at line %i.\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}


	if (sz == cap) {
		fprintf(stderr, "ERROR: Found a BUG. Call to estimateLargerThreshold with full heap in file %s, at line %i.\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}

	if (sz < 300) {
		//estimate is likely very bad
		float ct = current_threshold;

		float tentative_t = 50.0 / N;

		if (ct < tentative_t)
			ct = tentative_t;

		//fprintf(stderr, "ct(%5.3f), sz(%5.3f), mfe(%5.3f)\n",
		//		ct, sz, mfe );
		return (2 * ct) + mfe;
	}

	//find adequate step size
	float step_size = 0.005;
	size_t cnt = 0;
	do {
		cnt = 0;
		step_size *= 2.0;
		for (size_t i = 2; i < (sz + 2); ++i) {
			if (heap[i] >= heap[2] - step_size) {
				++cnt;
			}
		}
	} while (cnt < sz * 0.1); //fraction of the population selected
	if (cnt > 0.125 * sz) {
		step_size = current_threshold / 4.0;
	}

	//count items close to current_threshold
	float above_01 = 0;
	float above_02 = 0;
	for (size_t i = 2; i < (sz + 2); ++i) {
		if (heap[i] >= (heap[2] - step_size)) {
			above_01 += 1.0;
		}
		if (heap[i] >= (heap[2] - (2 * step_size))) {
			above_02 += 1.0;
		}
	}

	double previous_step = (above_02 - above_01);
	double dd = above_01 / previous_step;


	previous_step = above_01;
	float toget = (cap - sz);
	float steps = 0;

	while (toget > 0) {
		float step = dd * previous_step;
		if (step < 10 + iota) {
			break;
		}
		toget -= step;
		previous_step = step;
		steps += 1.0;
	}
	// this is f35 bug fix
	if (steps * step_size < iota) {
		return current_threshold;
	}

	return mfe + current_threshold + (steps * step_size);
}

/////////////////////////////////////////////////////////////////////////////////////////////// BACKTRACK

//returns mfe
float backtrack(const int *seqAsInt, const int N, const float delta_threshold,
	float(*outputSolution)(struct status_t *, const int, const float, const int, const int, float, const int*),
	const struct array_t F, const struct array_t L, const struct array_t R, const struct array_t MB,
	int verbose, float calculated_threshold, const char *seq) {

	//------------------------------ INIT

	//benchmarking
	SOLUTIONS = 0;
	COMPONENTS = 0;
	BACKTRACK_NODES = 0;

	int I = 0, J = N - 1;

	//Get MFE value appropriate for the MODE
	float mfe = a2dh(MB, I, J);
	if (DUPLEX_MODE == SIMPLE) {
		mfe = a2dh(R, I, J);
	}
	if (DUPLEX_MODE == NCM) {
		mfe = a2dh(L, I, J);
	}
	if (DUPLEX_MODE == NCM4) {
		mfe = LARGEFLOAT;
		float thisHinge = a2d(hinges, seqAsInt[I] - 1, seqAsInt[J] - 1);
		float nucleotideTax = reactivity[I] + reactivity[J];
		thisHinge += nucleotideTax;
		for (int k = 0; k < NCMs; ++k) {
			float newScore = a3dh(F, I, J, k) + thisHinge;
			mfe = mfe < newScore ? mfe : newScore;
		}
	}

	//adjust threshold to account for numerical errors
	float precAdj = -mfe * iota;
	precAdj = precAdj < iota ? iota : precAdj;
	precAdj += N * N*iota;
	float threshold = mfe + precAdj + delta_threshold;

	if (calculated_threshold < LARGEFLOAT) {
		threshold = calculated_threshold + precAdj;
	}
	else if (DUPLEX_MODE != NONE || N < 50) {
		//fix threshold for duplex mode and otherwise short sequences (will need to adjust this)
		float delta = 1e-4;
		threshold = threshold < (mfe + delta) ? (mfe + delta) : threshold;
	}

	float delta = 1e-3;
	threshold = threshold < (mfe + delta) ? (mfe + delta) : threshold;


	if (verbose) {
		if (DUPLEX_MODE > 0) {
			fprintf(stdout, "mfe(%5.5f), th(%5.5f) (INFO: Duplex mode solutions are not sorted.)\n", mfe, threshold);
		}
		else if (outputSolution != &heapSolutionOutputer) {
			if (verbose > 1) {
				fprintf(stdout, "mfe(%5.5f), th(%5.5f)\n", mfe, threshold);
			}
		}
	}

	if (NONE != DUPLEX_MODE) {
		if (mfe >= LARGEFLOAT) {
			printf("No solution.\n");
			exit(NORMAL_EXIT);
		}
	}

	if (DUPLEX_MODE != NONE && verbose == 0) {
		printf("%+5.3f", mfe);
		exit(NORMAL_EXIT);
	}


	size_t stackSize = 2 * N;// Amazingly N suffices when ltl is off
	struct status_t *S = (struct status_t*) alloc_aligned(sizeof(struct status_t) * stackSize);
	if (NULL == S) {
		fprintf(stderr, "ERROR: Found a BUG. Can't allocate stack of statuses in: %s, at: %i\n", __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}

	int bsp = 0, //bottom stack pointer
		tsp = (int)stackSize; //top stack pointer

		//start node
	setStatus(S + bsp, I, J, I, J, -1, -1, -1, BRANCH, -1, NORMALBEHAVIOUR, 0, 0);
	if (DUPLEX_MODE > SIMPLE) {
		setStatus(S + bsp, I, J, I, J, -1, -1, -1, LEFT, -1, NORMALBEHAVIOUR, 0, 0);
	}
	S[bsp].forward = 0;
	S[bsp].score = 0;

	if (mfe > LARGEFLOAT - 1000.0) {
		fprintf(stdout, "no solutions.\n");
		exit(NORMAL_EXIT);
	}

	//------------------------------ GO BACKTRACK

	while (bsp >= 0) { //stack not empty

		++BACKTRACK_NODES; //benchmarking

		if (tsp <= bsp) {
			fprintf(stderr, "ERROR: Found a BUG: Crossing stacks.\n");
			exit(ABORT_EXIT);
		}

		//showStatus(S, bsp, tsp);
		//dotBracketSolutionOutputer(S,bsp, S[bsp].score,I,J,threshold, seqAsInt);

		switch (S[bsp].type) {

		case RIGHT:                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case RIGHT

			if (S[bsp].subtype != BLOCKED_STATE && S[bsp].j - S[bsp].i > 1) {

				float forward = a2dh(R, S[bsp].i, S[bsp].j);
				float Rscore = forward + S[bsp].reverse + S[bsp].promise;
				if (Rscore <= threshold) {
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! push a new node to LEFT
					setStatus(S + bsp + 1, S[bsp].i, S[bsp].j, S[bsp].i, S[bsp].j, -1, -1, -1,
						LEFT, -1, NORMALBEHAVIOUR,
						S[bsp].reverse, S[bsp].promise);
					S[bsp + 1].score = Rscore;
					S[bsp + 1].forward = forward;
					//ICI
					if (mfree(S[bsp].j)) {
						--S[bsp].j;//inc this state before pushing
					}
					else {
						S[bsp].subtype = BLOCKED_STATE;
					}
					++bsp;
				}
				else {
					//if( 0 && mfree(S[bsp].j) ){ //ici
					//    --S[bsp].j;//inc this state
					//} else {
					pop(S, &bsp, &tsp);
					//}
				}
			}
			else {
				pop(S, &bsp, &tsp);
			}

			break;                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case RIGHT



		case TLRIGHT:                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case TLRIGHT

			if (S[bsp].subtype != BLOCKED_STATE && ALLOW_LONG_TERMINAL_LOOP) {

				float forward = 0.0;
				float score = forward + S[bsp].reverse + S[bsp].promise;
				float reverse = S[bsp].reverse;

				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Terminate stem by a long zero cost loop

				if (tsp == stackSize) { //one solution found

					if (score <= threshold) {
						S[bsp].forward = forward;
						S[bsp].score = score;
						++SOLUTIONS; //benchmarking
						threshold = (*outputSolution)(S, bsp, score, I, J, threshold, seqAsInt);
					}
					pop(S, &bsp, &tsp);
				}
				else {
					if (score + S[bsp].promise <= threshold) {
						copyStatus(S, tsp++, bsp + 1);
						S[bsp + 1].reverse = reverse + forward;
						S[bsp + 1].forward = forward;
						S[bsp + 1].score = score;
						S[bsp].subtype = BLOCKED_STATE;
						++bsp;
					}
					else {
						pop(S, &bsp, &tsp);
					}
				}
			}
			else {
				pop(S, &bsp, &tsp);
			}



			break;                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case TLRIGHT






		case LEFT:                // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> case LEFT

			//go to K if
			if (S[bsp].subtype != BLOCKED_STATE && S[bsp].j - S[bsp].i > 1) {
				float forward = a2dh(L, S[bsp].i, S[bsp].j);
				float Lscore = forward + S[bsp].reverse + S[bsp].promise;
				if (Lscore <= threshold) {
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! push a new node to K with (i,j) check
					setStatus(S + bsp + 1, S[bsp].i, S[bsp].j, S[bsp].i, S[bsp].j, -1, -1, -1,
						K, -1, NORMALBEHAVIOUR,
						S[bsp].reverse, S[bsp].promise);
					S[bsp + 1].score = Lscore;
					S[bsp + 1].forward = forward;
					//ICI
					if (mfree(S[bsp].i)) {
						++S[bsp].i;//inc this state before pushing
					}
					else {
						S[bsp].subtype = BLOCKED_STATE;
					}
					++bsp;
				}
				else {
					// if( 0 && mfree(S[bsp].i) ){ //ici
					//     ++S[bsp].i;//inc this state
					// } else {
					pop(S, &bsp, &tsp);
					//}
				}
			}
			else {
				pop(S, &bsp, &tsp);
			}

			break;                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case LEFT









		case K:                   // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> case K

			if (!mpair(S[bsp].i, S[bsp].j, seqAsInt)) {
				pop(S, &bsp, &tsp);
				break;
			}

			if (S[bsp].k < 0) { //Start looping on k
				if (ALLOW_STEMS_OF_ONE) {
					S[bsp].k = -1;
				}
				else {
					S[bsp].k = FIRST_TWO_STRANDS_NCM - 1;
				}
			}

			if (++S[bsp].k >= NCMs) { //done looping NCMs
				pop(S, &bsp, &tsp);
			}
			else {
				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if threshold obtained, push stack with (i,j,k) in KK state check
				float thisHinge = a2d(hinges, seqAsInt[S[bsp].i] - 1, seqAsInt[S[bsp].j] - 1);
				float nucleotideTax = reactivity[S[bsp].i] + reactivity[S[bsp].j];
				thisHinge += nucleotideTax;
				float forward = a3dh(F, S[bsp].i, S[bsp].j, S[bsp].k);
				float score
					= S[bsp].reverse
					+ S[bsp].promise
					+ forward
					+ thisHinge;

				if (isnan(score) || !isfinite(score)) score = LARGEFLOAT;
				if (score <= threshold) {
					setStatus(S + bsp + 1, S[bsp].i, S[bsp].j, S[bsp].i, S[bsp].j, S[bsp].k, -1, -1,
						KK, -1, NORMALBEHAVIOUR,
						S[bsp].reverse + thisHinge, S[bsp].promise);
					S[bsp + 1].score = score;
					S[bsp + 1].forward = forward;
					++bsp;//push stack
				} //end sc >= threshold
			} //end sc >= threshold

			break;                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case K









		case KK:                    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> case KK


			if (S[bsp].k >= NCMs) {
				pop(S, &bsp, &tsp);
				break;
			}
			if (S[bsp].k < FIRST_TWO_STRANDS_NCM &&
				((S[bsp].j - S[bsp].i + 1) == a2d(strands, S[bsp].k, 0))
				) {
				//this is a terminal NCM

				//check that all nucleotides in the terminal loop are allowed to not pair according to user mask
					//if( ! (S[bsp].i+1 < S[bsp].j-1 && mfree_loop(S[bsp].i+1, S[bsp].j-1) ) ){
				if (!S[bsp].i + 1 < S[bsp].j - 1 && !mfree_loop(S[bsp].i + 1, S[bsp].j - 1)) {
					pop(S, &bsp, &tsp); break;
				}

				if (!mfree(S[bsp].i + 1)) {
					pop(S, &bsp, &tsp); break;
				}

				//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				float forward = a3dh(F, S[bsp].i, S[bsp].j, S[bsp].k);
				float reverse = S[bsp].reverse;

				float thisHinge = a2d(hinges, seqAsInt[S[bsp].i] - 1, seqAsInt[S[bsp].j] - 1);
				float nucleotideTax = reactivity[S[bsp].i] + reactivity[S[bsp].j];
				thisHinge += nucleotideTax;
				float score = reverse + forward + S[bsp].promise + thisHinge;

				if (tsp == stackSize) { //one solution found ?
					if (score <= threshold) { //yes
						S[bsp].forward = forward;
						S[bsp].score = score;
						++SOLUTIONS; //benchmarking
						if (DUPLEX_MODE != NONE) {
							score = reverse;
							SIMPLEDUPLEX_SolutionOutputer(S, bsp, score, I, J, threshold, seqAsInt);
						}
						else {
							threshold = (*outputSolution)(S, bsp, score, I, J, threshold, seqAsInt);
						}
					}
					pop(S, &bsp, &tsp);
				}
				else {
					S[bsp].k = NCMs + 1; //block reccursion in this node
					if (score <= threshold) {
						copyStatus(S, tsp++, bsp + 1);
						S[bsp + 1].reverse = reverse + forward;
						S[bsp + 1].forward = forward;
						S[bsp + 1].score = score;
						++bsp;
					}
					else {
						pop(S, &bsp, &tsp);
					}
				}
				break;
			}

			//k is an internal NCM
			int p = S[bsp].i + (a2d(strands, S[bsp].k, 0) - 1); //the 5' arm of NCM k,
			int q = S[bsp].j - (a2d(strands, S[bsp].k, 1) - 1); //the 3' arm of NCM k.
			//      if( (q-p) >= (MIND) && p>S[bsp].i && q<S[bsp].j
			if ((q - p) > 1 && p > S[bsp].i && q < S[bsp].j
				&& mpair(p, q, seqAsInt)
				&& (p - S[bsp].i <= 0 || mfree_loop(S[bsp].i + 1, p - 1))
				&& (S[bsp].j - q <= 0 || mfree_loop(q + 1, S[bsp].j - 1))
				) { //can this NCM be realized or is is too close do diagonal or refused by user masks?

				if (++S[bsp].kk < NCMs) {

					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! If threshold obtained: push KK.

					//score for stem ellongation
					float reverse = S[bsp].reverse;
					float promise = S[bsp].promise;
					float forward = a3dh(F, p, q, S[bsp].kk);
					float kNCMe = a2d(energies, seq2idx(seqAsInt, S[bsp].i, p, q, S[bsp].j), S[bsp].k);
					float hinge = a2d(hinges, seqAsInt[p] - 1, seqAsInt[q] - 1);
					float nucleotideTax = reactivity[p] + reactivity[q];

					float transition = a4d(transitions, S[bsp].kk, S[bsp].k, seqAsInt[p] - 1, seqAsInt[q] - 1);
					float junction = a2d(junctions, S[bsp].k, S[bsp].kk);
					float join = ((transition + hinge) > (float)1.0 ? (float)1.0 : transition + hinge) + junction; //min is hidden rule 1;

					float kkreverse = kNCMe + join + reverse + nucleotideTax;

					float score = kkreverse + forward + promise;

					if (DUPLEX_MODE != NONE && S[bsp].kk < FIRST_TWO_STRANDS_NCM) {
						score = reverse + hinge + kNCMe + nucleotideTax;
						kkreverse = score;
						if (isnan(score) || !isfinite(score)) score = LARGEFLOAT;
						if (score <= threshold) {
							setStatus(S + bsp + 1, p, q, p, q, S[bsp].kk, -1, -1,
								KK, -1, NORMALBEHAVIOUR,
								kkreverse, S[bsp].promise);
							S[bsp + 1].forward = forward;
							S[bsp + 1].score = score;
							++bsp;
						}
					}
					else {
						if (isnan(score) || !isfinite(score)) score = LARGEFLOAT;
						if (score <= threshold) {

							setStatus(S + bsp + 1, p, q, p, q, S[bsp].kk, -1, -1,
								KK, -1, NORMALBEHAVIOUR,
								kkreverse, S[bsp].promise);
							S[bsp + 1].forward = forward;
							S[bsp + 1].score = score;
							++bsp;
						}
					}
				}
				else if (DUPLEX_MODE < NCM
					&&  S[bsp].kk == NCMs
					&& (q - p) > 1 //MIN_MB
					&& p + 1 > S[bsp].i
					&& q - 1 < S[bsp].j) {
					//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! push to BRANCH.
					float reverse = S[bsp].reverse;
					float promise = S[bsp].promise;
					float hinge = a2d(hinges, seqAsInt[p] - 1, seqAsInt[q] - 1);//cap the stem
					float nucleotideTax = reactivity[p] + reactivity[q];
					hinge += nucleotideTax;
					float kNCMe = a2d(energies, seq2idx(seqAsInt, S[bsp].i, p, q, S[bsp].j), S[bsp].k);
					float forward = a2dh(MB, p + 1, q - 1);

					float BRANCHreverse = reverse + kNCMe + hinge; //hinge was missing before f5
					if (BRANCHreverse + promise + forward <= threshold) {
						setStatus(S + bsp + 1, p + 1, q - 1, p + 1, q - 1, -1, -1, -1,
							BRANCH, -1, NORMALBEHAVIOUR,
							BRANCHreverse, promise);
						S[bsp + 1].forward = 0;
						S[bsp + 1].score = 0;

						++bsp;
					}
					else { //done cycling on kk and done with branching
						pop(S, &bsp, &tsp);
					}
				}
				else { //done cycling on kk and done with branching
					pop(S, &bsp, &tsp);
				}
			} //if( p-q>=MIND )
			else { //cannot materialize the internal NCM k
				pop(S, &bsp, &tsp);
			}
			break;                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case KK






		case BRANCH:                  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> case BRANCH

			if (S[bsp].s < 0) { //init loop
				S[bsp].s = S[bsp].i + MIND - 1;
				//S[bsp].s = S[bsp].i + 1;
			}
			++S[bsp].s; //increment loop
			if ((S[bsp].s + 1) >= (S[bsp].j - MIND)) { //terminate loop. (a) push a right node (b) push terminal long loop if ok (c) pop
				if (bsp == 0 || S[bsp].behaviour == RIGHTBEHAVIOUR) { //no pairing (i-1,j+1) so no constraints on min length of internal loop
					if (++S[bsp].subtype == 0) {
						setStatus(S + bsp + 1, S[bsp].i, S[bsp].j, S[bsp].i, S[bsp].j, -1, -1, -1,
							RIGHT, -1, NORMALBEHAVIOUR,
							S[bsp].reverse, S[bsp].promise);
						S[bsp + 1].forward = 0;
						S[bsp + 1].score = 0;
						++bsp;
					}
					else if (S[bsp].subtype == 1 && mfree_loop(S[bsp].i, S[bsp].j)) {
						setStatus(S + bsp + 1, S[bsp].i, S[bsp].j, S[bsp].i, S[bsp].j, -1, -1, -1,
							TLRIGHT, -1, NORMALBEHAVIOUR,
							S[bsp].reverse, S[bsp].promise);
						S[bsp + 1].forward = 0;
						S[bsp + 1].score = 0;
						++bsp;
					}
					else {
						pop(S, &bsp, &tsp);
					}
				}
				else {
					// pairing (i-1,j+1): constraints on min internal loop length
					//this node is actually an extension of a stem, don't replicate known NCMs
					if (++S[bsp].subtype == 0) {
						//push 2_6 on a RIGHT node
						if (S[bsp].i + 0 < S[bsp].j - 3 && mfree_loop(S[bsp].j - 3, S[bsp].j)) {
							setStatus(S + bsp + 1, S[bsp].i + 0, S[bsp].j - 4, S[bsp].i + 0, S[bsp].j - 4, -1, -1, -1,
								RIGHT, -1, NORMALBEHAVIOUR, S[bsp].reverse, S[bsp].promise);
							S[bsp + 1].forward = 0;
							S[bsp + 1].score = 0;
							++bsp;
						}
					}
					else if (S[bsp].subtype == 1) {
						//push 6_2 LEFT break
						if (S[bsp].i + 3 < S[bsp].j && mfree_loop(S[bsp].i, S[bsp].i + 3)) {
							setStatus(S + bsp + 1, S[bsp].i + 4, S[bsp].j, S[bsp].i + 4, S[bsp].j, -1, -1, -1,
								LEFT, -1, NORMALBEHAVIOUR, S[bsp].reverse, S[bsp].promise);
							S[bsp + 1].forward = 0;
							S[bsp + 1].score = 0;
							++bsp;
						}
					}
					else if (S[bsp].subtype == 2) {
						//push 6_3 LEFT break
						if (S[bsp].i + 3 < S[bsp].j - 1 && mfree_loop(S[bsp].i, S[bsp].i + 3) && mfree(S[bsp].j)) {
							setStatus(S + bsp + 1, S[bsp].i + 4, S[bsp].j - 1, S[bsp].i + 4, S[bsp].j - 1, -1, -1, -1,
								LEFT, -1, NORMALBEHAVIOUR, S[bsp].reverse, S[bsp].promise);
							S[bsp + 1].forward = 0;
							S[bsp + 1].score = 0;
							++bsp;
						}
					}
					else if (S[bsp].subtype == 3) {
						//push 5_4 LEFT break
						if (S[bsp].i + 3 < S[bsp].j - 2 && mfree_loop(S[bsp].i, S[bsp].i + 2) && mfree_loop(S[bsp].j - 1, S[bsp].j)) {
							setStatus(S + bsp + 1, S[bsp].i + 3, S[bsp].j - 2, S[bsp].i + 3, S[bsp].j - 2, -1, -1, -1,
								LEFT, -1, NORMALBEHAVIOUR, S[bsp].reverse, S[bsp].promise);
							S[bsp + 1].forward = 0;
							S[bsp + 1].score = 0;
							++bsp;
						}
					}
					else if (S[bsp].subtype == 4) {
						//push 4_5 LEFT break;
						if (S[bsp].i + 2 < S[bsp].j - 3 && mfree_loop(S[bsp].i, S[bsp].i + 1) && mfree_loop(S[bsp].j - 2, S[bsp].j)) {
							setStatus(S + bsp + 1, S[bsp].i + 2, S[bsp].j - 3, S[bsp].i + 2, S[bsp].j - 3, -1, -1, -1,
								LEFT, -1, NORMALBEHAVIOUR, S[bsp].reverse, S[bsp].promise);
							S[bsp + 1].forward = 0;
							S[bsp + 1].score = 0;
							++bsp;
						}
					}
					else if (S[bsp].subtype == 5) {
						//long terminal loop
						if (S[bsp].j - S[bsp].i + 1 >= 5 && mfree_loop(S[bsp].i, S[bsp].j)) {
							setStatus(S + bsp + 1, S[bsp].i, S[bsp].j, S[bsp].i, S[bsp].j, -1, -1, -1,
								TLRIGHT, -1, NORMALBEHAVIOUR, S[bsp].reverse, S[bsp].promise);
							S[bsp + 1].forward = 0;
							S[bsp + 1].score = 0;
							++bsp;
						}
					}
					else {
						pop(S, &bsp, &tsp);
					}
				}
				break;
			}

			if (DUPLEX_MODE != NONE)
				break;

			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if threshold obtained: push left on bottom, push MB on top

			float LEFT_forward = a2dh(L, S[bsp].i, S[bsp].s);
			float RIGHT_forward = a2dh(MB, S[bsp].s + 1, S[bsp].j);

			float reverse = S[bsp].reverse;
			float promise = S[bsp].promise;

			float score = LEFT_forward + RIGHT_forward + reverse + promise;

			if (isnan(score) || !isfinite(score)) score = LARGEFLOAT;
			if (score <= threshold) {
				//push left on bs
				setStatus(S + bsp + 1, S[bsp].i, S[bsp].s, S[bsp].i, S[bsp].s, -1, -1, -1,
					LEFT, -1, LEFTBEHAVIOUR,
					reverse, promise + RIGHT_forward);
				S[bsp + 1].forward = LEFT_forward;
				S[bsp + 1].score = score;
				//push MB on ts
				setStatus(S + tsp - 1, S[bsp].s + 1, S[bsp].j, S[bsp].s + 1, S[bsp].j, -1, -1, -1,
					BRANCH, -1, RIGHTBEHAVIOUR,
					0, promise);
				S[tsp - 1].forward = RIGHT_forward;
				S[tsp - 1].score = score;
				++bsp;
				--tsp;
			}

			break;                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< case BRANCH




		default:                      // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DEFAULT case

			fprintf(stderr, "ERROR: Found a BUG: ended up in default case during backtrack.\n");
			exit(ABORT_EXIT);

			break;                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< DEFAULT case

		} // switch
	}//while

	return mfe;
}//backtrack()




/////////////////////////////////////////////////////////////////////////////////////////////// FLASHFOLD

/*
 * fold wrapper
 * load energy tables if not already done
 * call forward on 0..N-1 with window = N
 * call backtrack on 0..N-1 with adequate solution outputer(s)
 */

 //if set, then print F,L,RO,R,MBO,MB to stdout after call to forward
int DEBUG_PRINT_FORWARD_TABLES = 0;

void fold(const char *seq,
	float delta_threshold,
	float explore,
	float(*solutionOutputer)(struct status_t *, const int, const float, const int, const int, float, const int*),
	int verbose,
	int sort_output,
	int fullOutput,
	int thinOutput) {

	if (verbose > 2)
		fprintf(stderr, "INFO: reading files\n");

	if (!LOADEDENERGYTABLES) { //move this check to loadEnergyTables()
		loadEnergyTables();
	}

	if (verbose > 2)
		fprintf(stderr, "INFO: done reading files\n");

	//length of input sequence
	int N = (int)strlen(seq);

	//convert the input nucleotides [acgu] -> [1234]
	int *seqAsInt = toSeqAsInt(seq);
	if (seqAsInt == 0) {
		fprintf(stderr, "ERROR: Found a BUG: bad sequence asFlash_INTat line %i in file %s.\n", __LINE__, __FILE__);
		exit(ABORT_EXIT);
	}


	struct array_t F, L, R, RO, MB, MBO;

	square3Halloc(&F, N, NCMs);
	square2Halloc(&L, N);
	square2Halloc(&R, N);
	square2Halloc(&RO, N);
	square2Halloc(&MB, N);
	square2Halloc(&MBO, N);



	if (verbose > 2)
		fprintf(stderr, "INFO: forward computing\n");



#ifndef WindowsCompilation
	getrusage(RUSAGE_SELF, &ticks[BEFORE_FORWARD]);
#endif // !WindowsCompilation

	//run forward algo on full sequence without windowing
	forward(seqAsInt, N, N, F, L, R, RO, MB, MBO);
#ifndef WindowsCompilation
	getrusage(RUSAGE_SELF, &ticks[AFTER_FORWARD]);
#endif // !WindowsCompilation

	float mfe = a2dh(MB, 0, N - 1);


	if (explore > 0) {
		explore /= 100.0; //from %
		delta_threshold = -1.0 * explore * mfe;
		//fprintf(stdout, "mfe (%f), exp (%f), Th(%f)\n", mfe, explore, delta_threshold);
	}
	else {
		;//fprintf(stdout, "ELSE: mfe (%f), exp (%f), Th(%f)\n", mfe, explore, delta_threshold);
	}


	if (DEBUG_PRINT_FORWARD_TABLES > 0) {
		printf("\n\nF:\n");
		showUpper3DArray(F, N, NCMs);
		printf("\n\nL:\n");
		showUpper2DArray(L, N);
		printf("\n\nRO:\n");
		showUpper2DArray(RO, N);
		printf("\n\nR:\n");
		showUpper2DArray(R, N);
		printf("\n\nMBO:\n");
		showUpper2DArray(MBO, N);
		printf("\n\nMB:\n");
		showUpper2DArray(MB, N);
		printf("\n");
		print_mask_as_table(N, seqAsInt);
	}

	if (a2dh(MB, 0, N - 1) >= LARGEFLOAT) {
		fprintf(stderr, "INFO: No solution. Is your mask too restrictive?.\n");
		exit(NORMAL_EXIT);
	}

	if (verbose > 2)
		fprintf(stderr, "INFO: backtracking\n");

	//
	// some cases require two separate calls to backtrack
	// others don't
	// in some cases we need to fiddle with the solutionOutputer function
	// This should be re-engineered but for now calls are done on a case by case spaghetti
	//

	if (sort_output) {
		if (solutionOutputer == &dotBracketSolutionOutputer) {
			switch (sort_output) {
			case 1:
				solutionOutputer = &sortedDotBracketSolutionOutputer;
				break;
			case 2:
			case 3:
				solutionOutputer = &graphSolutionOutputer;
				break;
			default:
				fprintf(stderr, "ERROR: Found a BUG. Unknown sort option (%i) in file %s at line %i\n", sort_output, __FILE__, __LINE__);
				exit(ABORT_EXIT);
			}
		}
	}

	//backtrack on full sequence without windowing
#ifndef WindowsCompilation
	getrusage(RUSAGE_SELF, &ticks[BEFORE_BACKTRACK]);
#endif // !WindowsCompilation

	float calculated_threshold = LARGEFLOAT;
	do {
		float mfe = backtrack(seqAsInt, N, delta_threshold, solutionOutputer, F, L, R, MB, verbose, calculated_threshold, seq);
		if ((solutionOutputer == &heapSolutionOutputer) && (heap[1] < heap[0])) {
			calculated_threshold = estimateLargerThreshold(delta_threshold, mfe, N);
			delta_threshold = calculated_threshold - mfe;
			if (verbose > 1) {
				fprintf(stderr, "\nINFO: Found %i solutions.\nTrying with new threshold of %+5.3f (delta=%+5.3f)\n", (int)heap[1], calculated_threshold, calculated_threshold - mfe);
			}
			heap[1] = 0; //reset the heap
			if (calculated_threshold >= LARGEFLOAT)
				break;
		}
		else {
			break;
		}
	} while (1);
#ifndef WindowsCompilation
	getrusage(RUSAGE_SELF, &ticks[AFTER_BACKTRACK]);

	getrusage(RUSAGE_SELF, &ticks[BEFORE_LAST_BACKTRACK]);
#endif // !WindowsCompilation    
	if (solutionOutputer == &heapSolutionOutputer) {
		//do second pass using set threshold and dotbracketOutputer
		float threshold = heap[2];
		switch (sort_output) {
		case 0:
			backtrack(seqAsInt, N, delta_threshold, &dotBracketSolutionOutputer, F, L, R, MB, verbose, threshold, seq);
			break;
		case 1:
			backtrack(seqAsInt, N, delta_threshold, &sortedDotBracketSolutionOutputer, F, L, R, MB, verbose, threshold, seq);
			break;
		case 2:
		case 3:
			backtrack(seqAsInt, N, delta_threshold, &graphSolutionOutputer, F, L, R, MB, verbose, threshold, seq);
			break;
		default:
			fprintf(stderr, "ERROR: Found a BUG: unknown sort option (%i) in file %s at line %i\n", sort_output, __FILE__, __LINE__);
			exit(ABORT_EXIT);
		}
	}
#ifndef WindowsCompilation
	getrusage(RUSAGE_SELF, &ticks[AFTER_LAST_BACKTRACK]);

	getrusage(RUSAGE_SELF, &ticks[BEFORE_OUTPUT]);

#endif // !WindowsCompilation

	// if needed, output solutions.
	switch (sort_output) {
	case 0:
		break;
	case 1:
		outputSortedDBSolutions();
		break;
	case 2:
		outputGraphDBSolutions(fullOutput, thinOutput, seqAsInt, 0);
		break;
	case 3:
		outputGraphDBSolutions(fullOutput, thinOutput, seqAsInt, 1);
		break;
	default:
		fprintf(stderr, "ERROR: Found a BUG: unknown sort option (%i) in file %s at line %i\n", sort_output, __FILE__, __LINE__);
		exit(ABORT_EXIT);
	}
#ifndef WindowsCompilation 
	getrusage(RUSAGE_SELF, &ticks[AFTER_OUTPUT]);
#endif // !WindowsCompilation


	//cleanup outputers allocated data
	free_sorted();
	free_graph();

	//we don't free the energy tables (just in case where we need them again.
	//computed stuff, we trash
	free(seqAsInt);
	free(F.data);
	free(L.data);
	free(R.data);
	free(MB.data);
	free(RO.data);
	free(MBO.data);
}

//////////////////////////////////////////////////////////////////////////// PROGRAM OPTIONS AND MAIN()

//memory allocated bellow this line is not freed.

typedef enum opt_type { Flash_INT = 0, Flash_FLOAT = 1, Flash_CHAR_PTR = 2 } opt_type_t;
typedef union opt_value {
	int   INT_value;
	float FLOAT_value;
	char* CHAR_PTR_value;
} opt_value_t;

//order in this enum must match order of declarations in options array ! (else change the code in main)
#define MAX_PARAMETER_TYPES 50
#define MAX_OPTIONS_PER_PARAMETER 10
enum my_options {
	SEQ = 2, NAME, THRESHOLD, EXPLORE, MAX, VERBOSE, VV, COUNT, NOSORT, GRAPH, GRAPH_SUMMARY, THIN_GRAPH_SUMMARY, EDGE_TYPE, MAX_ERRORS,
	DEBUG_OUTPUTER, SHOW_PARSE, PROFILE_TIMES, TIMES_HEADER, SHOW_TABLES,
	MASK, UMASK, OMASK, OUMASK, ALTDB,
	SDUPLEX, ZIP, ZIPZIP, STEMOFONE, USESIMD, LTL,
	TABLES_PATH, REACTIVITY_MASK, REACTIVITY_PER_POSITION, UNBOUNDEDRM, REACTIVITY_MASK_AS_INTS
};
struct options_t {
	int parameters;                               //number of parameters for this option
	char *names[MAX_OPTIONS_PER_PARAMETER];       //names for option: YOU ARE RESPONSIBLE FOR NOT DUPLICATING NAMES IN OTHER OPTIONS
	opt_value_t **values;                         //values for option for each detected parameter : {{one1,two1,tree1},{one2,two2,tree3},NULL}
	int encounters;                               //increased on each encounter of the argument on the line
	int max_encounters;                           //-1 or parameter set no more than this many times
	opt_type_t types[MAX_OPTIONS_PER_PARAMETER];  //types of options (one INT,FLOAT,CHAR_PTR)
	int mandatory;                                //if != 0 the parsing will verify that option is set
	int incompatible_with[MAX_PARAMETER_TYPES]; //set of incompatible options specified using values from my_options
	const char *description;                      //describe the option to the end user.
} options[] = {
	//help is [0], version  is [1]
	{ 0, {"help","h","HELP","H", NULL }, NULL, 0, -1, { -1 }, 0, { -1 }, "Print short help." },
	{ 0, {"V","version", NULL}, NULL, 0, -1, { -1 }, 0, { -1 }, "Print software version." },

	//the rest is ordered as you like
	{ 1, {"seq","s", NULL},                            NULL, 0,  1, { Flash_CHAR_PTR },           1, { -1 }, "Set sequence data (-s ACTGGACA)." },
	{ 1, {"name", "n", NULL},                          NULL, 0,  1, { Flash_CHAR_PTR },           0, { -1 }, "Set sequence name."},
	{ 1, {"threshold","t", NULL},                      NULL, 0,  1, { Flash_FLOAT },              0, { -1 }, "Maximum energy difference of output structures to MFE (kcal/mol). (-t 0.5)\n NOTE: that MFE is defined as the most stable structure outputed.\n So the MFE structure may be different from native if -m or -um are used." },
	{ 1, {"explore","e", NULL},                        NULL, 0,  1, { Flash_FLOAT },              0, { THRESHOLD, -1 }, "Maximum energy difference of output structures to MFE in % of MFE. (-e 1)\n NOTE: that MFE is defined as the most stable structure outputed.\n So the MFE structure may be different from native if -m or -um are used." },
	{ 1, {"floating_threshold", "ft", NULL},           NULL, 0,  1, { Flash_INT},                0, { SDUPLEX, ZIP, ZIPZIP, -1 }, "Adapt threshold so as to output about this many structures. (-ft 1000)" },
	{ 0, {"verbose","v", NULL},                        NULL, 0, -1, { -1 },                 0, { -1 }, "Print additional information." },
	{ 0, {"v2","vv", NULL },                           NULL, 0, -1, { -1 },                 0, { -1 }, "Print additional information with increased verbosity." },
	{ 0, {"count","c", NULL},                          NULL, 0, -1, { -1 },                     0, { VERBOSE, NOSORT, MAX, GRAPH, GRAPH_SUMMARY, SDUPLEX, ZIP, ZIPZIP, -1 }, "Output only the number of structures that would be printed otherwise (Not duplex modes)." },
	{ 0, {"no_sort","ns", NULL},                       NULL, 0, -1, { -1 },                     0, { MAX, GRAPH, GRAPH_SUMMARY, EDGE_TYPE, SDUPLEX, ZIP, ZIPZIP, -1 }, "Do not sort but output solutions as they are found.\n NOTE that output masks are not matched and shapes are note computed in the no_sort mode.\n In this mode memory use is constant relative to output size and is thus suitable for very large outputs." },
	{ 0, {"graph","g", NULL},                          NULL, 0, -1, { -1 },                     0, { SDUPLEX, ZIP, ZIPZIP, GRAPH_SUMMARY, -1 }, "Calculate the disjoint components." },
	{ 0, {"graph_summary","gs", NULL},                 NULL, 0, -1, { -1 },                     0, { SDUPLEX, ZIP, ZIPZIP, -1 }, "Calculate the disjoint components and output only the summary." },
	{ 0, {"thin_graph_summary","tgs", NULL},           NULL, 0, -1, { -1 },                     0, { GRAPH_SUMMARY, GRAPH, SDUPLEX, ZIP, ZIPZIP, -1 }, "Like -gs but slightly thiner output." },
	{ 1, {"edge_type", "et", NULL},                    NULL, 0,  1, {Flash_INT},                0, { -1 }, "Type of operations to build graph (0 or 1). (-et 1)" },
	{ 1, {"MAX_ERRORS", "xe", NULL},                   NULL, 0,  1, {Flash_INT},                0, { -1 }, "Number of allowed steps between pairs of nodes in the graph. (-xe 10)." },
	{ 0, {"debug_outputer","d", NULL},                 NULL, 0, -1, { -1 },                     0, { -1 }, "Used for debuging purposes." },
	{ 0, {"show_parse", "sp", NULL},                   NULL, 0, -1, { -1 },                     0, { GRAPH_SUMMARY, GRAPH, SDUPLEX, ZIP, ZIPZIP, -1 }, "Print parse before dot bracket. Use only with -no_sort" },
	{ 0, {"times", NULL},                              NULL, 0, -1, { -1 },                     0, { -1 }, "Output CPU seconds consumed for each phase and other statistics." },
	{ 0, {"times-header", NULL},                       NULL, 0, -1, { -1 },                     0, { -1 }, "Header for times and statistics." },
	{ 0, {"show_tables","st", NULL},                   NULL, 0, -1, { -1 },                     0, { -1 }, "Print tables used internally (not normally used)." },
	{ 1, {"mask","m", NULL},                           NULL, 0,  1, { Flash_CHAR_PTR },           0, { -1 }, "Specify a balanced mask (-m 'mask')." },
	{ 1, {"unbalanced_mask","um", NULL},               NULL, 0,  1, { Flash_CHAR_PTR },           0, { -1 }, "Specify an unbalanced mask (-um 'mask')." },
	{ 2, {"output_mask","om", NULL},                   NULL, 0, -1, { Flash_CHAR_PTR, Flash_CHAR_PTR }, 0, { -1 }, "Set a named balanced output mask (-om 'name' 'mask')." },
	{ 2, {"output_unbalanced_mask","oum", NULL},       NULL, 0, -1, { Flash_CHAR_PTR, Flash_CHAR_PTR }, 0, { -1 }, "Set a named unbalanced output mask (-oum 'name' 'mask')." },
	{ 0, {"alternate_dot_brackets","alt", NULL},       NULL, 0, -1, { -1 },                     0, { -1 }, "Use alternate dot brackets format for output." },
	{ 1, {"simple_duplex","sd", NULL},                 NULL, 0,  1, { Flash_CHAR_PTR },           0, { ZIP, ZIPZIP, -1 }, "Consult documentation under SIMPLE DUPLEX MODE." },
	{ 1, {"NCM_DUPLEX","zip_duplex", "zd", NULL},      NULL, 0,  1, { Flash_CHAR_PTR },           0, { ZIPZIP, -1 }, "Consult documentation under NCM ZIP MODE." },
	{ 1, {"NCM4_DUPLEX","zipzip_duplex", "zzd", NULL}, NULL, 0,  1, { Flash_CHAR_PTR },           0, { -1 }, "Consult documentation under NCM4 ZIP MODE." },
	{ 0, {"STEM_OF_ONE", "soo", NULL},                 NULL, 0,  1, { -1 },                     0, { -1 }, "Allow terminal stems of one base pair." },
	{ 0, {"SIMD", NULL},                               NULL, 0,  1, { -1 },                     0, { -1 }, "Use SIMD operations if available in graph computation." },
	{ 0, {"NO_LONG_TERMINAL_LOOP","nltl",NULL},        NULL, 0,  1, { -1 },                     0, { -1 }, "Disallow stems to terminate in long ( >4 nucleotides) loops of unpaired nucleotides."},
	{ 1, {"tables", NULL},                             NULL, 0,  1, { Flash_CHAR_PTR },         0, { -1 }, "Set path for tables directory." },
	{ 1, { "reactivity_mask","rm", NULL },             NULL, 0,  1, { Flash_CHAR_PTR },         0, { REACTIVITY_MASK_AS_INTS, -1 }, "Nucleotides energy taxes mask. Use quotes: ie. \"0 1.2 0\" for a 3 nucleotides sequence." },
	{ 2, { "nucleotide_reactivity","nr", NULL },       NULL, 0, -1, { Flash_INT, Flash_FLOAT }, 0, { -1 }, "Energy tax for single nucleotides. Repeat as needed. (-nr POS VALUE) POS is position in seq, Value is Energy." },
	{ 0, { "UNBOUNDEDRM","ubrm",NULL},                 NULL, 0, -1, { -1 },                     0, { -1 }, "Length of reactivity mask may be different than that of sequence." },
	{ 1, { "reactivity_mask_as_integer","rmi",NULL },  NULL, 0,  1, { Flash_CHAR_PTR },         0, { REACTIVITY_MASK, -1 }, "Nucleotides energy taxes mask as a sequence of integers (one per nucleotide, no quote). Ie: -rmi 0000111222..." },

	{ -1,{ (char*)-1 },                                NULL, 0, -1, { -1 },                     0, { -1 }, "" } //<- this entry is required as it marks the end of the array
};


//ERROR reporting for bad options
void badoption(char*opt) {
	fprintf(stderr, "ERROR: option '%s' not recognized at current location in command line.\nERROR: Try --help.\n", opt);
	exit(ABORT_EXIT);
}

//if full == 0 output only end user useful information
void describe_options(int full) {

	for (int i = 0; options[i].parameters != -1; ++i) {
		if (full)
			printf("\nparamNum (%i), encounters (%i)\n", options[i].parameters, options[i].encounters);
		else
			printf("\n");

		char ** ns = options[i].names;
		//printf("alias:");
		for (int nmidx = 0; NULL != ns[nmidx]; ++nmidx) {
			char *n = ns[nmidx];
			printf("-%s ", n);
		}
		printf("\n");

		if (full)
			for (int p = 0; p < options[i].parameters; ++p) {
				static const char *typeNames[] = { "INT", "FLOAT", "CHAR_PTR" };
				printf(" %s", typeNames[options[i].types[p]]);
			}
		if (full)
			if (options[i].parameters > 0)
				printf("\n");
		if (full)
			if (NULL != options[i].values) {
				for (int v = 0; NULL != options[i].values[v]; ++v) {
					for (int p = 0; p < options[i].parameters; ++p) {
						switch (options[i].types[p]) {
						case Flash_INT:
							printf(" %i", options[i].values[v][p].INT_value);
							break;
						case Flash_FLOAT:
							printf(" %+5.3f", options[i].values[v][p].FLOAT_value);
							break;
						case Flash_CHAR_PTR:
							printf(" %s", options[i].values[v][p].CHAR_PTR_value);
							break;
						default:
							fprintf(stderr, "ERROR: Found a BUG. Unknown options' parameter type (%i,i:%i,p:%i) in file %s at line %i.", options[i].types[p], i, p, __FILE__, __LINE__);
							exit(ABORT_EXIT);
						}
						printf("\n");
					}
				}
				printf("\n");
			}

		if (options[i].mandatory) {
			printf(" This parameter is mandatory.\n");
		}

		printf(" %s\n", options[i].description);
	}
	printf("\n");
}


//order of parameters in input are irrelevant
//duplicate parameters with different values are ok
//'-opt' or '--opt' are the same
void parse_options(int argc, char * const argv[]) {

	for (int arg = 1; arg < argc; ++arg) {

		//tolerate options to start with '-' or '--'
		char *p = argv[arg];
		if (p[0] == '-') { ++p; }
		else { badoption(argv[arg]); }
		if (p[0] == '-') { ++p; }

		//try all options to see if one matches
		int found_match_for_p = 0;
		for (int opt = 0; -1 != (options[opt].parameters); ++opt) {
			int param_name_match = -1;
			int param_name_idx = 0;
			while (NULL != options[opt].names[param_name_idx] &&
				(param_name_match = strcmp(options[opt].names[param_name_idx], p))
				) {
				++param_name_idx;
			}
			if (0 == param_name_match) { //one of the names matches the argv[arg]
				found_match_for_p = 1;
				//help
				if (opt == 0) {
					printf("Parameters for %s:\n", argv[0]);
					describe_options(0);
					exit(NORMAL_EXIT);
				}
				//version
				if (opt == 1) {
#ifndef NOTIMESTAMP	  
					// printf( "This is MCFlashfold named %s,\nfrom file %s compiled on %s\n", argv[0], __FILE__, __TIMESTAMP__ );
					printf("This is MCFlashfold (called as %s) version %s compiled on %s\n", argv[0], SOFTWARE_VERSION, __DATE__" "__TIME__);
#else
					printf("This is MCFlashfold (called as %s) version %s\n", argv[0], SOFTWARE_VERSION);
#endif
					exit(NORMAL_EXIT);
				}
				++options[opt].encounters;
				if ((-1 != options[opt].max_encounters) && (options[opt].max_encounters < options[opt].encounters)) {
					fprintf(stderr, "ERROR: parameter %s can only be set %i time%s.\n", p, options[opt].max_encounters, (options[opt].max_encounters > 1 ? "s" : ""));
					exit(ABORT_EXIT);
				}



				//allocate ram for void* of that will receive the options parameter's values
				options[opt].values = (opt_value_t**)realloc(options[opt].values, sizeof(opt_value_t*)*(options[opt].encounters + 1));
				check_out_of_memory((void*)options[opt].values, __FILE__, __LINE__);
				options[opt].values[options[opt].encounters - 1] = (opt_value_t*)alloc_aligned(sizeof(opt_value_t) * (options[opt].parameters));
				check_out_of_memory((void*)options[opt].values[options[opt].encounters - 1], __FILE__, __LINE__);
				options[opt].values[options[opt].encounters] = NULL;

				for (int param_value_idx = 0; param_value_idx < options[opt].parameters; ++param_value_idx) {
					if (++arg >= argc) {
						fprintf(stderr, "ERROR: missing argument for option %s at command line.\n", options[opt].names[0]);
						exit(ABORT_EXIT);
					} //if missing parameter for option at end of argv

					switch (options[opt].types[param_value_idx]) {
					case Flash_INT:
						errno = 0;
						char *endptr = NULL;
						long val = (int)strtol(argv[arg], &endptr, 10);

						if ((endptr != NULL && 0 != endptr[0]) ||
							((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0))) {
							//unable to parse full value, type are unrecognizable char in argv[arg]
							fprintf(stderr, "ERROR: unable to recognize an INTEGER value (%s) at command line.\n", argv[arg]);
							exit(ABORT_EXIT);
						}

						if (INT_MAX < val || INT_MIN > val) {
							fprintf(stderr, "ERROR: value too large or too small (%li) at command line.\n", val);
							exit(ABORT_EXIT);
						}

						options[opt].values[options[opt].encounters - 1][param_value_idx].INT_value = (int)val;

						break;

					case Flash_FLOAT:

						errno = 0;
						endptr = NULL;
						double dval = strtod(argv[arg], &endptr);

						if ((endptr != NULL && 0 != endptr[0]) ||
							((errno == ERANGE && (dval == DBL_MAX || dval == DBL_MIN)) || (errno != 0 && dval == 0))) {
							//unable to parse full dvalue, type are unrecognizable char in argv[arg]
							fprintf(stderr, "ERROR: unable to recognize an NUMERIC value (%s) at command line.\n", argv[arg]);
							exit(ABORT_EXIT);
						}

						if (FLT_MAX < dval || -FLT_MAX > dval) {
							fprintf(stderr, "ERROR: numeric value too large or too small (%+5.3f) at command line.\n", dval);
							exit(ABORT_EXIT);
						}

						options[opt].values[options[opt].encounters - 1][param_value_idx].FLOAT_value = (float)dval;
						break;

					case Flash_CHAR_PTR:
						options[opt].values[options[opt].encounters - 1][param_value_idx].CHAR_PTR_value = argv[arg];
						break;

					default:
						fprintf(stderr, "ERROR: Found a BUG: trying to parse an unknown options' parameter type (%i) in file %s at line %i.",
							options[opt].types[param_value_idx],
							__FILE__, __LINE__);
						exit(ABORT_EXIT);
						break;
					}//switch on option type

				}//parse each parameter value

				continue; //found a matching parameter continue the while
			} //if find a matching parameter
		} //for options[opt].parameters
		if (0 == found_match_for_p) {
			//found no match for arg
			badoption(argv[arg]);
		}

	} //for each arg


	//check if all mandatory options are set
	int problem = 0;
	for (int opt = 0; -1 != (options[opt].parameters); ++opt) {
		if (0 != options[opt].mandatory && 0 == options[opt].encounters) {
			++problem;
			fprintf(stderr, "ERROR: parameter %s must be set.\n", options[opt].names[0]);
		}
	}

	//check if selected options are compatible
	for (int opt = 0; -1 != (options[opt].parameters); ++opt) {
		if (0 < options[opt].encounters) {
			for (int i = 0; -1 != options[opt].incompatible_with[i]; ++i) {
				int oopt = options[opt].incompatible_with[i];
				if (options[oopt].encounters > 0) {
					fprintf(stderr, "ERROR: Set either option '%s' or '%s' but not both.\n", options[opt].names[0], options[oopt].names[0]);
					++problem;
				}
			}
		}
	}

	if (0 != problem) {
		exit(ABORT_EXIT);
	}

} //function parse_options



//return false if fn does not point to a directory (directly of via a link)
int is_dir(char * fn) {
	struct stat  buf;
	if (stat(fn, &buf) == -1) {
		return 0; //not a valid file name
	}
#ifndef WindowsCompilation
	return S_ISDIR(buf.st_mode);
#else
	return (_S_IFDIR & buf.st_mode);
#endif
}

int main(int argc, char * const argv[]) {


#ifndef WindowsCompilation
	getrusage(RUSAGE_SELF, &ticks[AT_BEGIN]);
#endif // !WindowsCompilation


	parse_options(argc, argv);

	if (options[DEBUG_OUTPUTER].encounters > 0) {
		//if debug output asked then also debug parameters
		describe_options(1);
	}

	if (options[SHOW_PARSE].encounters > 0) {
		PRINT_PARSE_BEFORE_DOT_BRACKET_SOLUTIONS = 1;
	}

	if (options[TIMES_HEADER].encounters > 0) {
		//just output header line for times benchmarking and exit
		report_times_header();
		exit(NORMAL_EXIT);
	}

	float threshold = 0;
	if (options[THRESHOLD].encounters > 0) {
		threshold = options[THRESHOLD].values[0][0].FLOAT_value;
		if (threshold < 0) {
			fprintf(stderr, "ERROR: negative threshold values are meaningless.\n");
			exit(ABORT_EXIT);
		}
	}

	float explore = -1.0;
	if (options[EXPLORE].encounters > 0) {
		explore = options[EXPLORE].values[0][0].FLOAT_value;
		if (explore < 0) {
			fprintf(stderr, "ERROR: negative explore values are meaningless.\n");
			exit(ABORT_EXIT);
		}
	}


	char *seq = strdup(options[SEQ].values[0][0].CHAR_PTR_value);

	char *name = NULL;
	if (options[NAME].encounters > 0)
		name = options[NAME].values[0][0].CHAR_PTR_value;

	char *duplex_seq = NULL;
	if (options[SDUPLEX].encounters > 0)
		duplex_seq = options[SDUPLEX].values[0][0].CHAR_PTR_value;
	if (options[ZIP].encounters > 0)
		duplex_seq = options[ZIP].values[0][0].CHAR_PTR_value;
	if (options[ZIPZIP].encounters > 0)
		duplex_seq = options[ZIPZIP].values[0][0].CHAR_PTR_value;

	int verbose = options[VERBOSE].encounters > 0;
	if (options[VV].encounters > 0) {
		verbose = 2;
	}
	


	ALLOW_STEMS_OF_ONE = options[STEMOFONE].encounters;
	ALLOW_LONG_TERMINAL_LOOP = options[LTL].encounters == 0;

#ifndef ___SIMD_COMPILATION__
	if (options[USESIMD].encounters) {
		fprintf(stderr, "WARNING: this compiled version does not support use of SIMD extensions.\n");
		options[USESIMD].encounters = 0;
	}
#endif
	db_differences_fptr = options[USESIMD].encounters ? db_differences_SIMD : db_differences;

	int sort_output = 1; //if 0 dont sort, if 1 sort. if 2 do the graph thing with type0 metric, if 3 do the graph thing with type1 metric.
	int fullOutput = 0;  //if 1 output all structures when computing the graph if 0 then just the summary
	int thinOutput = 0; //if 1 output no shape and truncate some verbose output during -gs.
	if (options[NOSORT].encounters > 0) {
		sort_output = 0;
	}
	if (options[GRAPH_SUMMARY].encounters > 0 || options[THIN_GRAPH_SUMMARY].encounters > 0) {
		if (strlen(seq) >= USHRT_MAX) {
			fprintf(stderr, "ERROR: For efficiency reasons, graph components calculations are limited to sequences shorter than %i.\n", USHRT_MAX);
			exit(ABORT_EXIT);
		}
		sort_output = 2;
		fullOutput = 0;
		if (options[THIN_GRAPH_SUMMARY].encounters > 0) {
			thinOutput = 1;
		}
	}
	if (options[GRAPH].encounters > 0) {
		if (strlen(seq) >= USHRT_MAX) {
			fprintf(stderr, "ERROR: For efficiency reasons, graph components calculations are limited to sequences shorter than %i.\n", USHRT_MAX);
			exit(ABORT_EXIT);
		}
		sort_output = 2;
		fullOutput = 1;
	}

	if (options[EDGE_TYPE].encounters > 0) {
		if (sort_output != 2) {
			fprintf(stderr, "ERROR: You must use -gs or -g to use -et. Try -help.\n");
			exit(ABORT_EXIT);
		}
		int edgetype = options[EDGE_TYPE].values[0][0].INT_value;
		if (edgetype == 1) {
			sort_output = 3;
		}

		//now set user specified number of allowed errors in computation of graph components.
		if (options[MAX_ERRORS].encounters > 0) {
			int max_errors = options[MAX_ERRORS].values[0][0].INT_value;
			if (max_errors > 0) {
				MAX_ERRORS_ALLOWED = max_errors;
			}
			else {
				fprintf(stderr, "ERROR: Value for max errors (%i) is out of range.\n", max_errors);
				exit(ABORT_EXIT);
			}
		}
	}

	if (MAX_ERRORS_ALLOWED > 1 && sort_output < 3) {
		fprintf(stderr, "WARNING: Value for max errors is not compatible with the graph exploration model.\n");
	}


	float(*solutionOutputer)(struct status_t *, const int, const float, const int, const int, float, const int*) =
		&dotBracketSolutionOutputer;

	if (options[COUNT].encounters > 0) {
		solutionOutputer = &countSolutionOutputer;
	}

	if (options[DEBUG_OUTPUTER].encounters > 0) {
		solutionOutputer = &debugSolutionOutputer;
	}

	if (options[PROFILE_TIMES].encounters > 0) {
		KEEP_TIME = 1;
	}


	if (options[SHOW_TABLES].encounters > 0) {
		DEBUG_PRINT_FORWARD_TABLES = 1;
	}

	int return_one_solution = (threshold == 0.0 && explore <= 0.0);
	if (options[MAX].encounters > 0 || return_one_solution) { //number of elts to return floating_threshold
		int solCount = 1;
		if (options[MAX].encounters > 0)
			solCount = options[MAX].values[0][0].INT_value;
		solutionOutputer = &heapSolutionOutputer;
		heapSetup(solCount);
		check_out_of_memory(heap, __FILE__, __LINE__);
	}

	//set the tables directory path
	if (options[TABLES_PATH].encounters > 0) {
		//user supplies the path, trust it
		TABLES_SOURCE_PATH = options[TABLES_PATH].values[0][0].CHAR_PTR_value;
		if (TABLES_SOURCE_PATH[strlen(TABLES_SOURCE_PATH) - 1] != '/') {
			TABLES_SOURCE_PATH = mergeToNewString(options[TABLES_PATH].values[0][0].CHAR_PTR_value, "/");
		}

	}
	else {
		//user did not supply the path on the command line.
		//check if an environment variable is set and trust it.
		//Otherwise, look in ./tables or in ~/MC-Flashfold/tables or else in /usr/local/MC-Flashfold/tables

		const char * tablesDir = getenv("MCFTABLES");
		if (tablesDir != NULL) {
			TABLES_SOURCE_PATH = strdup(tablesDir);
		}
		else {
			TABLES_SOURCE_PATH = strdup("./tables/");
			if (!is_dir(TABLES_SOURCE_PATH)) {
				free(TABLES_SOURCE_PATH);
				TABLES_SOURCE_PATH = mergeToNewString(getenv("HOME"), "/MC-Flashfold/tables/");
				if (!is_dir(TABLES_SOURCE_PATH)) {
					free(TABLES_SOURCE_PATH);
					TABLES_SOURCE_PATH = strdup("/usr/local/MC-Flashfold/tables/");
					if (!is_dir(TABLES_SOURCE_PATH)) {
						free(TABLES_SOURCE_PATH);
						fprintf(stderr, "Could not find a valid tables directory");
						exit(0);
					}
				}
			}
		}
	}
	if (verbose > 1) {
		fprintf(stderr, "INFO: Using tables directory: %s\n", TABLES_SOURCE_PATH);
	}


	DUPLEX_MODE = NONE;

	if (options[SDUPLEX].encounters > 0) {
		DUPLEX_MODE = SIMPLE;
		duplex_seq = options[SDUPLEX].values[0][0].CHAR_PTR_value;
	}
	if (options[ZIP].encounters > 0) {
		DUPLEX_MODE = NCM;
		duplex_seq = options[ZIP].values[0][0].CHAR_PTR_value;
	}
	if (options[ZIPZIP].encounters > 0) {
		DUPLEX_MODE = NCM4;
		duplex_seq = options[ZIPZIP].values[0][0].CHAR_PTR_value;
	}

	if (DUPLEX_MODE != NONE) {
		solutionOutputer = SIMPLEDUPLEX_SolutionOutputer;
		DUPLEX_HAIRPIN_START = (int)strlen(seq);
		DUPLEX_HAIRPIN_LENGTH = 2; //length of "AA"

		if (DUPLEX_MODE == NCM4) {
			//can't run ZIPZIP mode on sequences of unequal lengths.
			if (strlen(seq) != strlen(duplex_seq)) {
				fprintf(stderr, "ERROR: NCM4 zip mode (-zzd) must be used with sequences of equal lenghts.\n       (Maybe you mean -sd or -zd).\n");
				exit(ABORT_EXIT);
			}
		}

		seq = (char*)realloc(seq, sizeof(char) * (DUPLEX_HAIRPIN_LENGTH + strlen(seq) + strlen(duplex_seq) + 1)); // "AA" and 1 for null
		check_out_of_memory(seq, __FILE__, __LINE__);
		strcat(seq, "AA");
		strcat(seq, duplex_seq);
	}


	reactivity = alloc_aligned(strlen(seq) * sizeof(float));
	for (int i = 0; i < strlen(seq); ++i)
		reactivity[i] = 0;

	//these options are not meant to be used together. Their co-occurence is filtered at parameter options parsing.
	if (options[REACTIVITY_MASK].encounters > 0 || options[REACTIVITY_MASK_AS_INTS].encounters > 0) {

		int got = 0; //count of entries
		if (options[REACTIVITY_MASK].encounters > 0) {
			char * values = options[REACTIVITY_MASK].values[0][0].CHAR_PTR_value;
			got = splitToFloat(values, reactivity, strlen(seq));
		}
		else {
			char * values = options[REACTIVITY_MASK_AS_INTS].values[0][0].CHAR_PTR_value;
			got = strlen(values);
			int parseErrorCount = 0;
			const int MAXPARSEERRORCOUNTTOREPORT = 3;
			for (int i = 0; i < got; ++i) {
				if (values[i] > '9' || values[i] < '0') {
					reactivity[i] = 0.0;
					if (parseErrorCount++ < MAXPARSEERRORCOUNTTOREPORT) {
						fprintf(stderr, "WARNING: The reactivity mask as integer has bad value (%c) at position %i. Reseting to 0.0.\n", values[i], i);
					}
					else if (parseErrorCount == 1 + MAXPARSEERRORCOUNTTOREPORT) {
						fprintf(stderr, "WARNING: There are more errors in the reactivity mask.\n");
					}
				}
				else
					reactivity[i] = (float)(values[i] - '0');
			}
		}
		if (got != strlen(seq)) {
			if (options[UNBOUNDEDRM].encounters == 0) {
				fprintf(stderr, "ERROR: The reactivity mask length is different than that of the sequence (seen %i values while sequence has %i nucleotides).\n", got, strlen(seq));
				fprintf(stderr, "ERROR: Fix mask or use -ubrm\n");
				exit(ABORT_EXIT);
			}
			if (verbose) {
				fprintf(stderr, "WARNING: The reactivity mask length is different than that of the sequence (seen %i values while sequence has %i nucleotides).\n", got, strlen(seq));
			}
		}
	}

	int otherReactCount = options[REACTIVITY_PER_POSITION].encounters;
	while (otherReactCount-- > 0) {
		int posInSeq = options[REACTIVITY_PER_POSITION].values[otherReactCount][0].INT_value;
		posInSeq = posInSeq < 1 ? 1 : posInSeq >= strlen(seq) ? strlen(seq) : posInSeq; //clamp user value base 1
		float energyTaxValue = options[REACTIVITY_PER_POSITION].values[otherReactCount][1].FLOAT_value;
		reactivity[posInSeq - 1] += energyTaxValue;
	}


	int fmlen = 0;
	int umlen = 0;
	FULL_MASK = NULL;
	UNBALANCED_MASK = NULL;
	//if there is a full mask then we must conciliate it with the unbalanced one before parsing either
	//if no balanced one specified by user then create one with only xxxxxx... and place there all unbalaced constraints found in the balanced one if any
	if (options[MASK].encounters > 0 || DUPLEX_MODE != NONE) {
		char *bmask = NULL;
		if (DUPLEX_MODE != NONE && options[MASK].encounters == 0) {
			bmask = strdup(seq);
			for (int i = 0; i < strlen(seq); ++i) {
				bmask[i] = 'x';
			}
		}
		else {
			bmask = strdup(options[MASK].values[0][0].CHAR_PTR_value);
		}
		if (DUPLEX_MODE != NONE) {
			for (int i = DUPLEX_HAIRPIN_START; i < DUPLEX_HAIRPIN_START + DUPLEX_HAIRPIN_LENGTH; ++i) {
				bmask[i] = '.';
			}
			bmask[DUPLEX_HAIRPIN_START - 2] = '(';
			bmask[DUPLEX_HAIRPIN_START + DUPLEX_HAIRPIN_LENGTH + 1] = ')';
		}
		fmlen = (int)strlen(bmask);
		char *umask = NULL;
		if (options[UMASK].encounters > 0) {
			umask = strdup(options[UMASK].values[0][0].CHAR_PTR_value);
			umlen = (int)strlen(umask);
			if (umlen != fmlen) {
				fprintf(stderr, "ERROR: -um and -m masks are not of the same length.\n");
				exit(ABORT_EXIT);
			}
		}
		else {
			umask = strdup(bmask);
			umlen = fmlen;
			for (int i = 0; i < fmlen; ++i)
				umask[i] = 'x';
		}
		if (fmlen == strlen(seq)) //otherwise, the error will be caught later on before call to fold. (with an informative message.)
			conciliate_full_masks((int)strlen(seq), bmask, umask);
		if (verbose) {
			fprintf(stderr, "INFO: using bmask: %s\n", bmask);
			fprintf(stderr, "INFO: using umask: %s\n", umask);
		}

		FULL_MASK = string2balanded_idb(bmask);
		UNBALANCED_MASK = string2UMASK_idb(umask);
		free(bmask);
		free(umask);

	}
	else {
		//if there is only an unbalanced mask and the Duplex_mode is NONE then no conciliation is needed
		if (options[UMASK].encounters > 0) {
			const char *mask = options[UMASK].values[0][0].CHAR_PTR_value;
			umlen = (int)strlen(mask);
			UNBALANCED_MASK = string2UMASK_idb(mask);
		}
	}

	//output user masks
	for (int omoc = 0; omoc < options[OMASK].encounters; ++omoc) {
		char *mask = options[OMASK].values[omoc][1].CHAR_PTR_value;
		char *mask_name = options[OMASK].values[omoc][0].CHAR_PTR_value;
		if (strlen(mask) != strlen(seq)) {
			fprintf(stderr, "ERROR: Balanced output mask (%s) is not of the correct length.\n", mask_name);
			exit(ABORT_EXIT);
		}
		add_output_user_mask(mask_name, mask, OUTPUT_MASK_BALANCED);
	}

	//unbalanced output user masks
	for (int omoc = 0; omoc < options[OUMASK].encounters; ++omoc) {
		char *mask = options[OUMASK].values[omoc][1].CHAR_PTR_value;
		char *mask_name = options[OUMASK].values[omoc][0].CHAR_PTR_value;
		if (strlen(mask) != strlen(seq)) {
			fprintf(stderr, "ERROR: Unbalanced output mask (%s) is not of the correct length.\n", mask_name);
			exit(ABORT_EXIT);
		}
		add_output_user_mask(mask_name, mask, OUTPUT_MASK_UNBALANCED);
	}

	//output structures using alternate symbols for base pairs that are not cannonical. ie: () for GC and AU but <> for others
	if (options[ALTDB].encounters > 0) {
		specialBracketsFor = &iscanon;
		USING_ALTERNATE_OUTPUT_DB = 1; //have to set this so that output mask matching knows how to interpret () in output (either it is a base pair or it is a canonical base pair)
	}



	if (NULL != FULL_MASK) {
		if (fmlen != strlen(seq)) {
			fprintf(stderr, "ERROR: balanced mask length differs from length of sequence.\n");
			if (fmlen > strlen(seq)) {
				fprintf(stderr, "       mask too long by %li.\n", fmlen - strlen(seq));
			}
			else {
				fprintf(stderr, "       mask too short by %li.\n", strlen(seq) - fmlen);
			}
			exit(ABORT_EXIT);
		}
	}
	if (NULL != UNBALANCED_MASK) {
		if (umlen != strlen(seq)) {
			fprintf(stderr, "ERROR: unbalanced mask length differs from length of sequence.\n");
			if (umlen > strlen(seq)) {
				fprintf(stderr, "       mask too long by %li.\n", umlen - strlen(seq));
			}
			else {
				fprintf(stderr, "       mask too short by %li.\n", strlen(seq) - umlen);
			}
			exit(ABORT_EXIT);
		}
	}

	if (verbose) {
		printf("Explored\n>%s\n%s\n", name, seq);
		if (verbose > 1) { fprintf(stderr, "input threshold %+5.5f, explore( %+5.5f)\n", threshold, explore); }
	}

	//benchmarking
	INPUT_LENGTH = (int)strlen(seq);

	fold(seq, threshold, explore, solutionOutputer, verbose, sort_output, fullOutput, thinOutput);

	if (solutionOutputer == &countSolutionOutputer) {
		printf("%i\n", numberOfSolutions);
	}

#ifndef WindowsCompilation
	getrusage(RUSAGE_SELF, &ticks[AT_END]);
#endif // !WindowsCompilation


	if (KEEP_TIME) {
		report_times();
	}


	//getchar();
}
