#include <stdio.h>
#include <stdlib.h>

//all heap memory alligned on 16 bytes boundaries if SIMD else, malloc
void *alloc_aligned( ssize_t sz ){
    void * ptr = malloc(sz);
    if( NULL == ptr ){
        fprintf( stderr, "ERROR: Out of memory.\n");
        exit( 0 );
    }
    return ptr;
}

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
int seq2idx( const int *seqAsInt, const int i, int k, const int j, int l ) {
    int val = 0;
    int posInWord = 1;
    if( l>0) while( l>=j ) { val += seqAsInt[l] * (float)trunc( (float)pow(4,posInWord-1) ); ++posInWord; --l; }
    while( k>=i ) { val += seqAsInt[k] * (float)trunc( (float)pow(4,posInWord-1) ); ++posInWord; --k; }
    if(val<21){
        printf("ERROR in seq2idx ikjl(%i,%i,%i,%i)\n", i,k,j,l);
        volatile int *N=0;
        *N=1;
    }
    return val-20-1;
}

/*
 * convert a RNA as character string to integer string where A->1, C->2, G->3, U->4, T->4
 * crash program if another character found in stream.
 */
int *toSeqAsInt( const char *seq ){
    if( !seq ){
        printf("ERROR: no sequence buffer passed at line %i in %s\n.", __LINE__, __FILE__ );
        exit( 0 );
    }
    size_t seqLen = strlen(seq);
    if(seqLen == 0){
        return 0;
    }
    int *seqAsInt = (int*)alloc_aligned( sizeof(int) * (seqLen+1) );
    if( seqAsInt == NULL ){
        printf("ERROR: OUT OF MEMORY at line %i in %s\n", __LINE__, __FILE__ );
        exit( 0 );
    }
    
    for( int i = 0; i < seqLen; ++i ){
        switch( seq[i] ){
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
                exit( 0 );
        } //switch
    } //for i
    return seqAsInt;
}//toSeqAsInt()

main()
{
  int *sai = toSeqAsInt("GAUCGCG");
  int idx = seq2idx( sai, 0, 3, 4, 6 ) + 1;
  printf("%i\n", idx);
}
