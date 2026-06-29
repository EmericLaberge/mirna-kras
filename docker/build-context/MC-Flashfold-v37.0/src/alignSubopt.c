//
//  alignSubopt.c
//
//  Created by Paul Dallaire on 5/14/13.
//  Copyright (c) 2013 Paul Dallaire. All rights reserved.
//
/*
 * alignSubopt.c
 * Performs sequence alignment with affine gaps and outputs sup-optimal alignments 
 * whose scores are not worse than the best score minus some threshold.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <unistd.h>

//strdup has been DEPRECATED FROM standard unix because of security risks
//so we reintroduce a naive version here
char *strdup(const char *str)
{
    ssize_t n = strlen(str) + 1;
    char *dup = malloc(n);
    if( !dup ){
        fprintf(stderr,"Out of memory\n");
        exit(1);
    }
    
    if(dup)
    {
        strcpy(dup, str);
    }
    return dup;
}

//maximum of two floats
float MAXF( float a, float b ){ return a>b?a:b; }

//values for the scoring of alignment
float
DEFAULT_BAD = -1E10, //some small number
DEFAULT_GOOD = 0,
GAP = -1,  //gap ellongation cost
INIT = -5.0, //gap initiation cost
MATCH_VALUE = 2.0,
MISSMATCH_VALUE = -3.0;

//refine arbitrarily as needed.
//we just check for identity of characters
float MATCH(char a, char b){
    if( a==b ) return MATCH_VALUE;
    return MISSMATCH_VALUE;
}


/*
 * User supplied information
 */
char *SA, *SB; //the sequences
int n, m;      //length of the sequences
float threshold; //Distance from BAS to output solutions for (T)
int verbose = 0; //dump some internal data to stdout

float *A = NULL, *M = NULL, *H = NULL, *V = NULL, *X = NULL; //dynamic programming tables

//accessor for the dynamic programming tables
// (as bad practice that simplifies code, we hardcode the width
//  of the tables using the global m wich is the length of the sequence SB )
float *a2d( float* M, int row, int col ){
    return M + ( (row * (m+1)) + col );
}

//fill the dynamic programming table
int forward(){
    
    //alloc matrices for dynamic prog
    //free is at end of main (at least after the call to backtrack())
    A = malloc( (n+1) * (m+1) * sizeof(float) );
    M = malloc( (n+1) * (m+1) * sizeof(float) );
    H = malloc( (n+1) * (m+1) * sizeof(float) );
    V = malloc( (n+1) * (m+1) * sizeof(float) );
    X = malloc( (n+1) * (m+1) * sizeof(float) );

    if( !A || !M || !H || !V || !X ){
        fprintf(stderr,"Out of memory\n");
        exit(1);
    }
    
    //initialize tables
    for( int i=0;i<=n;++i )
        for(int j=0;j<=m;++j ){
            (*a2d(A,i,j)) = DEFAULT_BAD;
            (*a2d(M,i,j)) = DEFAULT_BAD;
            (*a2d(H,i,j)) = DEFAULT_BAD;
            (*a2d(V,i,j)) = DEFAULT_BAD;
            (*a2d(X,i,j)) = DEFAULT_BAD;
        }

    //epsilon cases
    (*a2d(M,0,0))=DEFAULT_GOOD;
    
    (*a2d(V,0,0))=INIT;
    (*a2d(H,0,0))=INIT;
    
    for( int i=1; i<=n; ++i )
        (*a2d(V,i,0)) = (*a2d(V,i-1,0)) + GAP;
    for( int j=1; j<=m; ++j )
        (*a2d(H,0,j)) = (*a2d(H,0,j-1)) + GAP;
    
 
    //main dynamic programming loop
    for( int i=1;i<=n;++i){
        for( int j=1;j<=m;++j){
            
            
            (*a2d(M,i,j)) =
            MAXF( MAXF(
                       MATCH(SA[i-1],SB[j-1]) + (*a2d(M,i-1,j-1)),
                       MATCH(SA[i-1],SB[j-1]) + (*a2d(V,i-1,j-1))
                       ),
                 MATCH(SA[i-1],SB[j-1]) + (*a2d(H,i-1,j-1))
                 );
            
            (*a2d(V,i,j)) =
            MAXF(
                 INIT + GAP + (*a2d(M,i-1,j)),
                 GAP + (*a2d(V,i-1,j))
                 );
            
            (*a2d(H,i,j)) =
            MAXF( MAXF(
                       INIT + GAP + (*a2d(M,i,j-1)),
                       GAP + (*a2d(H,i,j-1))
                       ),
                 GAP + (*a2d(X,i,j-1)) //DEFAULT_BAD
                 );
            
            (*a2d(A,i,j)) =
            MAXF( MAXF(
                       (*a2d(M,i,j)),
                       (*a2d(V,i,j))
                       ),
                 (*a2d(H,i,j))
                 );
            
            (*a2d(X,i,j)) = GAP + (*a2d(V,i-1,j));
                
        }
    }
    
    return 0;
}

/*
 * DEBUG: print matrices to screen after dynamic programming routine.
 */
float cap( float f ){ //limit size of low values so that they are easily printed on screen
    if( f < -99.9 ) return -99.0;
    return f;
}
void showMatrices(){
    printf("A:%s\nB:%s\n",SA,SB);
    printf("table A:\n");
    for(int i=0;i<=n;++i){
        for( int j=0; j<=m; ++j ){
            printf(" %+3.0f", cap(*a2d(A,i,j)));
        }
        printf( "\n" );
    }
    printf("\ntable M:\n");
    for(int i=0;i<=n;++i){
        for( int j=0; j<=m; ++j ){
            printf(" %+3.0f", cap(*a2d(M,i,j)));
        }
        printf( "\n" );
    }
    printf("\ntable V:\n");
    for(int i=0;i<=n;++i){
        for( int j=0; j<=m; ++j ){
            printf(" %+3.0f", cap(*a2d(V,i,j)));
        }
        printf( "\n" );
    }
    printf("\ntable H:\n");
    for(int i=0;i<=n;++i){
        for( int j=0; j<=m; ++j ){
            printf(" %+3.0f", cap(*a2d(H,i,j)));
        }
        printf( "\n" );
    }

    printf("\ntable X:\n");
    for(int i=0;i<=n;++i){
        for( int j=0; j<=m; ++j ){
            printf(" %+3.0f", cap(*a2d(X,i,j)));
        }
        printf( "\n" );
    }
}

/*
 * definition of backtracking nodes and stack
 */
//the following two need to be kept in sync
const char *NODE_NAMES [] = { "START_MATCH", "START_V", "START_H", "EXTENSION", "CLOSING_V", "CLOSING_H", "INITGAP_V", "INITGAP_H", "EXT_V", "EXT_H", "EXT_HX", "EXT_XV", "DONE" };
enum { START_MATCH = 0, START_V, START_H, EXTENSION, CLOSING_V, CLOSING_H, INITGAP_V, INITGAP_H, EXT_V, EXT_H, EXT_HX, EXT_XV, DONE };

//the states contain the following stateful values
typedef struct {
    int i,j;
    float g_score;
    char emit_a, emit_b, match_type;
    int NODE_IDENTIFIER, CNODE;
} node_t;

node_t *stack;

//initialize a stack state using given values
void init_node( int sp, int i, int j, float g_score, int NODE_IDENTIFIER ){
    stack[sp].i = i;
    stack[sp].j = j;
    stack[sp].g_score = g_score;
    stack[sp].emit_a = ' ';
    stack[sp].emit_b = ' ';
    stack[sp].match_type = ' ';
    stack[sp].NODE_IDENTIFIER = NODE_IDENTIFIER;
}

//debug printing of states on the stack
void printnode( int sp ){
    char ea = stack[sp].emit_a;
    char eb = stack[sp].emit_b;
    if( ea == 0 ) ea = '0';
    if( eb == 0 ) eb = '0';
    printf("(%i) %s, (%i,%i) %5.3f [%c %c]\n", sp, NODE_NAMES[stack[sp].CNODE], stack[sp].i, stack[sp].j, stack[sp].g_score, ea, eb );
}

//called everytime the stack contains a complete alignment (ie: i==-1 and j==-1)
void printSolution( int sp ){
    char * solna = malloc( sizeof(char) * n+m+1 );
    char * solnb = malloc( sizeof( char ) * n+m+1);
    char * mtype = malloc( sizeof( char ) * n+m+1 );
    
    if( !solna || !solnb || !mtype ){
        fprintf(stderr,"Out of memory\n");
        exit(1);
    }
    bzero(solna,n+m+1);
    bzero(solnb,n+m+1);
    bzero(mtype,n+m+1);
    
    float g_score = stack[sp].g_score;
    
    int i = 0;
    //the entry at the top of the stack is useless (except for score value), so is the one on the bottom
    //still we output them and they look like spaces.
    //they can be useful in verbose mode
    while( sp >= 0 ){
        if(verbose) printnode(sp);
        solna[i]=stack[sp].emit_a;
        solnb[i]=stack[sp].emit_b;
        mtype[i]=stack[sp].match_type;
        if(solna[i] == 0) solna[i]='-';
        if(solnb[i] == 0) solnb[i]='-';
        ++i;
        --sp;
    }
    printf("%5.3f\n\t%s          \n\t%s\n\t%s\n",g_score, solna, mtype, solnb);
    free(solna);
    free(solnb);
    free(mtype);
}

//enumerate sub-optimal solutions
int backtrack(){
    
    //this much should suffice amply
    const int stackHeight = (m + n) * 2;
    
    //free at end of backtrack()
    stack = malloc( sizeof(node_t) * stackHeight);
    if( !stack ){
        fprintf(stderr,"Out of memory\n");
        exit(1);
    }
    
    int tos = -1; //Top Of Stack
    
    //push start node on stack
    init_node( ++tos, n-1, m-1, 0, START_MATCH);
    
    float BAS = (*a2d(A,n,m)); //Best Acheivable Score
    printf("BAS: %f\n",BAS);
    while( tos >= 0 ){
        assert( tos < stackHeight );
        //printnode(tos);
        
        int i = stack[tos].i;
        int j = stack[tos].j;
        
        float g = stack[tos].g_score;

        stack[tos].CNODE = stack[tos].NODE_IDENTIFIER;
        
        if( i < 0 && j < 0 ){
            //output solution and then pop
            printf("\n");
            if( stack[tos-1].CNODE == EXT_V
               || stack[tos-1].CNODE == EXT_H
               || stack[tos-1].CNODE == EXT_HX
               || stack[tos-1].CNODE == EXT_XV
               )
                stack[tos].g_score += INIT;
            
            printSolution(tos);
            --tos;
        }
        
        switch( stack[tos].NODE_IDENTIFIER ){
            case START_MATCH : {
                float BSP = *(a2d(M,i+1, j+1));
                //when this node is backtracked into, it will apply Rule 2 (A -> G)
                stack[tos].NODE_IDENTIFIER = START_V;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = ' ';
                    stack[tos].emit_b = ' ';
                    stack[tos].match_type = ' ';
                    init_node( tos+1, i, j, 0.0, EXTENSION );
                    //push it
                    ++tos;
                }
                break;
            }
            case START_V :{
                float BSP = *(a2d(V,i+1, j+1));
                // we pop() by setting this node to DONE
                stack[tos].NODE_IDENTIFIER = START_H;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = ' ';
                    stack[tos].emit_b = ' ';
                    stack[tos].match_type = ' ';
                    init_node( tos+1, i, j, 0.0, EXT_V);
                    ++tos;
                } 
                break;
            }

            case START_H :{
                float BSP = *(a2d(H,i+1, j+1));
                // we pop() by setting this node to DONE
                stack[tos].NODE_IDENTIFIER = DONE;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = ' ';
                    stack[tos].emit_b = ' ';
                    stack[tos].match_type = ' ';
                    init_node( tos+1, i, j, 0.0, EXT_H);
                    ++tos;
                }
                break;
            }

                
            case EXTENSION :{
                if( !(i>=0 && j>=0) ){
                    //can't extend.
                    --tos;
                    continue;
                }
                
                float h = (*a2d(M,i-1+1,j-1+1));
                float BSP = g + MATCH(SA[i],SB[j]) + h;
                
                //On backtrack, go to next rule
                stack[tos].NODE_IDENTIFIER = CLOSING_V;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = SA[i];
                    stack[tos].emit_b = SB[j];
                    stack[tos].match_type = SA[i]==SB[j]?'|':'X';
                    init_node( tos+1, i-1, j-1, BSP - h, EXTENSION );
                    ++tos;
                } 
                break;
            }

            case CLOSING_V :{
                
                if( !(i>=0 && j>=0) ){
                    //can't close.
                    --tos;
                    continue;
                }
                
                float h = (*a2d(V,i-1+1,j-1+1));
                float BSP = g + MATCH(SA[i],SB[j]) + h;
                
                //
                stack[tos].NODE_IDENTIFIER = CLOSING_H;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = SA[i];
                    stack[tos].emit_b = SB[j];
                    stack[tos].match_type = SA[i]==SB[j]?'|':'X';
                    init_node( tos+1, i-1, j-1, BSP - h, EXT_V );
                    ++tos;
                }
                break;
            }
                
            case CLOSING_H :{
                
                if( !(i>=0 && j>=0) ){
                    //can't close.
                    --tos;
                    continue;
                }
                
                float h = (*a2d(H,i-1+1,j-1+1));
                float BSP = g + MATCH(SA[i],SB[j]) + h;
                
                //this node should not be backtracked into so we replace this node by the one we push
                //or we pop()
                stack[tos].NODE_IDENTIFIER = DONE;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = SA[i];
                    stack[tos].emit_b = SB[j];
                    stack[tos].match_type = SA[i]==SB[j]?'|':'X';
                    init_node( tos+1, i-1, j-1, BSP - h, EXT_H );
                    ++tos;
                }
                break;
            }
                
            case INITGAP_V:{
                if( !(i>=0) ){
                    //can't extend on A.
                    stack[tos].NODE_IDENTIFIER = DONE;
                    continue;
                }
                                
                float h = (*a2d(M,i-1+1,j+1));
                float BSP = g + INIT + GAP + h;
                
                //this node should not be backtracked into so we replace this node by the one we push
                //or we pop()
                stack[tos].NODE_IDENTIFIER = DONE;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = SA[i];
                    stack[tos].emit_b = 0;
                    stack[tos].match_type = '+';                    
                    init_node( tos+1, i-1, j, BSP - h, EXTENSION );
                    ++tos;
                } 
                break;
            }
            case INITGAP_H:{
                
                if( !(j>=0) ){
                    //can't extend on B.
                    stack[tos].NODE_IDENTIFIER = DONE;
                    continue;
                }
                

                float h = (*a2d(M,i+1,j-1+1));
                float BSP = g + INIT + GAP + h;
                
                //this node should not be backtracked into so we replace this node by the one we push
                //or we pop()
                stack[tos].NODE_IDENTIFIER = EXT_HX;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = 0;
                    stack[tos].emit_b = SB[j];
                    stack[tos].match_type = '-';                 
                    init_node( tos+1, i, j-1, BSP - h, EXTENSION );
                    ++tos;
                }
                break;
            }
                
            case EXT_V :{
                
                if( !(i>=0) ){
                    //can't extend on A.
                    stack[tos].NODE_IDENTIFIER = DONE;
                    continue;
                } 

                float h = (*a2d(V,i-1+1,j+1));
                float BSP = g + GAP + h;
                
                //this node should not be backtracked into so we replace this node by the one we push
                //or we pop()
                stack[tos].NODE_IDENTIFIER = INITGAP_V;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = SA[i];
                    stack[tos].emit_b = 0;
                    stack[tos].match_type = '+';              
                    init_node( tos+1, i-1, j, BSP - h, EXT_V );
                    ++tos;
                }
                break;
            }

            case EXT_H :{
                if( !(j>=0) ){
                    //can't extend on B.
                    stack[tos].NODE_IDENTIFIER = EXT_HX;
                    continue;
                }

                float h = (*a2d(H,i+1,j-1+1));
                float BSP = g + GAP + h;
                
                //this node should not be backtracked into so we replace this node by the one we push
                //or we pop()
                stack[tos].NODE_IDENTIFIER = INITGAP_H;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = 0;
                    stack[tos].emit_b = SB[j];
                    stack[tos].match_type = '-';                 
                    init_node( tos+1, i, j-1, BSP - h, EXT_H );
                    ++tos;
                } 
                break;
            }
                
                
                
            case EXT_HX :{
                
                if( !(i>=0) ){
                    //can't extend on A.
                    stack[tos].NODE_IDENTIFIER = DONE;
                    continue;
                }
                
                float h = (*a2d(X,i+1,j-1+1));
                float BSP = g + GAP + h;
                
                //this node should not be backtracked into so we replace this node by the one we push
                //or we pop()
                stack[tos].NODE_IDENTIFIER = DONE;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = 0;
                    stack[tos].emit_b = SB[j];
                    stack[tos].match_type = '-';
                    init_node( tos+1, i, j-1, BSP - h, EXT_XV );
                    ++tos;
                }
                break;
            }

                //this is exactly the same definition as EXT_V
            case EXT_XV :{
                
                if( !(i>=0) ){
                    //can't extend on A.
                    stack[tos].NODE_IDENTIFIER = DONE;
                    continue;
                }
                
                float h = (*a2d(V,i-1+1,j+1));
                float BSP = g + GAP + h;
                
                //this node should not be backtracked into so we replace this node by the one we push
                //or we pop()
                stack[tos].NODE_IDENTIFIER = DONE;
                if( BSP >= BAS - threshold ){
                    //set next node on stack
                    stack[tos].emit_a = SA[i];
                    stack[tos].emit_b = 0;
                    stack[tos].match_type = '+';
                    init_node( tos+1, i-1, j, BSP - h, EXT_V );
                    ++tos;
                }
                break;
            }


            case DONE :
            {
                --tos;
                break;
            }
            default :
                fprintf(stderr, "Caught undefined NODE_IDENTIFIER during backtracking.\n");
                break;

        } //end switch(stack[tos].TypeOfNode)
        
    }//end while(tos>=0)
    
    free(stack);
    return 0;
}

/*
 * parse programm arguments and set main variables
 * no error checking so watch out
 */
int parseArgs( int argc, const char **argv ){

    if( ( argc == 5 && strcmp(argv[argc-1], "DEBUG") )
       ||
       ( argc != 4 && argc !=5) )
    {
        printf( "Specify threshold followed by both sequences to match on command line.\n");
        printf( "Add keyword \"DEBUG\" at end of line to see dynamic programming tables and backtracking stack of solutions.\n");
        return 1;
    }
    
    if( argc == 5 && !strcmp(argv[argc-1], "DEBUG") )
       verbose = 1;
    
    //set sequences data
    //free at end of main()
    SA = strdup(argv[2]);
    SB = strdup(argv[3]);

    //find values for length of sequences
    n = (int)strlen(SA);
    m = (int)strlen(SB);
    
    threshold = atof(argv[1]);
    
    return 0;
}


int main(int argc, const char * argv[])
{
    if( parseArgs( argc, argv ) != 0 ){
        return(1);
    }
    
    forward();
    if( verbose ) showMatrices();
    backtrack();
    
    //free matrices memory
    free(M);
    free(A);
    free(V);
    free(H);
    free(X);

    //free sequences memory
    free(SA);
    free(SB);
    
    return 0;
}

