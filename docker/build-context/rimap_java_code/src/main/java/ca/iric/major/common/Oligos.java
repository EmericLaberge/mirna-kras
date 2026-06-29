/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.ArrayList;

/**
* Oligos is a class for generating RNA oligos
*
* @license     MIT
* @author      Francois Major, Université de Montréal
* @version     %I%, %G%
* @since       1.0
*/

/** -------------------------------------------
      Oligos
        is a class to produce various RNA oligo-nucleotide constructs
*/

public class Oligos {

    static public char[] nucleotides = { 'A', 'C', 'G', 'U' };
    static public int alphabetLength = 4;

    public static void allKmer( int k, ArrayList<String> kmers ) {
	allKmerRec( nucleotides, "", alphabetLength, k, kmers );
    }

    private static void allKmerRec( char[] set, String prefix, int n, int k, ArrayList<String> kmers ) {
	// Base case: k is 0, print prefix
	if(k == 0) {
	    kmers.add( prefix );
	    return;
	}
	// One by one add all characters
	// from set and recursively
	// call for k equals to k-1
	for( int i = 0; i < n; ++i ) {
	    // Next character of input added
	    String newPrefix = prefix + set[i];
         
	    // k is decreased, because
	    // we have added a new character
	    allKmerRec( set, newPrefix, n, k - 1, kmers );
	}
    }
}
