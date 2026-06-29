/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Antisense genertor produces a String generator for antisense with fixed 1-4 and 13-15 positions
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 1.0 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Predicate;
import java.util.Set;
import java.util.HashSet;

public class AntisenseGeneratorOld implements Iterator<String> {
    private static final char[] BASES = {'A', 'C', 'G', 'U'};
    // A cytosine, C, at the first position of the guide strand (g1) of an siRNA can increase its stability and potency
    // by promoting a stronger interaction between the guide and passenger strands, which leads to more efficient
    // incorporation into the RNA-induced silencing complex (RISC).

    // Potent siRNAs frequently contain A at g10, a feature associated with enhanced knockdown (Raynolds et al. 2004 Nat Biotech; Vert et al. 2006 BMC Bioinformatics)

    // Reynolds' and Tuschl's rules noted a general preference for an A or U residue at positions around the cleavage site (e.g., position 10) in potent siRNAs.
    private final char[] currentKmer = { 'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'U', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', 'U', 'U' };
    private final int[] variablePositions;
    private final int[] indices;
    private final Predicate<String> filter;
    private String nextValid = null;
    private boolean computedNext = false;
    private Set<String> seedsToAvoid = new HashSet<>();
    private Set<String> suppsToAvoid = new HashSet<>();

    public AntisenseGeneratorOld( String fixed1to4, String fixed13to15, Predicate<String> filter ) {
	if( fixed1to4.length() != 4 || fixed13to15.length() != 3 ) {
            throw new IllegalArgumentException( "Expected 4 bases for positions 1–4 and 3 bases for positions 13–15" );
        }
	
	this.filter = filter;

        // set fixed bases
        for( int i = 0; i < 4; i++ ) this.currentKmer[ 1 + i] = fixed1to4.charAt( i );
        for( int i = 0; i < 3; i++ ) this.currentKmer[12 + i] = fixed13to15.charAt( i );

	// track variable poisitions
        this.variablePositions = new int[9];
        int v = 0;
        for( int i = 0; i < 21; i++ ) {
            if( ( i > 4 && i < 12 ) || ( i > 14 && i < 17 ) ) {
                variablePositions[v++] = i;
            }
        }

        this.indices = new int[variablePositions.length]; // all start at 0
	advance(); // preload the first valid sequence
    }

    private void advance() {
        this.nextValid = null;
        while( true ) {
	    // if index combinations exhausted, break
	    if( !incrementIndices() ) {
                break;
            }

	    // update the currentKmer
            for( int i = 0; i < variablePositions.length; i++ ) {
                this.currentKmer[variablePositions[i]] = BASES[indices[i]];
            }
	    
	    String candidate = new String( currentKmer );
	    // make sure candidate does not include a seed to be avoided
	    if( !this.seedsToAvoid.contains( candidate.substring( 1, 8 ) ) ) {
		if( this.filter.test( candidate ) ) {
		    this.nextValid = candidate;
		    break;
		}
	    }
        }
        this.computedNext = true;
	//System.out.println( "next: " + this.nextValid );
    }

    private boolean incrementIndices() {
        for( int i = indices.length - 1; i >= 0; i-- ) {
            if( indices[i] < BASES.length - 1 ) {
                indices[i]++;
                for( int j = i + 1; j < indices.length; j++ ) {
                    indices[j] = 0;
                }
                return true;
            }
        }
        return false; // no more combinations
    }

    public void addSeedToAvoid( String seed ) {
	this.seedsToAvoid.add( seed );
    }

    @Override
    public boolean hasNext() {
        if( !this.computedNext ) advance();
	return nextValid != null;
    }

    @Override
    public String next() {
        if( !hasNext() ) throw new NoSuchElementException();
	String result = this.nextValid;
	this.computedNext = false;
	this.advance();
        return result;
    }
}
