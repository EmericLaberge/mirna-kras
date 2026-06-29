/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * Class for an RNA sequence based on a String implementation
 *
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2021 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

public class StringSequence implements Sequence {

    // Sequence utilities
    public static char[] nucleotides = { 'A', 'C', 'G', 'U' };
    public static int alphabetLength = 4;

    // GU base pairs
    public final static String guBps = "GUG";

    // Canonical base pairs
    public final static String canonicalBps = "AUA CGC";

    // Canonical + GU base pairs
    public final static String canonicalGUBps = "AUA GCGUG";

    public static Character complement( Character c ) {
	switch( c ) {
	case 'A': return 'U';
	case 'C': return 'G';
	case 'U': return 'A';
	case 'G': return 'C';
	default: return ' ';
	}
    }

    public static Character wobble( Character c ) {
	switch( c ) {
	case 'U': return 'G';
	case 'G': return 'U';
	default: return ' ';
	}
    }

    public static char[]  allowedNts( char c, boolean WCComplement ) {
	char comp = complement( c );
	if( WCComplement ) return new char[]{ comp };

	char[] result = new char[3];
	int idx = 0;
	for( char nt: nucleotides )
	    if( nt != comp ) result[idx++] = nt;
	return result;
    }

    public static int countBPs( String sense, String antisense ) { // assume sense and antisense same size
	String antiantisense = Utils.reverse( antisense ); // reverse the antisense, as it is given 5'->3'
	int numBPs = 0;
	for( int i = 0; i < sense.length(); i++ )
	    if( canonicalGUBps.contains( Character.toString( sense.charAt( i ) ) + Character.toString( antiantisense.charAt( i ) ) ) ) numBPs++;
	return numBPs;
    }

    private String sequence = "";

    public StringSequence( String sequence ) { this.sequence = sequence; }
    
    // getters from Sequence interface
    public String getSequence() { return this.sequence; }
    
    public String getSequence( int start, int end ) {
	// substring [start..end[
	if( !( inSequence( start, this ) && inSequence( end - 1, this ) ) ) notInSequence( start, end, this );
	return this.sequence.substring( start, end );
    }
    
    public char getNucleotide( int position ) {
	if( !inSequence( position, this ) ) notInSequence( position, this );
	return this.sequence.charAt( position );
    }
    
    public int length() { return this.sequence.length(); }
    public void setSequence( String sequence ) { this.sequence = sequence; }

    public String toString() { return this.sequence; }

    // Utilities
    
    public static boolean inSequence( int position, Sequence sequence ) {
	return position >= 0 && position < sequence.length();
    }

    public static boolean inSequence( int start, int end, Sequence sequence ) {
	return start <= end && start >= 0 && end < sequence.length();
    }

    public static void notInSequence( int position, Sequence sequence ) throws IllegalArgumentException {
	throw new IllegalArgumentException( position + " is not in sequence range [0.." + (sequence.length() - 1) + "]" );
    }

    public static void notInSequence( int start, int end, Sequence sequence ) throws IllegalArgumentException {
	throw new IllegalArgumentException( "[" + start + ".." + end + "[ not in sequence range [0.." + (sequence.length() - 1) + "]" );
    }

    public static String toFasta( Sequence sequence ) {
	String fastaSequence = "";
	int i = 0;
	for( i = 0; i < sequence.length(); i = i + 60 )
	    fastaSequence += sequence.getSequence( i, Math.min( i+60, sequence.length() ) ) + "\n";
	return fastaSequence.substring( 0, fastaSequence.length() - 1 );  // remove the last '\n'
    }

    public static boolean isRNA( String sequence ) {
        // Iterate over each character in the sequence
        for( int i = 0; i < sequence.length(); i++ ) {
            char c = sequence.charAt( i );
            if( c != 'A' && c != 'C' && c != 'G' && c != 'U' ) {
                return false; // Found an invalid character
            }
        }
        return true; // No invalid characters found
    }

    public static String reverseComplement( String strand ) {
	String antiStrand = "";
	for( int k = 0; k < strand.length(); k++ )
	    antiStrand += complement( strand.charAt( k ) );
	return Utils.reverse( antiStrand );
    }

    public static List<String> reverseComplementWithAtMostOneWobble( String strand ) {
	// partials0: anti-strands built so far with 0 wobble used
	// partials1: anti-strands built so far with 1 wobble used
	List<StringBuilder> partials0 = new ArrayList<>();
	List<StringBuilder> partials1 = new ArrayList<>();
	partials0.add( new StringBuilder() );

	for( int i = 0; i < strand.length(); i++ ) {
	    char base  = strand.charAt( i );
	    char canon = complement( base );
	    char wob   = wobble( base );
	    boolean hasWobbleChoice = ( wob != ' ' ) && ( wob != canon );

	    List<StringBuilder> next0 = new ArrayList<>( partials0.size() ); // still 0 wobble used
	    List<StringBuilder> next1 = new ArrayList<>( partials1.size() + ( hasWobbleChoice ? partials0.size() : 0 ) ); // 1 wobble used

	    // Extend sequences that have used 0 wobble so far
	    for( StringBuilder sb : partials0 ) {
		// canonical keeps wobble count at 0
		next0.add( new StringBuilder( sb ).append( canon ) );
		// if wobble available, taking it consumes our single wobble
		if( hasWobbleChoice ) {
		    next1.add( new StringBuilder( sb ).append( wob ) );
		}
	    }

	    // Extend sequences that already used 1 wobble: only canonical allowed
	    for( StringBuilder sb : partials1 ) {
		next1.add( new StringBuilder( sb ).append( canon ) );
	    }

	    partials0 = next0;
	    partials1 = next1;
	}

	// Collect reverse complements (reverse the anti-strands)
	List<String> out = new ArrayList<>( partials0.size() + partials1.size() );
	for( StringBuilder sb : partials0 ) out.add( sb.reverse().toString() ); // 0 wobble
	for( StringBuilder sb : partials1 ) out.add( sb.reverse().toString() ); // 1 wobble
	return out;
    }

    // 3 or 4 in the seed, 3'-5' => 2 or 3
    public static List<String> reverseComplementWithOneDefaultAt2or3( String strand ) {
	// partials0: anti-strands built so far with 0 variant used
	// partials1: anti-strands built so far with 1 variant (wobble or mismatch) already used
	List<StringBuilder> partials0 = new ArrayList<>();
	List<StringBuilder> partials1 = new ArrayList<>();
	partials0.add(new StringBuilder());

	// Choose alphabet for mismatches
	final char[] MISMATCH = new char[]{'A','C','G','U'};

        for( int i = 0; i < strand.length(); i++ ) {
	    char base  = strand.charAt( i );
	    char canon = complement( base );  // your existing function
	    char wob   = wobble( base );      // your existing function (return ' ' when none)

	    boolean wobbleAvailable = ( wob != ' ' ) && ( wob != canon );
	    boolean variantAllowedHere = ( i == 2 || i == 3 );  // only positions 2 or 3 (0-based)

	    List<StringBuilder> next0 = new ArrayList<>( partials0.size() ); // still 0 variants used
	    List<StringBuilder> next1 = new ArrayList<>( partials1.size() + ( variantAllowedHere ? partials0.size()*3 : 0) ); // rough cap

	    // Extend sequences that have used 0 variants so far
	    for( StringBuilder sb : partials0 ) {
		// Always can add canonical without consuming the variant
		next0.add( new StringBuilder( sb ).append( canon ) );

		if( variantAllowedHere ) {
		    // (A) Use wobble here, if available (consumes the single variant)
		    if( wobbleAvailable ) {
			next1.add( new StringBuilder( sb ).append( wob ) );
		    }
		    // (B) Use a mismatch here (any base != canon). Also avoid duplicating wobble.
		    for( char b : MISMATCH ) {
			if( b != canon && ( !wobbleAvailable || b != wob ) ) {
			    next1.add( new StringBuilder( sb ).append( b ) );
			}
		    }
		}
	    }
	    // Extend sequences that already used their single variant: only canonical allowed now
	    for( StringBuilder sb : partials1 ) {
		next1.add( new StringBuilder( sb ).append( canon ) );
	    }

	    partials0 = next0;
	    partials1 = next1;
	}

	// Collect reverse complements (reverse the anti-strands)
	List<String> out = new ArrayList<>( partials1.size() );
	for( StringBuilder sb : partials1 ) out.add( sb.reverse().toString() ); // 1 variant (at index 3 or 4)
	//Utils.debug( out.toString() );
	return out;
    }

    public static List<String> reverseComplementWithDeletionAt2or3( String complement ) {
	List<String> out = new ArrayList<>(2);
	final int n = complement.length();

	// We allow a single deletion at 0-based positions 2 or 3 (if in range).
	List<Integer> delPositions = new ArrayList<>( List.of( 2, 3 ) );
	if( complement.charAt( 2 ) == complement.charAt( 3 ) ) delPositions.remove( 1 );
	
	for( int del : delPositions ) {
	    // Build the anti-strand forward (5'->3') by complementing every base
	    // except the deleted position, then reverse at the end to get the reverse-complement.
	    StringBuilder anti = new StringBuilder(n - 1);
	    for (int i = 0; i < n; i++) {
		if (i == del) continue;                 // skip (delete) this position
		char canon = complement(complement.charAt(i)); // your existing complement(base)
		anti.append(canon);
	    }
	    out.add(anti.reverse().toString());
	}
	return out;
    }

    public static List<String> reverseComplementWithInsertionAt3( String complement ) {
	List<String> out = new ArrayList<>();
	String prefix = "";
	for( int i = 0; i < 3; i++ ) prefix += complement( complement.charAt( i ) );
	String postfix = "";
	for( int i = 3; i < complement.length(); i++ ) postfix += complement( complement.charAt( i ) );
	
	String middle = prefix.toString() + "%s" + postfix.toString();
	for( char ins : new char[]{'A', 'C', 'G', 'U'} ) {
	    String s = String.format( middle, ins );
	    out.add( new StringBuilder( s ).reverse().toString() );
	}
	return out;
    }

    public static List<String> reverseComplementWithWobbles( String strand ) {
	// Build all anti-strands (not reversed) allowing wobble choices, then reverse at the end
	List<StringBuilder> partials = new ArrayList<>();
	partials.add( new StringBuilder() ); // start with empty

	for( int i = 0; i < strand.length(); i++ ) {
	    char base = strand.charAt( i );
	    char canon = complement( base );
	    char wob   = wobble( base );

	    boolean hasWobbleChoice = ( wob != ' ' ) && ( wob != canon );
	    List<StringBuilder> next = new ArrayList<>( hasWobbleChoice ? partials.size() * 2 : partials.size() );

	    for( StringBuilder sb : partials ) {
		// canonical option
		next.add( new StringBuilder( sb ).append( canon ) );
		// wobble option (if distinct)
		if( hasWobbleChoice ) {
		    next.add( new StringBuilder( sb ).append( wob ) );
		}
	    }
	    partials = next;
	}
	// Reverse each anti-strand to get reverse complements
	List<String> out = new ArrayList<>( partials.size() );
	for( StringBuilder sb : partials ) {
	    out.add( sb.reverse().toString() );
	}
	return out;
    }

    // generate all complementary sequences to argument 'complement'
    public static List<Integer> complementSites( String sequence, String complement, boolean wobble ) {
	List<Integer> result = new ArrayList<>();
	List<String> revComp;
	if( wobble ) revComp = reverseComplementWithWobbles( complement );
	else revComp = Arrays.asList( StringSequence.reverseComplement( complement ) );
	//System.out.println( "complement: " + complement + ", revComp: " + revComp );
	for( String rC: revComp ) {
	    int pos = 0;
	    int next = sequence.indexOf( rC, pos );
	    while( next != -1 ) {
		//System.out.println( "pos: " + pos + ", revComp: " + revComp );
		result.add( next );
		pos = next + 1;
		next = sequence.indexOf( rC, pos );
	    }
	}
	Collections.sort( result );
	return result;
    }

    // generate all complementary sequences to argument 'complement' with one wobble or mismatch at position 3 or 4
    public static List<Integer> complementSites( String sequence, String complement ) {
	List<Integer> result = new ArrayList<>();
	List<String> revComp = StringSequence.reverseComplementWithOneDefaultAt2or3( complement );
	//Utils.debug( "special complement: " + complement + ", revComp: " + revComp );
	for( String rC: revComp ) {
	    int pos = 0;
	    int next = sequence.indexOf( rC, pos );
	    while( next != -1 ) {
		result.add( next );
		pos = next + 1;
		next = sequence.indexOf( rC, pos );
	    }
	}
	Collections.sort( result );
	return result;
    }

    // generate all complementary sequences to argument 'complement' with one deletion at position 2 or 3
    public static List<Integer> complementWithDeletionAt2or3( String sequence, String complement ) {
	List<Integer> result = new ArrayList<>();
	List<String> revComp = StringSequence.reverseComplementWithDeletionAt2or3( complement );
	//Utils.debug( "delete complement: " + complement + ", revComp: " + revComp );
	for( String rC: revComp ) {
	    int pos = 0;
	    int next = sequence.indexOf( rC, pos );
	    while( next != -1 ) {
		result.add( next );
		pos = next + 1;
		next = sequence.indexOf( rC, pos );
	    }
	}
	Collections.sort( result );
	return result;
    }

    // generate all complementary sequences to argument 'complement' with one insertion between 3 and 4
    public static List<Integer> complementWithInsertionAt3( String sequence, String complement ) {
	List<Integer> result = new ArrayList<>();
	List<String> revComp = StringSequence.reverseComplementWithInsertionAt3( complement );
	//Utils.debug( "insert complement: " + complement + ", revComp: " + revComp );
	for( String rC: revComp ) {
	    int pos = 0;
	    int next = sequence.indexOf( rC, pos );
	    while( next != -1 ) {
		result.add( next );
		pos = next + 1;
		next = sequence.indexOf( rC, pos );
	    }
	}
	Collections.sort( result );
	return result;
    }

    // GC percent

    public static boolean isGCOkay( String sirna, double gcMin, double gcMax ) {
        // Count G and C nucleotides
        int gcCount = 0;
        for( char nucleotide : sirna.toCharArray() ) {
            if( nucleotide == 'G' || nucleotide == 'C' ) {
                gcCount++;
            }
        }
        // Calculate the percentage of GC content
	double gcPercentage = gcCount / (double) sirna.length();
	// Return true if GC content is between gcPercentMin and gcPercentMax, inclusive
        return  gcPercentage >= gcMin && gcPercentage <= gcMax;
    }

    // KMer utilities

    public static void generateKMers( int k, List<String> exclusions, List<String> kmers ) {
	generateKMersRec( nucleotides, "", alphabetLength, k, exclusions, kmers );
    }

    public static void generateKMers( int k, List<String> kmers ) {
	generateKMersRec( nucleotides, "", alphabetLength, k, new ArrayList<String>(), kmers );
    }

    private static boolean contains( String s, List<String> substrings ) {
	// check if contains one of the substrings
	for( String sub : substrings )
	    if( s.contains( sub ) ) return true;
	return false;
    }

    private static void generateKMersRec( char[] set, String prefix, int n, int k, List<String> exclusions, List<String> kmers ) {
	// Base case: k is  => add prefix to list
	if( k == 0 ) {
	    if( !contains( prefix, exclusions ) ) // if the prefix does not contain one of the exclusions
		kmers.add( prefix ); // add it to the list
	    return;
	}
	// Add each char from set and
	//     recursively call for k-1
	for( int i = 0; i < n; ++i ) {
	    // Add current char
	    String newPrefix = prefix + set[i];
	    generateKMersRec( set, newPrefix, n, k-1, exclusions, kmers ); // Decrease k and call recursively
	}
    }

}
