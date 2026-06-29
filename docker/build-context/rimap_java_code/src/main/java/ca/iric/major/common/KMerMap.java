/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * KMerMap stores kmer positions, common to a set of transcripts
 * The position of a kmer is the position of its first character in a transcript
 * The positions are restricted to the region (5' and/or 3'UTR, and/or CDS) received as an argument in the constructor
 *
 * @version 3.0
 * @author Francois Major
 * @copyright 1.0 2021 - MajorLab, IRIC, Universite de Montreal
 * @copyright 2.0 2024 - MajorLab, IRIC, Universite de Montreal
 * @copyright 3.0 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.lang.IllegalArgumentException;

import java.util.Set;
import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

public class KMerMap {

    // Method to encode a k-mer string into an integer
    public static int kmerToInt( String kmer ) {
        int index = 0;
        for( int i = 0; i < kmer.length(); i++ ) {
            char nucleotide = kmer.charAt(i);
            int value = 0;
            switch( nucleotide ) {
	    case 'A': value = 0; break;
	    case 'C': value = 1; break;
	    case 'G': value = 2; break;
	    case 'U': value = 3; break;
	    default: throw new IllegalArgumentException("Invalid nucleotide: " + nucleotide);
            }
            index = ( index << 2 ) | value; // Shift left by 2 bits (base 4) and add the value
        }
        return index;
    }

    // Method to decode an integer back into a k-mer string
    public static String intToKmer( int index, int k ) {
        StringBuilder kmer = new StringBuilder();
        for( int i = 0; i < k; i++ ) {
            int value = index & 3; // Get the last 2 bits (base 4)
            char nucleotide = 'A'; // Default to 'A'
            switch (value) {
	    case 0: nucleotide = 'A'; break;
	    case 1: nucleotide = 'C'; break;
	    case 2: nucleotide = 'G'; break;
	    case 3: nucleotide = 'U'; break;
            }
            kmer.insert(0, nucleotide); // Insert at the beginning
            index = index >> 2; // Shift right by 2 bits
        }
        return kmer.toString();
    }

    public static List<Integer> findKmerPositions( CodingTranscript t, String kmer, int region ) {
        List<Integer> positions = new ArrayList<>();
	String sequence = t.getStringSequence();
        int index = sequence.indexOf( kmer ); // find first instance
        while( index >= 0 ) {
	    if( t.validateRegion( region, index, index + kmer.length() -1 ) )
		positions.add( index );
            index = sequence.indexOf( kmer, index + 1 ); // continue searching from the next character
        }
        return positions;
    }

    public List<String> getExclusive( Set<ProteinCodingTranscript> inSet, Set<ProteinCodingTranscript> notInSet ) {
	List<String> kmerList = new ArrayList<>();
	for( int i = 0; i < Math.pow( 4, this.k); i++ ) { // for each kmer
	    String kmer = intToKmer( i, this.k ); // kmer String
	    if( this.isInAll( i, inSet ) && isNotIn( kmer, notInSet ) )
		kmerList.add( kmer );
	}
	return kmerList;
    }

    public List<String> getExclusive( Set<ProteinCodingTranscript> inSet ) {
	List<String> kmerList = new ArrayList<>();
	// search kmers
	for( int i = 0; i < Math.pow( 4, this.k); i++ ) { // for each kmer
	    String kmer = intToKmer( i, this.k ); // kmer String
	    if( this.isInAll( i, inSet ) )
		kmerList.add( kmer );
	}
	return kmerList;
    }

    protected List<ProteinCodingTranscript> transcripts;               // transcript set
    protected Map<ProteinCodingTranscript,Integer> transcriptsIndices; // transcript indices for the positions' lists
    protected int region;                                              // region of the transcript considered
    protected List<String> exclusions;                                 // excluded kmers
    private List<Integer>[][] positions;                               // positions of the kmers in all transcripts, first dimension kmer, second dimension transcript index
    private int k;                                                     // length of the kmer
    
    // Constructor
    public KMerMap( List<ProteinCodingTranscript> transcripts, int region, List<String> exclusions, int k ) {

	// initialize the attributes
	this.transcripts = transcripts;
	this.region = region;
	this.exclusions = exclusions;
	this.k = k;
	this.transcriptsIndices = new HashMap<>();

	// initialize the arrays of positions
	int depth = (int)Math.pow( 4, this.k );
	this.positions = new ArrayList[depth][this.transcripts.size()];
	for( int i = 0; i < depth; i++ )
	    for( int j = 0; j < this.transcripts.size(); j++ )
		positions[i][j] = new ArrayList<>();

	// initialize transcript indices
	for( int i = 0; i < this.transcripts.size(); i++ ) { // iterate the transcripts
	    this.transcriptsIndices.put( this.transcripts.get( i ), i );
	}

	// find positions
	for( int i = 0; i < depth; i++ ) { // for each kmer
	    String kmer = intToKmer( i, this.k ); // kmer String
	    for( int j = 0; j < this.transcripts.size(); j++ ) // iterate the transcripts
		if( !containsAnySubstring( kmer, exclusions ) )
		    this.positions[i][j] = findKmerPositions( this.transcripts.get( j ), kmer, region );
	}
    }

    public int getNumberOfTranscripts() { return this.transcripts.size(); }
    public int getK() { return this.k; }
    public List<ProteinCodingTranscript> getTranscripts() { return this.transcripts; }
    // get the positions of kmers in a given transcript
    public List<Integer> getPositions( String kmer, ProteinCodingTranscript pct ) {
	return this.positions[kmerToInt( kmer )][this.transcriptsIndices.get( pct )];
    }
    // get the positions of a kmer in a transcript by its index
    public List<Integer> getPositions( String kmer, int pct ) { return this.positions[kmerToInt( kmer )][pct]; }
    public int getRegion() { return this.region; }
    public List<String> getExclusions() { return this.exclusions; }    
    public ProteinCodingTranscript getTarget( int i ) { return this.transcripts.get( i ); } // get transcript at index i

    // Utilities

    private static boolean containsAnySubstring( String kmer, List<String> substrings ) {
        for( String substring : substrings ) {
            if( kmer.contains( substring ) ) {
                return true; // return true as soon as a match is found
            }
        }
        return false; // return false if no matches are found
    }

    private static boolean isNotIn( String kmer, Set<ProteinCodingTranscript> excluded ) {
	for( ProteinCodingTranscript t : excluded )
	    if( t.getTargetableSequence().contains( kmer ) ) return false;
	return true;
    }

    private boolean isInAll( int kmerIndex, Set<ProteinCodingTranscript> transcripts ) {
	for( ProteinCodingTranscript t : transcripts ) if( this.positions[kmerIndex][this.transcriptsIndices.get( t )].size() == 0 ) return false;
	return true;
    }

    public String toString() {
	String result = "";
	// visualize map
	for( int i = 0; i < Math.pow( 4, this.k ); i++ ) { // for each kmer
	    String kmer = intToKmer( i, this.k ); // kmer String
	    for( int j = 0; j < this.transcripts.size(); j++ ) // iterate the transcripts
		if( !this.positions[i][j].isEmpty() ) result += kmer + " in " + this.transcripts.get( j ).getName() + ": " + this.positions[i][j] + "\n";
	}
	return result;
    }
}
