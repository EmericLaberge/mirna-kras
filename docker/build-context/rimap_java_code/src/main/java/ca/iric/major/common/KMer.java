/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Interface for KMers (composed of a Transcript and Sequence, extracted from the corresponding Transcript)
 *
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2021 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.lang.IllegalArgumentException;

public class KMer {

    // composed of Transcript and Sequence
    
    private Transcript transcript;  // KMer source
    private Sequence sequence;      // Result of this.transcript.getSequence( start, end )
    private int start, end;         // KMer range in source [start..end[
    private int occurrences;        // number of KMer occurrences

    // getters
    public Transcript getTranscript() { return this.transcript; }
    public Sequence getSequence()     { return this.sequence; }
    public int getOccurrences()       { return this.occurrences; }
    public int getStart()             { return this.start; }
    public int getEnd()               { return this.end; }
	
    // setters
    public void setOccurrences( int n ) { this.occurrences = n; }
    
    public void setStart( int start ) {
	if( !StringSequence.inSequence( start, this.sequence ) ) StringSequence.notInSequence( start, this.sequence );
	this.start = start;
    }

    public void setEnd( int end ) {
	if( !StringSequence.inSequence( end, this.sequence ) ) StringSequence.notInSequence( end, this.sequence );
	this.end = end;
    }

    public void setTranscript( Transcript transcript ) {
	this.transcript = transcript;
	// set start and end of the KMer in transcript
    }
    
    public void setSequence( Sequence sequence ) {
	if( this.transcript == null ) this.sequence = sequence;
	// make sure the KMer is within transcript and modify start and end accordingly
    }

    @Override
    public String toString() {
	// Returns the Ensembl/GeneCode corresponding header
	// >10
	return ">" +
	    this.getOccurrences() + "\n" +
	    StringSequence.toFasta( this.sequence );
    }
    
    // Constructors
    public KMer() {
	this.sequence = null;
	this.occurrences = -1;
	this.transcript = null;
    }
}
