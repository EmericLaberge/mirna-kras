/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * class for protein coding transcripts (composed of a Sequence)
 *
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2021 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.util.List;
import java.util.ArrayList;

import java.lang.IllegalArgumentException;

public class ProteinCodingTranscript extends CodingTranscript implements Transcript {

    // Composed of an Sequence

    private String geneId;          // ENSGXXXXXXXXXXX.X
    private String geneHId;         // OTTHUMGXXXXXXXXXXX.X
    private String transcriptHId;   // OTTHUMTXXXXXXXXXXX.X
    //private String name;            // HGNC-based nomenclature
    private String family;          // HGNC-based nomenclature

    // utilities

    public String getAntisense( int tlast, int t1, int bridge, int region, int length ) { // perfect antisense [tlast..tlast+7,tlast+bridge+1..t1]
	if( t1 < 0 || tlast >= this.sequence.length() ) return null;
	if( validateRegion( region, t1, tlast ) ) { // make sure the siRNA fits in region
	    int beginABox = tlast - 10; // (10 = 7 + 3; seed + box A lengths
	    int endSuppRg = beginABox - bridge + 2; // bridge includes box A and is therefore of minimum 3
	    if( ( tlast - beginABox + 1) + ( endSuppRg - t1 + 1) != length ) Utils.stop( "siRNA length invalid t1: " + t1 + " tlast: " + tlast + " bridge: " + bridge + " beginABox: " + beginABox + " endSuppRg: " + endSuppRg, 0 );
	    return
		StringSequence.reverseComplement( this.getSequence().getSequence( beginABox, tlast + 1 ) ) + // seed + A-Box
		StringSequence.reverseComplement( this.getSequence().getSequence( t1, endSuppRg + 1 ) ); // Supp (BC-Boxes) + D-Box
	}
	return null;
    }

    public List<Integer> seedMatch( String seed ) {
	List<Integer> result = new ArrayList<>();
        if( seed.isEmpty() ) return result;

	String transcript = this.getStringSequence();
        int index = 0;
        while( ( index = transcript.indexOf( seed, index ) ) != -1 )
	    result.add( index++ );
        return result;
    }

    public boolean isFound( String kmer ) {
	String transcript = this.getTargetableSequence();
	return transcript.indexOf( kmer ) != -1;
    }

    // getters
    public String getGeneId()         { return this.geneId; }
    public String getTranscriptId()   { return this.getId(); }
    public String getGeneHId()        { return this.geneHId; }
    public String getTranscriptHId()  { return this.transcriptHId; }
    public String getFamily()         { return this.family; }

    // setters
    public void setGeneId( String id )         { this.geneId = id; }
    public void setGeneHId( String id )        { this.geneHId = id; }
    public void setTranscriptHId( String id )  { this.transcriptHId = id; }
    public void setFamily( String family )     { this.family = family; }
    
    @Override
    public String toString() { return this.getName(); }
    public String toString2() {
	// Returns the Ensembl/GeneCode corresponding header and sequence
	// >ENST00000641515.2|ENSG00000186092.7|OTTHUMG00000001094.4|OTTHUMT00000003223.4|OR4F5-201|OR4F5|2618|UTR5:1-60|CDS:61-1041|UTR3:1042-2618|
	String utr5 = this.UTR5End > -1 ? "UTR5:" + (this.getUTR5Start() + 1) + "-" + (this.getUTR5End() + 1) + "|" : "";
	String utr3 = this.UTR3Start > -1 ? "UTR3:" + (this.getUTR3Start() + 1) + "-" + (this.getUTR3End() + 1) + "|\n" : "\n";
	return ">" +
	    this.getGeneId() + "|" +
	    this.getTranscriptId() + "|" +
	    this.getGeneHId() + "|" +
	    this.getTranscriptHId() + "|" +
	    this.getName() + "|" +
	    this.getFamily() + "|" +
	    this.sequence.length() + "|" +
	    utr5 +
	    "CDS:" + (this.getCDSStart() + 1) + "-" + (this.getCDSEnd() + 1) + "|" +
	    utr3 +
	    StringSequence.toFasta( this.sequence );
    }

    // Transcript interface methods
    
    // Constructor
    public ProteinCodingTranscript() {
	super();
    }
}
