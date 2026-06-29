/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Interface for reporter transcripts (composed of a Sequence)
 *
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2024 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.lang.IllegalArgumentException;
import java.util.List;
import java.util.ArrayList;

public abstract class CodingTranscript implements Transcript {

    // Composed of an Sequence

    protected Sequence sequence;      // Sequence holder
    protected Biotype biotype;
    protected String name;            // HGNC-based nomenclature
    protected String id;              // ENST-based id
    protected String length;          // Read from Gencode file (a String)
    protected int UTR5Start = -1, UTR5End = -1;
    protected int CDSStart = -1, CDSEnd = -1;
    protected int UTR3Start = -1, UTR3End = -1;

    // utilities
    public boolean inUTR5( int position ) { return position >= 0 && position <= this.UTR5End; }
    public boolean inCDS( int position ) { return position >= this.CDSStart && position <= this.CDSEnd; }
    public boolean inUTR3( int position ) { return position >= this.UTR3Start && position <= this.UTR3End; }

    public static String printRegions( int region ) {
	switch( region ) {
	case 1: return "3UTR";
	case 2: return "CDS";
	case 3: return "CDS and 3UTR";
	case 4: return "5UTR";
	case 5: return "5UTR and 3UTR";
	case 6: return "5UTR and CDS";
	case 7: return "5UTR, 3UTR and CDS";
	default: Utils.stop( "invalid input region " + region, 0 );
	}
	return "";
    }

    public int regionsLength( int region ) {
	switch( region ) {
	case 1: return UTR3End - UTR3Start + 1; // "3'UTR";
	case 2: return CDSEnd - CDSStart + 1; // "CDS";
	case 3: return UTR3End == -1 ? CDSEnd - CDSStart + 1 : UTR3End - CDSStart + 1; // "CDS and 3'UTR";
	case 4: return UTR5End - UTR5Start + 1; // "5'UTR";
	case 5: return (UTR5End - UTR5Start + 1 ) + (UTR3End - UTR3Start + 1 ); //  "5' and 3'UTRs";
	case 6: return CDSEnd - UTR5Start + 1; // "5'UTR and CDS";
	case 7: return UTR3End - UTR5Start + 1; // "5' and 3'UTRs and CDS";
	default: Utils.stop( "invalid input region " + region, 0 );
	}
	return -1;
    }

    public String regionsSequence( int region ) {
	switch( region ) {
	case 1: return this.sequence.getSequence().substring( UTR3Start, UTR3End + 1 ); // "3'UTR";
	case 2: return this.sequence.getSequence().substring( CDSStart, CDSEnd + 1 ); // "CDS";
	case 3: return this.regionsSequence( 2 ) + this.regionsSequence( 1 ); // "CDS and 3'UTR";
	case 4: return this.sequence.getSequence().substring( UTR5Start, UTR5End + 1); // "5'UTR";
	case 5: return this.regionsSequence( 4 ) + this.regionsSequence( 1 );
	case 6: return this.regionsSequence( 4 ) + this.regionsSequence( 2 ); // "5'UTR and CDS";
	case 7: return this.sequence.getSequence(); // "5' and 3'UTRs and CDS";
	default: Utils.stop( "invalid input region " + region, 0 );
	}
	return "";
    }

    public boolean validateRegion( int region, int start, int end ) {
	if( region == 7 ) return true;
	switch( region ) {
	case 1: return inUTR3( start ) && inUTR3( end );
	case 2: return inCDS( start ) && inCDS( end );
	case 3: return validateRegion( 1, start, end ) || validateRegion( 2, start, end ) || ( inCDS( start ) && inUTR3( end ) );
	case 4: return inUTR5( start ) && inUTR5( end );
	case 5: return validateRegion( 4, start, end ) || validateRegion( 1, start, end );
	case 6: return validateRegion( 4, start, end ) || validateRegion( 2, start, end ) || ( inUTR5( start ) && inCDS( end ) );
	default: Utils.stop( "invalid input region " + region, 0 );
    }
	return false;
    }

    public int determineRegion( int start, int end ) {
	if( inUTR3( start ) && inUTR3( end ) ) return 1;
	if( inCDS( start ) && inCDS( end ) ) return 2;
	if( validateRegion( 1, start, end ) || validateRegion( 2, start, end ) || ( inCDS( start ) && inUTR3( end ) ) ) return 3;
	if( inUTR5( start ) && inUTR5( end ) ) return 4;
	if( validateRegion( 4, start, end ) || validateRegion( 1, start, end ) ) return 5;
	if( validateRegion( 4, start, end ) || validateRegion( 2, start, end ) || ( inUTR5( start ) && inCDS( end ) ) ) return 6;
	return 7;
    }

    private void notInSequence( int position ) throws IllegalArgumentException {
	throw new IllegalArgumentException( position + " is not in sequence[0.." + (this.sequence.length() - 1) + "]" );
    }

    private void notInSequence( int start, int end ) throws IllegalArgumentException {
	throw new IllegalArgumentException( "[" + start + ".." + end + "[ not in sequence[0.." + (this.sequence.length() - 1) + "]" );
    }

    public boolean inTranscript( int start, int end ) { return StringSequence.inSequence( start, end, this.sequence ); }

    // getters
    public String getName()               { return this.name; }
    public String getId()                 { return this.id; }
    public String getLength()             { return this.length; }
    public int    getIntLength()          { return Integer.parseInt( this.length ); }
    public int    getUTR5Start()          { return this.UTR5Start; }
    public int    getCDSStart()           { return this.CDSStart; }
    public int    getUTR3Start()          { return this.UTR3Start; }
    public int    getUTR5End()            { return this.UTR5End; }
    public int    getCDSEnd()             { return this.CDSEnd; }
    public int    getUTR3End()            { return this.UTR3End; }
    public String getUTR5()               { return this.sequence.getSequence( this.getUTR5Start(), this.getUTR5End()+1 ); }
    public String getCDS()                { return this.sequence.getSequence( this.getCDSStart(), this.getCDSEnd()+1 ); }
    public String getUTR3()               { return this.sequence.getSequence( this.getUTR3Start(), this.getUTR3End()+1 ); }
    public String getTargetableSequence() {
	// two cases: with and without 3'UTR
	if( this.getUTR3End() > -1 ) // there is a 3,UTR
	    return this.sequence.getSequence( this.getUTR5Start(), this.getUTR3End()+1 ); // CDS+3'UTR
	else // no 3,UTR
	    return this.sequence.getSequence( this.getUTR5Start(), this.getCDSEnd()+1 ); // CDS only
    }
    
    public int    getRegion( int start, int end ) {
	// Returns region: 1(UTR5), 2(CDS), 4(UTR3), 3(overlap UTR5/CDS), 6(overlap CDS/UTR3), 7(overlap all)
	if( !( StringSequence.inSequence( start, this.sequence ) && StringSequence.inSequence( end, this.sequence ) ) ) notInSequence( start, end );
	if( inUTR5( start ) ) {
	    if( inUTR5( end ) ) return 4;
	    if( inCDS( end ) ) return 6;
	    return 7;
	}
	if( inCDS( start ) ) {
	    if( inCDS( end ) ) return 2;
	    return 3;
	}
	return 1;
    }
    public String stringRegion( int start, int end ) {
	return this.printRegions( this.getRegion( start, end ) );
    }
	
    // setters
    public void setName( String name )         { this.name = name; }
    public void setId( String id )             { this.id = id; }
    public void setLength( String length ) {
	if( !Integer.toString( this.sequence.length() ).equals( length ) )
	    Utils.stop( "Sequence length " + this.sequence.length() + " does not match Gencode data " + length + "!", 1 );
	this.length = length;
    }
    
    public void setUTR5( int start, int end )  { this.setUTR5Start( start ); this.setUTR5End( end ); }
    public void setCDS( int start, int end )   { this.setCDSStart( start ); this.setCDSEnd( end ); }
    public void setUTR3( int start, int end )  { this.setUTR3Start( start ); this.setUTR3End( end ); }
	
    public void setUTR5Start( int start ) {
	if( !StringSequence.inSequence( start, this.sequence ) ) notInSequence( start );
	this.UTR5Start = start;
    }
    public void setCDSStart( int start ) {
	if( !StringSequence.inSequence( start, this.sequence ) ) notInSequence( start );
	this.CDSStart = start;
    }
    public void setUTR3Start( int start ) {
	if( !StringSequence.inSequence( start, this.sequence ) ) notInSequence( start );
	this.UTR3Start = start;
    }
    public void setUTR5End( int end ) {
	if( !StringSequence.inSequence( end, this.sequence ) ) notInSequence( end );
	this.UTR5End = end;
    }
    public void setCDSEnd( int end ) {
	if( !StringSequence.inSequence( end, this.sequence ) ) notInSequence( end );
	this.CDSEnd = end;
    }
    public void setUTR3End( int end ) {
	if( !StringSequence.inSequence( end, this.sequence ) ) notInSequence( end );
	this.UTR3End = end;
    }

    // Interface Transcript
    public Biotype getBiotype()                  { return this.biotype; }
    public Sequence getSequence()                { return this.sequence; }
    public String getStringSequence()            { return this.sequence.getSequence(); }
    public void setBiotype( Biotype type )       { this.biotype = type; }
    public void setSequence( Sequence sequence ) { this.sequence = sequence; }

    @Override
    public String toString() {
	// Returns the header and sequence
	// >NAME|355|UTR5:1-60|CDS:61-301|UTR3:302-355|
	String utr5 = this.UTR5End > -1 ? "UTR5:" + (this.getUTR5Start() + 1) + "-" + (this.getUTR5End() + 1) + "|" : "";
	String utr3 = this.UTR3Start > -1 ? "UTR3:" + (this.getUTR3Start() + 1) + "-" + (this.getUTR3End() + 1) + "|\n" : "\n";
	return ">" +
	    this.getName() + "|" +
	    this.sequence.length() + "|" +
	    utr5 +
	    "CDS:" + (this.getCDSStart() + 1) + "-" + (this.getCDSEnd() + 1) + "|" +
	    utr3 +
	    StringSequence.toFasta( this.sequence );
    }

    @Override
    public boolean equals( Object o ) {
	if( this == o ) return true;
	if( !(o instanceof CodingTranscript ) ) return false;
	CodingTranscript that = (CodingTranscript) o;
	return this.name.equals( that.name );
    }

    @Override
    public int hashCode() {
	return this.name.hashCode();
    }

    // Constructor
    public CodingTranscript() {
	this.biotype = Biotype.reporter;
	this.length = ""; // initialized to the empty string in case not read from Gencode file
	this.UTR5Start = this.UTR5End = this.UTR3Start = this.UTR3End = -1; // initialize to -1 in case no UTRs
    }
}
