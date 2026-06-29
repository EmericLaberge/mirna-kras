/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.util.Scanner;
import java.util.Iterator;
import java.util.List;

public class RNASeqReads implements Iterator<Read> {

    String fileName;
    FileReader reader;
    BufferedReader bReader;
    String header;

    public static void setHeader( String header, Read read ) {
	// *** WARNING ***
	//     The sequence of the Read must be assigned before calling setHeader!
	//
	// header example:
	//      @SRR18046904.1 1 length=79
	//
	String[] data = header.substring( 1, header.length() ).split( " " ); // remove '@' at the beginning of the header
	// data = occurrences
	read.setId( data[0] );
    }

    public RNASeqReads( String fileName ) {
	this.fileName = fileName;
	try {
	    this.reader = new FileReader( fileName );
	} catch( FileNotFoundException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
	try {
	    this.bReader = new BufferedReader( this.reader );
	    this.header = this.bReader.readLine();
	} catch( IOException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
    }

    public RNASeqReads( String fileName, List<Read> list ) {
	this.fileName = fileName;
	long avgLength = 0;
	long numberOfSequences = 0;
	try {
	    this.reader = new FileReader( fileName );
	} catch( FileNotFoundException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
	try {
	    this.bReader = new BufferedReader( this.reader );
	    this.header = this.bReader.readLine();
	    while( this.header != null ) {
		Read read = new Read();
	    
		// Read the sequence
		String thisSequenceHeader = this.header;
		String sequence = "";
		try {
		    String line = bReader.readLine();
		    while( line != null && ! line.contains( "+" ) ) {
			sequence += line;
			line = bReader.readLine();
		    }
		    // line contains last sequence line or the next header if starts with >
		    if( line == null ) {
			this.header = null;
		    }
		    else
			if( line.contains( "+" ) ) { // announce qualtity values for the sequence, skip this one and the next one
			    line = bReader.readLine();
			    line = bReader.readLine();
			}
		    else
			if( line.contains( "@" ) ) {
			    this.header = line;
			}
		} catch( IOException e ) {
		    ReaderUtilities.processException( e, fileName );
		}
		// find the longest sequence free of N
		String[] splitSeq = sequence.split( "N" );
		String longest = "";
		for( String s : splitSeq )
		    if( s.length() > longest.length() )
			longest = s;
		read.setSequence( new StringSequence( longest.replace( 'T', 'U' ) ) );
		avgLength += longest.length();
		numberOfSequences++;
		setHeader( thisSequenceHeader, read);
		list.add( read );
	    }
	} catch( IOException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
	System.out.println( "average length = " + ( (double)avgLength / (double)numberOfSequences ) );
    }

    public Read next() {
	// header in this.header
	// Process header...
	Read read = new Read();
	
	// Read the sequence
	String thisSequenceHeader = this.header;
	String sequence = "";
	try {
	    String line = bReader.readLine();
	    while( line != null && ! line.contains( ">" ) ) {
		sequence += line;
		line = bReader.readLine();
	    }
	    // line contains last sequence line or the next header if starts with >
	    if( line == null ) {
		this.header = null;
	    }
	    else
		if( line.contains( ">" ) ) {
		    this.header = line;
		}
	} catch( IOException e ) {
	    ReaderUtilities.processException( e, this.fileName );
	}
	read.setSequence( new StringSequence( sequence.replace( 'T', 'U' ) ) );
	setHeader( thisSequenceHeader, read );
	return read;
    }

    public boolean hasNext() {
	if( this.header != null ) return true;
	try {
	    this.bReader.close();
	} catch( IOException e ) {
	    ReaderUtilities.processException( e, this.fileName );
	}
	return false;
    }
}
    
