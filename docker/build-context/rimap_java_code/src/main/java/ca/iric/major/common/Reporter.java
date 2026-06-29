/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.util.Scanner;
import java.util.Iterator;
import java.util.List;

public class Reporter implements Iterator<ReporterTranscript> {

    String fileName;
    FileReader reader;
    BufferedReader bReader;
    String header;

    public static void setHeader( String header, ReporterTranscript transcript ) {
	// *** WARNING ***
	//     The sequence of the transcript must be assigned before calling setHeader!
	//
	// header example:
	//      * 
	//    split at '|' makes 3 to 5 Strings, depends on regions (can be 3 or 1 if CDS only); assume order is UTR5, CDS, UTR3
	// >NAME|355|UTR5:1-60|CDS:61-301|UTR3:302-355|
	String[] data = header.substring( 1, header.length() ).split( "\\|" ); // remove '>' at the beginning of the header
	// data[0] = name
	// data[1] = length
	// data[2] = utr5 (must split at :)
	// data[2] = cds (must split at :); if UTR5 is absent
	// data[3] = cds (must split at :); if UTR5 is present
	// data[4] = utr3 (must split at :);if UTR5 is present
	transcript.setName( data[0] );
	transcript.setLength( data[1] );
	// There is at least a CDS; parse "REGION:start-end"; Gencode sequences start at 1
	String[] range;
	if( data[2].contains( "UTR5" ) ) {
	    range = data[2].split( ":", 2 )[1].split( "-", 2 );
	    transcript.setUTR5( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	}
	else if( data[2].contains( "CDS" ) ) {
	    range = data[2].split( ":", 2 )[1].split( "-", 2 );
	    transcript.setCDS( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	}
	if( data.length > 3 ) {
	    if( data[3].contains( "CDS" ) ) {
		range = data[3].split( ":", 2 )[1].split( "-", 2 );
		transcript.setCDS( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	    }
	    else if( data[3].contains( "UTR3" ) ) {
		range = data[3].split( ":", 2 )[1].split( "-", 2 );
		transcript.setUTR3( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	    }
	}
	if( data.length > 4 ) { // contains the UTR3
	    range = data[4].split( ":", 2 )[1].split( "-", 2 );
	    transcript.setUTR3( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	}
    }

    public Reporter( String fileName ) {
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

    public Reporter( String fileName, List<ReporterTranscript> list ) {
	this.fileName = fileName;
	try {
	    this.reader = new FileReader( fileName );
	} catch( FileNotFoundException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
	try {
	    this.bReader = new BufferedReader( this.reader );
	    this.header = this.bReader.readLine();
	    while( this.header != null ) {
		ReporterTranscript transcript = new ReporterTranscript();
	    
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
		    ReaderUtilities.processException( e, fileName );
		}
		transcript.setSequence( (Sequence) new StringSequence( sequence.replace( 'T', 'U' ) ) );
		setHeader( thisSequenceHeader, transcript );
		list.add( transcript );
	    }
	} catch( IOException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
    }

    public ReporterTranscript next() {
	// header in this.header
	// Process header...
	ReporterTranscript transcript = new ReporterTranscript();
	
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
	transcript.setSequence( (Sequence) new StringSequence( sequence.replace( 'T', 'U' ) ) );
	setHeader( thisSequenceHeader, transcript );
	return transcript;
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
    
    
