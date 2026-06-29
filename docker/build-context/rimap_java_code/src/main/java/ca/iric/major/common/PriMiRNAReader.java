/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.util.Scanner;
import java.util.Iterator;
import java.util.List;

public class PriMiRNAReader implements Iterator<PriMiRNA> {

    String fileName;
    FileReader reader;
    BufferedReader bReader;
    String header;
    char splitter;

    public void setHeader( String header, PriMiRNA miRNA, char splitter ) {
	// *** WARNING ***
	//     The sequence of the miRNA must be assigned before calling setHeader!
	//
	// header example:
	//      * 
	//    split at ' ' > hsa-let-7a-5p makes 2 data pieces or >hsa-let-7a-5p
	//
	// hsa-let-7a-5p
	String[] data = header.split( ""+this.splitter ); // remove '>' at the beginning of the header
	// data[1] == name
	miRNA.setName( data[1] );
    }

    public PriMiRNAReader( String fileName, char splitter ) {
	this.fileName = fileName;
	this.splitter = splitter;
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

    public PriMiRNAReader( String fileName, List<PriMiRNA> list, char splitter ) {
	this.fileName = fileName;
	this.splitter = splitter;
	try {
	    this.reader = new FileReader( fileName );
	} catch( FileNotFoundException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
	try {
	    this.bReader = new BufferedReader( this.reader );
	    this.header = this.bReader.readLine();
	    while( this.header != null ) {
		PriMiRNA miRNA = new PriMiRNA();
	    
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
		miRNA.setSequence( (Sequence) new StringSequence( sequence.replace( 'T', 'U' ) ) );
		setHeader( thisSequenceHeader, miRNA, this.splitter );
		list.add( miRNA );
	    }
	} catch( IOException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
    }

    public PriMiRNA next() {
	// header in this.header
	// Process header...
	PriMiRNA miRNA = new PriMiRNA();
	
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
	miRNA.setSequence( (Sequence) new StringSequence( sequence.replace( 'T', 'U' ) ) );
	setHeader( thisSequenceHeader, miRNA, this.splitter );
	return miRNA;
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
    
    
