/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.util.Scanner;
import java.util.Iterator;
import java.util.List;

public class GencodeKMer implements Iterator<KMer> {

    String fileName;
    FileReader reader;
    BufferedReader bReader;
    String header;

    public static void setHeader( String header, KMer kmer ) {
	// *** WARNING ***
	//     The sequence of the kmer must be assigned before calling setHeader!
	//
	// header example:
	//      * 
	//    
	// >10
	String data = header.substring( 1, header.length() ); // remove '>' at the beginning of the header
	// data = occurrences
	kmer.setOccurrences( Integer.parseInt( data ) );
    }

    public GencodeKMer( String fileName ) {
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

    public GencodeKMer( String fileName, List<KMer> list ) {
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
		KMer kmer = new KMer();
	    
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
		kmer.setSequence( new StringSequence( sequence.replace( 'T', 'U' ) ) );
		setHeader( thisSequenceHeader, kmer );
		list.add( kmer );
	    }
	} catch( IOException e ) {
	    ReaderUtilities.processException( e, fileName );
	}
    }

    public KMer next() {
	// header in this.header
	// Process header...
	KMer kmer = new KMer();
	
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
	kmer.setSequence( new StringSequence( sequence.replace( 'T', 'U' ) ) );
	setHeader( thisSequenceHeader, kmer );
	return kmer;
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
    
