/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.util.Scanner;
import java.util.Iterator;
import java.util.List;

public class MiRBaseMature implements Iterator<MatureMiRNA> {

    String fileName;
    FileReader reader;
    BufferedReader bReader;
    String header;

    public static void setHeader( String header, MatureMiRNA miRNA ) {
	// *** WARNING ***
	//     The sequence of the miRNA must be assigned before calling setHeader!
	//
	// header example:
	//      * 
	//    remove '>' and split at ' ' >hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p makes 6 data pieces
	//
	// hsa-let-7a-5p
	// MIMAT0000062
	// Homo
	// sapiens
	// let-7a-5p 
	String[] data = header.substring( 1, header.length() ).split( " " ); // remove '>' at the beginning of the header
	// data[0] = full name
	// data[1] = MIMAT code
	// data[2] = Homo
	// data[3] = sapiens
	// data[4] = name
	miRNA.setFullName( data[0] );
	miRNA.setMIMATCode( data[1] );
	miRNA.setOrganism( data[2] + " " + data[3] );
	miRNA.setName( data[4] );
    }

    public MiRBaseMature( String fileName ) throws IOException, FileNotFoundException {
	this.fileName = fileName;
	this.reader = new FileReader( fileName );
	this.bReader = new BufferedReader( this.reader );
	this.header = this.bReader.readLine();
    }

    public MiRBaseMature( String fileName, List<MatureMiRNA> list ) throws IOException, FileNotFoundException {
	this.fileName = fileName;
	this.reader = new FileReader( fileName );
	this.bReader = new BufferedReader( this.reader );
	this.header = this.bReader.readLine();
	while( this.header != null ) {
	    
	    // Read the sequence
	    String thisSequenceHeader = this.header;
	    String sequence = "";
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
	    MatureMiRNA miRNA = new MatureMiRNA( (Sequence) new StringSequence( sequence.replace( 'T', 'U' ) ) );
	    setHeader( thisSequenceHeader, miRNA );
	    list.add( miRNA );
	}
    }

    public MatureMiRNA next() {
	// header in this.header
	// Process header...
	
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
	MatureMiRNA miRNA = new MatureMiRNA( (Sequence) new StringSequence( sequence.replace( 'T', 'U' ) ) );
	setHeader( thisSequenceHeader, miRNA );
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
    
    
