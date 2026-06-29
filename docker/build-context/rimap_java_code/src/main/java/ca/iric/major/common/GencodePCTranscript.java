/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.util.Scanner;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;

public class GencodePCTranscript implements Iterable<ProteinCodingTranscript>{

    String fileName;
    FileReader reader;
    BufferedReader bReader;
    String header;
    List<ProteinCodingTranscript> transcripts;

    // getters
    public String getFileName() { return this.fileName; }
    public int size() { return this.transcripts.size(); }
    public ProteinCodingTranscript get( int i ) { return this.transcripts.get( i ); }
    public List<ProteinCodingTranscript> getList() { return this.transcripts; }
    public Iterator<ProteinCodingTranscript> iterator() { return this.transcripts.iterator(); }
    
    public static void setHeader( String header, ProteinCodingTranscript transcript ) {
	// *** WARNING ***
	//     The sequence of the transcript must be assigned before calling setHeader!
	//
	// header example:
	//      * 
	//    split at '|' makes 8 to 10 Strings, depends on regions (can be 3 or 1 if CDS only); assume order is UTR5, CDS, UTR3
	// >ENST00000212015.11|ENSG00000096717.12|OTTHUMG00000018340.3|OTTHUMT00000048296.2|SIRT1-201|SIRT1|4107|UTR5:1-66|CDS:67-2310|UTR3:2311-4107|
	String[] data = header.substring( 1, header.length() ).split( "\\|" ); // remove '>' at the beginning of the header
	// data[0] = transcriptId
	// data[1] = geneId
	// data[2] = geneId; Havana
	// data[3] = transcriptId; Havana
	// data[4] = transcript variant
	// data[5] = gene name
	// data[6] = length
	// data[7] = utr5 (must split at :)
	// data[7] = cds (must split at :); if UTR5 is absent
	// data[8] = cds (must split at :); if UTR5 is present
	// data[9] = utr3 (must split at :);if CDS is present
	transcript.setGeneId( data[0] );
	transcript.setId( data[1] );
	transcript.setGeneHId( data[2] );
	transcript.setTranscriptHId( data[3] );
	transcript.setName( data[4] );
	transcript.setFamily( data[5] );
	transcript.setLength( data[6] );
	// There is at least a CDS; parse "REGION:start-end"; Gencode sequences start at 1
	String[] range;
	if( data[7].contains( "UTR5" ) ) {
	    range = data[7].split( ":", 2 )[1].split( "-", 2 );
	    transcript.setUTR5( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	}
	else if( data[7].contains( "CDS" ) ) {
	    range = data[7].split( ":", 2 )[1].split( "-", 2 );
	    transcript.setCDS( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	}
	if( data.length > 8 ) {
	    if( data[8].contains( "CDS" ) ) {
		range = data[8].split( ":", 2 )[1].split( "-", 2 );
		transcript.setCDS( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	    }
	    else if( data[8].contains( "UTR3" ) ) {
		range = data[8].split( ":", 2 )[1].split( "-", 2 );
		transcript.setUTR3( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	    }
	}
	if( data.length > 9 ) { // contains the UTR3
	    range = data[9].split( ":", 2 )[1].split( "-", 2 );
	    transcript.setUTR3( Integer.parseInt( range[0] ) - 1, Integer.parseInt( range[1] ) - 1 );
	}
    }

    public GencodePCTranscript( String fileName ) throws IOException, FileNotFoundException {
	this.fileName = fileName;
	this.reader = new FileReader( fileName );
	this.bReader = new BufferedReader( this.reader );
	this.header = this.bReader.readLine();
	this.transcripts = new ArrayList<>();
	while( this.header != null ) {
	    ProteinCodingTranscript transcript = new ProteinCodingTranscript();
	    
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
	    transcript.setSequence( (Sequence) new StringSequence( sequence.replace( 'T', 'U' ) ) );
	    setHeader( thisSequenceHeader, transcript );
	    this.transcripts.add( transcript );
	}
    }
}
    
    
