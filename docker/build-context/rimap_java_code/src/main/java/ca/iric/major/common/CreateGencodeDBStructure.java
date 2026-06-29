/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Arrays;
import java.util.Set;
import java.util.HashSet;
import java.util.Collections;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.StandardOpenOption;
import java.nio.charset.StandardCharsets;

// reads pct transcripts in a gencode version, ex gencode.v46.pc_transcripts.fa and create a DB structure
//    Use the two first characters of a transcript as its sub directory, ex) gencode.v46/SI
//    Use the transcript variant id to store the sequence, ex) gencode.v46/SI/SIRT1-201
public class CreateGencodeDBStructure {

    public static void main( String[] args ) {

	// ***** PRESENTATION
	System.out.println( "CreateGencodeDBStructure v.1.0 August 16, 2024 - Major Lab, IRIC, Université de Montréal" );

	// ***** READ ARGUMENTS

	if( args.length != 1 ) {
	    System.out.println( "usage: CreateGencodeDBStructure <genecode protein coding transcripts fasta file>" );
	    System.out.println( "   example of input file: gencode.v46.pc_transcripts.fa" );
	    Utils.stop( "bye!" , 0 );
	}

	// gencodePath
	String workingPath = "./data/";
	String gencodeFile = args[0]; // no validation! file must contain the protein coding transcript entries
	// args[0] example: gencode.v46.pc_transcripts.fa
	String gencodePath = workingPath + gencodeFile;
	
	// get the top directory name from gencodeFile
	int firstDot = gencodeFile.indexOf( "." );
	int secondDot = gencodeFile.indexOf( ".", firstDot + 1 );
	String topDir = gencodeFile.substring( 0, secondDot );
	String topDirPath = workingPath + topDir;

	System.out.println( "processing " + gencodePath + " ..." );

	// create top directory if it does not exist
	try {
	    Files.createDirectories( Paths.get( topDirPath ) );
	    System.out.println( topDirPath + " created successfully" );
	} catch( IOException e ) {
	    Utils.stop( "Failed to create " + topDirPath, 0 );
	}

	// read the PCT from gencode file, create a file for each variant, and dump its sequence in it
	try {
	    GencodePCTranscript genPCT = new GencodePCTranscript( gencodePath ); // read all PCTs in gencode file
	    for( ProteinCodingTranscript pct: genPCT ) {
		// create the subdirectory for the transcript, if it does not exist
		String subDirName = topDirPath + "/" + pct.getName().substring( 0, 2 );
		Files.createDirectories( Paths.get( subDirName ) );
		
		// open a file for the transcript variant and dump its fasta sequence in it
		String variantFile = subDirName + "/" + pct.getName();
		Path variantPath = Paths.get( variantFile );
		// if the file already exists, it will be overwritten
		Files.write( variantPath, pct.toString().getBytes( StandardCharsets.UTF_8 ),
			     StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING );
	    }
	    // create the top directory at path
	    System.out.println( "created the DB structure for " + genPCT.size() + " protein coding transcripts from " + gencodePath );
	} catch( FileNotFoundException e ) {
	    e.printStackTrace();
	} catch( IOException e ) {
	    e.printStackTrace();
	}
    }
}
