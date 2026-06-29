/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.List;
import java.util.ArrayList;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileNotFoundException;

public class MakeTwoLetterGencodePCTranscripts {

    private static String path = "/Users/major/Dropbox/Dev/MiRBooking/data/Homo_sapiens.GRCh38.Gencode33.cdna";

    public static void writeGene( ProteinCodingTranscript pct ) {
	String geneName = pct.getName();
        if ( geneName == null || geneName.length() < 2 ) {
            Utils.stop( "Invalid gene name " + geneName, 0 );
            return;
        }

	// create a directory for it using the two first letters of its name, if it does not exist

        String dirName = path + "/" + geneName.substring( 0, 2 );
        File directory = new File( dirName );

        if( directory.exists() && directory.isDirectory() ) {
            System.out.println( "Directory already exists: " + dirName );
        } else {
            boolean isCreated = directory.mkdir();
            if( isCreated ) {
                System.out.println("Directory created: " + dirName);
            } else {
                System.out.println("Failed to create directory: " + dirName);
            }
        }

	// directory exists, write the pct to filename corresponding to its name
	
	String fileName = dirName + "/" + geneName;
	File file = new File( fileName );

	try( BufferedWriter writer = new BufferedWriter( new FileWriter( fileName ) ) ) {
            writer.write( pct.toString() );
            //System.out.println( "Content written to file: " + fileName );
        } catch( IOException e ) {
            Utils.stop( "An error occurred while writing to file: " + fileName, 0 );
            e.printStackTrace();
        }
    }

    public static void main( String[] args ) throws IOException, FileNotFoundException {
	List<ProteinCodingTranscript> myList = new ArrayList<>();
	GencodePCTranscript reader = new GencodePCTranscript( path + ".fa" );
	System.out.println( reader.size() + " Gencode entries" );
	for( ProteinCodingTranscript pct : reader )
	    writeGene( pct );
    }
}
