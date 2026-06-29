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

public class MakeMatureMiRNAs {

    private static String path = "/Users/major/Dropbox/Dev/MiRBooking/data/Homo_sapiens_mature";

    public static void writeMature( MatureMiRNA miRNA ) {
	String matureName = miRNA.getName();

	// directory exists: let or mir
	// name includes let-, or mir-, first data piece of a split using '-'

	String[] data = matureName.split( "-" );
	String dirName = path + "/" + data[0];
	
	String fileName = dirName + "/" + matureName;
	File file = new File( fileName );

	try( BufferedWriter writer = new BufferedWriter( new FileWriter( fileName ) ) ) {
            writer.write( miRNA.toString() );
            //System.out.println( "Content written to file: " + fileName );
        } catch( IOException e ) {
            Utils.stop( "An error occurred while writing to file: " + fileName, 0 );
            e.printStackTrace();
        }
    }

    public static void main( String[] args ) {
	List<MatureMiRNA> myList = new ArrayList<>();
	try {
	    MiRBaseMature reader = new MiRBaseMature( path + ".fa", myList );
	    System.out.println( myList.size() + " mature MiRNA entries" );
	} catch( FileNotFoundException e ) {
	    Utils.stop( "microRNA " + path + ".fa" + " not found", 0 );
	} catch( IOException e ) {
	    Utils.stop( "problem reading " + path + ".fa", 0 );
	}
	for( MatureMiRNA miRNA : myList )
	    writeMature( miRNA );
    }
}
