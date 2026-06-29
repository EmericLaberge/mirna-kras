/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import ca.iric.major.common.*;

import java.util.*;
import java.util.stream.Collectors;

import java.io.*;
import java.nio.file.*;

import java.text.ParseException;

public class RIMapRescueAccessibilityFiles {

    public static void main( String[] args ) {

	// ***** READ ARGUMENTS

	if( args.length == 0 || args.length > 1 ) {
	    System.out.println( "usage: RIMapRescueAccessibilityFiles <Gencode version>" );
	    System.out.println( "   Gencode version: ex) v46" );
	    Utils.stop( "bye!" , 0 );
	}

	char fileSystemSeparator = '/';
	// ***** READ OS and use appropriate file system separator
	String os = System.getProperty( "os.name" ).toLowerCase();
        if( os.contains( "win" ) ) {
            System.out.println( "running under Windows" );
	    fileSystemSeparator = '\\';
        } else if( os.contains( "mac" ) || os.contains( "nix" ) || os.contains( "nux" ) || os.contains( "aix" ) ) {
            System.out.println( "running under MacOS/Unix/Linux" );
        } else {
            System.out.println( "Unrecognized OS: " + os);
	    Utils.stop( "Unrecognized OS", 0 );
	}

	// ***** READ system variable for the paths
	String myVarValue = System.getenv( "DATA_PATH" ); // path to the data directory
	if( myVarValue != null ) {
	    System.out.println("Using path: " + myVarValue );
	} else {
	    myVarValue = fileSystemSeparator + "Users" + fileSystemSeparator + "major" + fileSystemSeparator + "web" + fileSystemSeparator + "data" + fileSystemSeparator + "rimap-risc" + fileSystemSeparator;
	}

	String version = args[0] + fileSystemSeparator;

	// ***** READ results

	String gencodePath = myVarValue + "gencode." + version;
	System.out.println( gencodePath );

        // Create a File object
        File directory = new File( gencodePath );

        // Check if the directory exists and is indeed a directory
        if( directory.exists() && directory.isDirectory() ) {
            // Get an array of File objects representing the files and directories
            File[] filesList = directory.listFiles();

            // Loop through the array and print the names of files and directories
            if( filesList != null ) {
                for( File file : filesList ) {
                    if( file.isDirectory() ) {
			System.out.println( "Main directory: " + file.getName() );
			String subDirectoryPath = gencodePath + file.getName();
			System.out.println( "subDirectoryPath = " + subDirectoryPath );
                        File subDirectory = new File( subDirectoryPath );
			File[] pctList = subDirectory.listFiles();
			if( pctList != null ) {
			    for( File pctFile : pctList ) {
				String pctSubDirectoryPath = subDirectoryPath + fileSystemSeparator + pctFile.getName();
				File pctSubDirectory = new File( pctSubDirectoryPath );
				File[] subPctList = pctSubDirectory.listFiles();
				if( subPctList != null ) {
				    if( pctFile.isDirectory() )
					if( pctFile.getName().contains( ".scan" ) ) {
					    System.out.println( ".scan directory: " + pctFile.getName() );
					    File[] resultList = pctFile.listFiles();
					    if( resultList != null ) { // contains results
						boolean containsAccessibilityFile = false;
						boolean lockedFile = false;
						int numberOfDataFiles = 0;
						for( File resultFile: resultList )
						    if( resultFile.getName().equals( "accessibilities" ) ) {
							containsAccessibilityFile = true;
							break; // exit loop
						    }
						    else if( resultFile.getName().contains( ".lock" ) ) { // lock file present
							lockedFile = true;
							break; // exit loop (someone else is working on it)
						    }
						    else // one data file
							numberOfDataFiles++;
						if( !lockedFile && !containsAccessibilityFile && numberOfDataFiles > 0 ) { // must generate accessibilities
						    System.out.println( "accessibilities must be generated!" );
						    // generate the accessibilities for this transcript
						    System.out.println( "pctFile.getName() = " + pctFile.getName() );
						    String transcript = "";
						    if( pctFile.getName().endsWith( ".scan" ) ) // should be the case
							transcript = pctFile.getName().substring( 0, pctFile.getName().length() - ".scan".length() );
						    String accessibilityFile = subDirectoryPath + fileSystemSeparator + pctFile.getName() + fileSystemSeparator + "accessibilities";
						    String dataPath = System.getenv( "DATA_PATH" ); // path to the data directory (includes the final file separator)
						    RIMapRISCIO riMapRISCIO = new RIMapRISCIO( dataPath, args[0] ); // gencode version without the "/"
						    ProteinCodingTranscript target = riMapRISCIO.getTranscript( transcript );
						    String sequence = target.getStringSequence();
						    System.out.println( "createAccessFile( " + accessibilityFile + ", " + sequence );
						    RIMapRISCIO.createAccessibilityFile( accessibilityFile, sequence );
						}
						// else if( numberOfDataFiles == 0 )
						//     System.out.println( "no data file present" );
						// else System.out.println( "accessibilities available!" );
					    }
					}
				}
			    }
			}
		    }
		}
            } else {
                System.out.println( "The directory is empty or an I/O error occurred." );
            }
        } else {
            System.out.println( "The specified path is not a directory or does not exist." );
        }
    }
}
