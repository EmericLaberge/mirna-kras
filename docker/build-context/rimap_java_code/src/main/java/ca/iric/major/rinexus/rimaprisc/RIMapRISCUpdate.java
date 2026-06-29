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
import java.text.SimpleDateFormat;
import java.util.Date;

public class RIMapRISCUpdate {

    public static String dropAfterDot( String input ) {
        int dotIndex = input.indexOf('.');
        if( dotIndex != -1 ) {
            return input.substring( 0, dotIndex );
        } else {
            return input; // No dot found, return the original string
        }
    }

    public static void main( String[] args ) {

	// ***** READ ARGUMENTS

	if( args.length == 0 || args.length > 2 ) {
	    System.out.println( "usage: MiRScanUpdate <Gencode version> [optional <date>]" );
	    System.out.println( "   Gencode version: ex) v46" );
	    System.out.println( "   date: 'erase' the results generated before date, in the format 'yyyy-MM-dd HH:mm:ss' (optional)" );
	    System.out.println( "          if not provided, the list of entries will be output" );
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

	SimpleDateFormat sdf = new SimpleDateFormat( "yyyy-MM-dd HH:mm:ss" );
	Date date = null;
	long dateMillis = 0;
	String inputDateString = "";
	String version = args[0] + fileSystemSeparator;
	if( args.length > 1 ) {
	    inputDateString = args[1];
	    try {
		// Parse the input date string into a Date object
		date = sdf.parse( inputDateString );
	    } catch( ParseException e ) {
		System.out.println( "Error: Unable to parse the input date. Please use the format 'yyyy-MM-dd HH:mm:ss'." );
		e.printStackTrace();
	    }
	    // Convert the Date object to milliseconds since the epoch
	    dateMillis = date.getTime();
	    String ppDate = sdf.format( date );
	    System.out.println( "Delete results generated before " + ppDate );
	}
	else
	    System.out.println( "List results" );


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
			String subDirectoryPath = gencodePath + fileSystemSeparator + file.getName();
                        File subDirectory = new File( subDirectoryPath );
			File[] pctList = subDirectory.listFiles();
			if( pctList != null ) {
			    for( File pctFile : pctList ) {
				String pctSubDirectoryPath = subDirectoryPath + fileSystemSeparator + pctFile.getName();
				File pctSubDirectory = new File( pctSubDirectoryPath );
				File[] subPctList = pctSubDirectory.listFiles();
				if( subPctList != null ) {
				    for( File subPctFile: subPctList ) {
					if( subPctFile.isFile()) {
					    // Get the last modified time in milliseconds
					    long lastModified = subPctFile.lastModified();

					    // Convert the milliseconds to a Date object
					    Date lastModifiedDate = new Date( lastModified );

					    // Format the date to a readable format
					    String formattedDate = sdf.format( lastModifiedDate );
					    System.out.print("File: " + file.getName() + fileSystemSeparator + pctFile.getName() + fileSystemSeparator + subPctFile.getName() + ", last modified: " + formattedDate );
					    if( dateMillis > 0 )
						if( lastModified < dateMillis ) {
						    subPctFile.delete();
						    System.out.println( " removed" );
						}
						else System.out.println( " kept" );
					    else System.out.println( "" );
					} else if( subPctFile.isDirectory() ) {
					    System.out.println( "Directory: " + subPctFile.getName() );
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
