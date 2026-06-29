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

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.*;
import com.fasterxml.jackson.databind.node.ObjectNode;

public class RIMapRISCUpdateTarBase {

    public static String dropAfterDot( String input ) {
        int dotIndex = input.indexOf('.');
        if( dotIndex != -1 ) {
            return input.substring( 0, dotIndex );
        } else {
            return input; // No dot found, return the original string
        }
    }

    private static final ObjectMapper MAPPER = new ObjectMapper();

    public static void updateNdjsonFile( File file, Set<Integer> posSet ) throws IOException {
        // Reader for one object per line
        ObjectReader reader = MAPPER.readerFor( JsonInteraction.class );

        // Create a temp file to write updated NDJSON
        Path tmp = Files.createTempFile( file.toPath().getParent(), "itb-", ".tmp" );

        try (
            InputStream in = Files.newInputStream( file.toPath() );
            BufferedWriter out = Files.newBufferedWriter( tmp )
        ) {
            MappingIterator<JsonInteraction> it = reader.readValues( in );
            boolean modifiedAny = false;

            while( it.hasNext() ) {
                JsonInteraction ji = it.next();
                // Set itb=true if pst is in TarBase set
                if( posSet.contains( ji.getPositionStart() ) && !ji.getInTarBase() ) {
                    ji.setInTarBase( true );
                    modifiedAny = true;
                }
                // Write one JSON object per line
                out.write( MAPPER.writeValueAsString( ji ) );
                out.newLine();
            }
            out.flush();
            if( modifiedAny ) {
                // Atomically replace original with updated temp file
                Files.move( tmp, file.toPath(), StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE );
                System.out.println( "Updated: " + file.getPath() );
            } else {
                Files.deleteIfExists( tmp );
            }
        }
    }

    public static void main( String[] args ) throws IOException {

	// ***** READ ARGUMENTS

	if( args.length < 2 || args.length > 2 ) {
	    System.out.println( "usage: MiRScanUpdateTarBase <Gencode version> <TarBaseSupportedSites>" );
	    System.out.println( "   Gencode version: ex) v46" );
	    System.out.println( "   TarBaseSupportedSites: file that contains the TarBase sites as a Map<String,Map<String, List<Integer>>>" );
	    Utils.stop( "bye!" , 0 );
	}

	char fileSystemSeparator = '/'; // default separator
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
	    System.out.println( "Environment variable DATA_PATH is not set." );
	    System.exit( 0 );
	}

	String version = args[0] + fileSystemSeparator;

	// ***** READ results

	String gencodePath = myVarValue + "gencode." + version;
	System.out.println( gencodePath );

	// ***** READ TarBase supported sites
	
	Path tarbaseJson = Path.of( args[1] );
	// structure { "UBR7-201": { "UGAGG...": [1167, 197, ...], ... }, ... }
	Map<String, Map<String, List<Integer>>> tarbase =
	    MAPPER.readValue( Files.newBufferedReader( tarbaseJson ),
			      new TypeReference<>() {} );

        // Create a Path to the gencode library
        Path directory = Paths.get(  gencodePath );

        // Check if the directory exists and is indeed a directory
        if( Files.exists( directory ) && Files.isDirectory( directory ) ) {
	    // loop thru the tarbase and update the json data file to set "itb" to true
	    for( var transcriptEntry : tarbase.entrySet() ) {
		String transcript = transcriptEntry.getKey();
		Map<String, List<Integer>> mirnaToPositions = transcriptEntry.getValue();

		// Derive subdir = first two letters of transcript
		String subdir = transcript.substring( 0, Math.min( 2, transcript.length() ) );

		for( var mirnaEntry : mirnaToPositions.entrySet() ) {
		    String sequence = mirnaEntry.getKey();
		    List<Integer> positions = mirnaEntry.getValue();

		    if( positions != null && !positions.isEmpty() ) {			
			//String path = gencodePath + subdir + "/" + transcript + ".scan/" + sequence;
			Path path = RIMapRISCIO.getInteractionsFilePath( directory, transcript, sequence );

			// get the corresponding file
			if( !Files.exists( path ) ) {
			    System.out.println( "file " + path + " does not exist!" );
			    continue;
			}

			// file exists
			File file = path.toFile();
			Set<Integer> posSet = new HashSet<>( positions );
			System.out.println( "updating " + file + ", " + posSet );
			updateNdjsonFile( file, posSet );
		    }
		}
	    }
        } else {
            System.out.println( "The specified path is not a directory or does not exist, " + directory + "." );
        }
    }
}
