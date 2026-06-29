/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import org.springframework.web.bind.annotation.*;
import org.springframework.http.ResponseEntity;
import org.springframework.http.MediaType;
import org.springframework.http.HttpStatus;

import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.*;
import java.util.*;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

@CrossOrigin( origins = "http://localhost:3000" ) // uncomment for local development
//@CrossOrigin( origins = "https://rimap-risc-backend.major.iric.ca", allowedHeaders = "*", methods = {RequestMethod.GET, RequestMethod.POST, RequestMethod.OPTIONS}) // uncomment for global development
@RestController
@RequestMapping( "/api" )
public class RISCScanController {

    private String miRNATargetPair = "";
    private String bindingLandscape = "";
    private List<Integer> transcriptRegions = new ArrayList<>();
    private List<Integer> bindingPositions = new ArrayList<>();
    private List<Double> transcriptAccessibility = new ArrayList<>();

    public static List<Integer> extractIntegersFromLine( String line ) {
        List<Integer> integers = new ArrayList<>();

        // Debug: Print the input line
        //System.out.println("Input line: " + line);

        // Regular expression to match integers
        Pattern pattern = Pattern.compile("\\d+");
	if( line == null ) return integers;
        Matcher matcher = pattern.matcher(line);

        // Debug: Check if matcher is created correctly
        //System.out.println("Pattern: " + pattern.pattern());

        // Find all matches
        while( matcher.find() ) {
            // Debug: Log each match found
            //System.out.println("Matched: " + matcher.group());

            // Parse the integer and add it to the list
            integers.add( Integer.parseInt(matcher.group() ) );
        }

        // Debug: Print the extracted integers
        //System.out.println("Extracted integers: " + integers);

        return integers;
    }

        public static List<Long> extractLongsFromLine( String line ) {
        List<Long> longs = new ArrayList<>();

        // Debug: Print the input line
        //System.out.println("Input line: " + line);

        // Regular expression to match integers
        Pattern pattern = Pattern.compile("\\d+");
        Matcher matcher = pattern.matcher(line);

        // Debug: Check if matcher is created correctly
        //System.out.println("Pattern: " + pattern.pattern());

        // Find all matches
        while( matcher.find() ) {
            // Debug: Log each match found
            //System.out.println("Matched: " + matcher.group());

            // Parse the integer and add it to the list
            longs.add( Long.parseLong(matcher.group() ) );
        }

        // Debug: Print the extracted integers
        //System.out.println("Extracted integers: " + integers);

        return longs;
    }

    public static List<Double> extractDoublesFromLine(String line) {
        List<Double> doubles = new ArrayList<>();

        // Debug: Print the input line
        //System.out.println("Input line: " + line);

        // Regular expression to match integers
        Pattern pattern = Pattern.compile("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?");
        Matcher matcher = pattern.matcher(line);

        // Debug: Check if matcher is created correctly
        //System.out.println("Pattern: " + pattern.pattern());

        // Find all matches
        while( matcher.find() ) {
            // Debug: Log each match found
            //System.out.println("Matched: " + matcher.group());

            // Parse the integer and add it to the list
            doubles.add( Double.parseDouble(matcher.group() ) );
        }

        // Debug: Print the extracted integers
        //System.out.println("Extracted integers: " + integers);

        return doubles;
    }

    @PostMapping( "/data" )
    public ResponseEntity<Map<String,Object>> process( @RequestBody RISCScanRequest request ) {
	//try {
	    //System.out.println("DEBUG: Received request = " + new ObjectMapper().writeValueAsString( request ) );
	//} catch( com.fasterxml.jackson.core.JsonProcessingException e ) {
	//   e.printStackTrace();
	//    return ResponseEntity.status( HttpStatus.INTERNAL_SERVER_ERROR )
	//	.body( Map.of( "error", "Internal processing error", "details", e.getMessage() ) );
	//}
        this.bindingLandscape = callRISCScanProgram(
					       request.getGencodeVersion(),
					       request.getTranscriptName(),
					       request.getMiRName(),
					       request.getMinSeed(),
					       request.getMaxKd(),
					       request.isUTR5() ? "true" : "false",
					       request.isCDS() ? "true" : "false",
					       request.isUTR3() ? "true" : "false",
					       request.isSupplementaryPaired() ? "true" : "false",
					       request.getSeedAccessible(),
					       request.getSupplementaryAccessible(),
					       request.getConserved() );

        // Prepare the response
        Map<String, Object> response = new HashMap<>();
	response.put( "miRNATargetPair", this.miRNATargetPair );
        response.put( "bindingLandscape", this.bindingLandscape );
        response.put( "transcriptRegions", this.transcriptRegions );
        response.put( "bindingPositions", this.bindingPositions );
	response.put( "transcriptAccessibility", this.transcriptAccessibility );

	return ResponseEntity.ok()
            .contentType( MediaType.APPLICATION_JSON )
            .body( response );
    }

    private String callRISCScanProgram(
				      String gencodeVersion,
				      String transcriptName,
				      String miRName,
				      String minSeed,
				      double maxKd,
				      String UTR5,
				      String CDS,
				      String UTR3,
				      String supplementaryPaired,
				      double seedAccessible,
				      double SupplementaryAccessible,
				      int conserved ) {
        try {
            // Define the path to the JAR and the working directory
	    String workingDirPath = System.getenv( "WORKING_DIR" );
	    if( workingDirPath == null ) {
		throw new IllegalStateException( "Environment variable WORKING_DIR is not set." );
	    }
	    File workingDir = new File( workingDirPath );
	    if( gencodeVersion == null ) gencodeVersion = "v46";
            if( transcriptName == null ) transcriptName = "SIRT1-201";
	    if( miRName == null ) miRName = "miR-34a-5p";
	    if( minSeed == null ) minSeed = "6mer";
	    if( maxKd == 0 ) maxKd = 1500;
            ProcessBuilder processBuilder = new ProcessBuilder(
							       "java",
							       "-Xms32g",
							       "-Xmx32g",
							       "-jar",
							       "RISCScan.jar",
							       gencodeVersion,
							       transcriptName,
							       miRName,
							       minSeed,
							       String.valueOf( maxKd ),
							       UTR5,
							       CDS,
							       UTR3,
							       supplementaryPaired,
							       String.valueOf( seedAccessible ),
							       String.valueOf( SupplementaryAccessible ),
							       String.valueOf( conserved ) );

            // Set the working directory
	    //processBuilder.redirectErrorStream( true ); // combine stderr with stdout if needed
            processBuilder.directory( workingDir );

            // Start the process
            Process process = processBuilder.start();

            // Read the output from the process
            BufferedReader reader = new BufferedReader( new InputStreamReader( process.getInputStream() ) );
            String line = "initString";
	    
	    // use the following loop to debut the input lines
	    // while( ( line = reader.readLine() ) != null ) {
	    // 	System.out.println( "OUT: " + line );  // 🔍 Log every output line
	    // 	if( line.trim().isEmpty() ) {
	    // 	    System.out.println( "⚠️ Skipping empty line" );
	    // 	    continue;
	    // 	}
	    // 	extractIntegersFromLine( line );  // will no longer crash if you added null check
	    // }
	    
            StringBuilder output = new StringBuilder();

	    // first line contains target:guide id
	    // second line contains transcript regions
	    // middle lines contain the binding landscape
	    // last line contains interaction positions (with explicit "interaction positions: ..." string)

	    line = reader.readLine(); // read first line
	    this.miRNATargetPair = line.trim().strip();;
	    line = reader.readLine(); // read second line
	    this.transcriptRegions = extractIntegersFromLine( line );

	    // the binding landscape; will be returned by function
	    int lineNo = 0;
	    String lastLine = "";
	    this.bindingPositions.clear();
	    this.transcriptAccessibility.clear();
            while( ( line = reader.readLine() ) != null ) {
		if( !line.contains( "controler" ) ) {
		    output.append( line ).append( "\n" );
		    if( ++lineNo % 3 == 0 ) output.append( "\n" );
		}
		// get the positions, line contains "interaction positions:"
		else if( line.contains( "positions" ) )	{
		    //System.out.println( "positions read" );
		    this.bindingPositions = extractIntegersFromLine( line ); // returns the List of positions
		}
		// get the accessibilities, line contains "accessibilities:"
		else if( line.contains( "accessibilities" ) ) {
		    //System.out.println( "accessibilities read" );
		    this.transcriptAccessibility = extractDoublesFromLine( line ); // returns the List of accessibilities
		}
            }

	    // Wait for the process to complete
            int exitCode = process.waitFor();
	    switch( exitCode ) {
	    case 0:
                this.bindingLandscape = output.toString(); // normal termination => process data
		break;
	    // treat abnormal termination codes for specific user message
	    case 1:
		this.bindingLandscape = "Illegal Gencode version.";
		break;  
	    case 2:
		this.bindingLandscape = "MicroRNA name is not a valid miRBase microRNA identifier.";
		break;
	    case 3:
		this.bindingLandscape = "Transcript not found in the protein-coding set of Gencode v46.";
		break;
	    case 4:
		this.bindingLandscape = "Please, provide valid miRBase microRNA name and gencode.v46 transcript name.";
		break;
	    case 5:
		this.bindingLandscape = "Error running bedtools.";
		break;
	    case 6:
		this.bindingLandscape = "I/O problem with transcript.";
		break;
	    case 7:
		this.bindingLandscape = "I/O problem with microRNA.";
		break;
	    case 8:
		this.bindingLandscape = "Accessibility file not found.";
		break;
	    case 9:
		this.bindingLandscape = "General I/O file system error.";
		break;
	    case 10:
		this.bindingLandscape = "Problem with file closure.";
		break;
	    case 11:
		this.bindingLandscape = "Accessibility file expected, but not found.";
		break;
	    case 12:
		this.bindingLandscape = "Error in folding: illegal symbol in dot-bracket.";
		break;
	    case 13:
		this.bindingLandscape = "Error in folding: unbalanced dot-bracket.";
		break;
	    case 14:
		this.bindingLandscape = "Error in dot-bracket building.";
		break;
	    case 15:
		this.bindingLandscape = "Error in duplex folding.";
		break;
	    case 16:
		this.bindingLandscape = "Error in writing Json guide.";
		break;
	    case 17:
		this.bindingLandscape = "Unix DATA_PATH variable not set.";
		break;
	    case 18:
		this.bindingLandscape = "IO error writing to accessibility file.";
		break;
	    case 19:
		this.bindingLandscape = "IO file system problem.";
		break;
	    case 20:
		this.bindingLandscape = "IO error writing interaction file";
		break;
	    case 21:
		this.bindingLandscape = "IO error deleting lock file.";
		break;
	    default:
		{
		    // BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		    // String errorLine;
		    // while ((errorLine = errorReader.readLine()) != null) {
		    //     System.out.println(errorLine);
		    // }
		    this.bindingLandscape = "Internal error, please contact the developer: francois.major@umontreal.ca; provide him with the transcript and microRNA names and this exit code: " + exitCode;
		}
	    }
	    
        } catch( IOException | InterruptedException e ) {
            e.printStackTrace();
            this.bindingLandscape = "Internal error, please contact the developer: francois.major@umontreal.ca; provide him with the transcript and microRNA names and this message: " + e.getMessage();
        }
        return this.bindingLandscape;
    }
}
