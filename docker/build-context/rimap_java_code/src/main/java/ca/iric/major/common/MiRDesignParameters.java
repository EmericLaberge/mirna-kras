/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class MiRDesignParameters {

    Map<String,List<String>> parameters = new HashMap<>();

    public MiRDesignParameters( String filename ) {
	parameters = readDataFromFile( filename );
    }

    // accessors
    public List<String> get( String key ) {
	return this.parameters.get( key );
    }

    public static Map<String, List<String>> readDataFromFile( String filePath ) {
        Map<String, List<String>> dataMap = new HashMap<>();
        // Pattern adjusted to use positive lookahead to capture until the next keyword or end of input
        Pattern pattern = Pattern.compile( "(size|required|excluded|optional|gencode|region|seed-prefix|supp-prefix|gc-percent-min|gc-percent-max)[:,;\\s]*(.*?)(?=\\b(size|required|excluded|optional|gencode|region|seed-prefix|supp-prefix|gc-percent-min|gc-percent-max)\\b|$)", Pattern.DOTALL );

        try( BufferedReader reader = new BufferedReader( new FileReader( filePath ) ) ) {
            StringBuilder fileContent = new StringBuilder();
            String line;
            // Read the entire file content into a single StringBuilder
            while( ( line = reader.readLine() ) != null ) {
                fileContent.append( line ).append( "\n" );  // Preserve new lines for correct parsing
            }
            
            // Use Matcher to find all keyword-data pairs in the file content
            Matcher matcher = pattern.matcher( fileContent.toString() );
            while( matcher.find() ) {
                String key = matcher.group( 1 );
                String value = matcher.group( 2 ).trim();
                List<String> items = new ArrayList<>(Arrays.asList( value.split( "[,;\\s\n]+" ) ) );
                items.removeIf( String::isEmpty );  // Remove empty strings
                dataMap.put( key, items );
            }
        } catch( IOException e ) {
            e.printStackTrace();
        }

        return dataMap;
    }
}
