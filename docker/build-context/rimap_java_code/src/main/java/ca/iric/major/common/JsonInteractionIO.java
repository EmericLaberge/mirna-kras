/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.module.SimpleModule;

import java.io.*;
import java.nio.file.*;
import java.util.List;

import java.util.stream.Collectors;

/*
 * interface a list of interactions saved in json formatted files (list of JsonInteraction)
 *     needs the path to the file
 */

public class JsonInteractionIO {

    private static final ObjectMapper mapper = new ObjectMapper();

    public static List<JsonInteraction> readJsonFile( Path path ) throws IOException {
	ObjectMapper mapper = new ObjectMapper();
	mapper.configure( DeserializationFeature.FAIL_ON_UNKNOWN_PROPERTIES, false );
	try( BufferedReader reader = Files.newBufferedReader( path ) ) {
	    return reader.lines()
                .map( line -> {
			try {
			    return mapper.readValue( line, JsonInteraction.class );
			} catch( IOException e ) {
			    throw new UncheckedIOException( "Failed to parse line: " + line, e );
			}
		    })
                .collect( Collectors.toList() );
	}
    }

    public static List<JsonInteraction> readJsonFile( Path path, long maxKd, int region ) throws IOException {
	ObjectMapper mapper = new ObjectMapper();
	try( BufferedReader reader = Files.newBufferedReader( path ) ) {
	    return reader.lines()
                .map( line -> {
			try {
			    return mapper.readValue( line, JsonInteraction.class );
			} catch( IOException e ) {
			    throw new UncheckedIOException( "Failed to parse line: " + line, e );
			}
		    })
                .filter( i -> i.getKd() <= maxKd && i.inRegion( region ) )
                .collect( Collectors.toList() );
	}
    }

    public static List<JsonInteraction> readJsonFile( Path path, long maxKd ) throws IOException {
	ObjectMapper mapper = new ObjectMapper();
	try( BufferedReader reader = Files.newBufferedReader( path ) ) {
	    return reader.lines()
                .map( line -> {
			try {
			    return mapper.readValue( line, JsonInteraction.class );
			} catch( IOException e ) {
			    throw new UncheckedIOException( "Failed to parse line: " + line, e );
			}
		    })
                .filter( i -> i.getKd() <= maxKd )
                .collect( Collectors.toList() );
	}
    }

    public static void writeJsonFile(Path path, List<JsonInteraction> interactions) throws IOException {

	ObjectMapper mapper = new ObjectMapper();

	// keep BufferedWriter open across calls
	mapper.getFactory().disable(JsonGenerator.Feature.AUTO_CLOSE_TARGET);

	// force 1 decimal for all doubles
	SimpleModule module = new SimpleModule();
	module.addSerializer(Double.class, new TwoDecimalDoubleSerializer());
	module.addSerializer(double.class, new TwoDecimalDoubleSerializer());
	mapper.registerModule(module);

	try (BufferedWriter writer = Files.newBufferedWriter(path)) {
	    for (JsonInteraction interaction : interactions) {
		mapper.writeValue(writer, interaction);
		writer.newLine();
	    }
	}
    }
}
