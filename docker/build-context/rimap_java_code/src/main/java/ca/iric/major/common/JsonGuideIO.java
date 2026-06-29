/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.*;
import java.nio.file.*;
import java.util.List;

/*
 * interface guides saved in json formatted files (individual JsonGuide)
 *     needs the path to the file
 */

public class JsonGuideIO {

    private static final ObjectMapper mapper = new ObjectMapper();

    public static JsonGuide readJsonFile( Path path ) throws IOException {
        try( InputStream in = Files.newInputStream( path ) ) {
            return mapper.readValue( in, new TypeReference<JsonGuide>() {} );
        }
    }

    public static void writeJsonFile( Path path, JsonGuide guide ) throws IOException {
        try( OutputStream out = Files.newOutputStream( path ) ) {
            mapper.writerWithDefaultPrettyPrinter().writeValue( out, guide );
        }
    }
}
