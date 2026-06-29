/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.nio.file.*;
import java.nio.charset.StandardCharsets;

// manage the I/O of guide sequences from 20 to 26 nucleotides
//     assume the information is stored in Json format
public class GuideStorageManager {

    private static final String ENV_VAR_NAME = "GUIDES_PATH"; // $GUIDES_PATH must be defined in the running shell
    private static final int MIN_GUIDE_LENGTH = 18; // some examples in miRBase are 16-nt long (18 said Fromm)
    private static final int MAX_GUIDE_LENGTH = 26; // some examples in miRBase are 28-nt long (26 said Fromm)
    private final Path root; // path to the guides root directory

    /**
     * use system environment variable $GUIDES_PATH.
     */
    public GuideStorageManager() {
        String envPath = System.getenv( ENV_VAR_NAME );
        if( envPath == null || envPath.isBlank() ) {
            throw new IllegalStateException( "Environment variable $" + ENV_VAR_NAME + " is not set." );
        }
        this.root = Paths.get( envPath );
    }

    // constructor by root directory
    public GuideStorageManager( String rootDirectory ) {
        this.root = Paths.get( rootDirectory );
    }

    // constructor by path
    public GuideStorageManager( Path rootDirectory ) {
        this.root = rootDirectory;
    }

    public Path getGuidePath( String guideSequence ) {
        if( guideSequence.length() < MIN_GUIDE_LENGTH || guideSequence.length() > MAX_GUIDE_LENGTH )
            throw new IllegalArgumentException( "Guide sequence must be between 20 and 26 nucleotides long, " + guideSequence );

        String g1 = guideSequence.substring( 0, 1 ); // g1
        String seed = guideSequence.substring( 1, 8 ); // seed
        String filename = guideSequence + ".json"; // file name is the guide sequence with extension json
        return root.resolve( Paths.get( g1, seed, filename ) );
    }

    public void saveGuide( String guideSequence, String guideData ) throws IOException {
        Path path = getGuidePath( guideSequence );
        Files.createDirectories( path.getParent() );
        Files.writeString( path, guideData, StandardCharsets.UTF_8, StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING );
    }

    public String loadGuide( String guideSequence ) throws IOException {
        Path path = getGuidePath( guideSequence );
        return Files.readString( path, StandardCharsets.UTF_8 );
    }

    public boolean guideExists( String guideSequence ) {
        return Files.exists( this.getGuidePath( guideSequence ) );
    }
}
