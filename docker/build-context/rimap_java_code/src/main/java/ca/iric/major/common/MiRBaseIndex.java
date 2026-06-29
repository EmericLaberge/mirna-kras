package ca.iric.major.common;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MiRBaseIndex implements Serializable {

    private static final long serialVersionUID = 1L;

    // Map: mature sequence -> miRNA name
    private final Map<String,String> seqToName;

    public MiRBaseIndex( Map<String, String> seqToName ) {
        this.seqToName = seqToName;
    }

    public Map<String,String> getSeqToNameMap() {
        return seqToName;
    }

    /**
     * Build the index from a miRBase directory containing "let" and "mir" subdirectories.
     *
     * miRBase/
     *   let/
     *     let-7a-2-3p
     *   mir/
     *     miR-34a-5p
     */
    public static MiRBaseIndex buildFromMiRBaseDir( Path miRBaseDir ) throws IOException {
        Map<String, String> map = new HashMap<>();

        // Process both "let" and "mir" subdirectories if they exist
        processSubDir( miRBaseDir.resolve( "let" ), map);
        processSubDir( miRBaseDir.resolve( "mir" ), map);

        return new MiRBaseIndex( map );
    }

    private static void processSubDir( Path dir, Map<String,String> map) throws IOException {
        if( !Files.exists( dir ) || !Files.isDirectory( dir ) ) {
            return; // silently skip if missing
        }

        try( DirectoryStream<Path> stream = Files.newDirectoryStream( dir ) ) {
            for( Path file : stream ) {
                if( Files.isDirectory( file ) ) {
                    continue;
                }

                String name = stripExtension( file.getFileName().toString() ); // e.g. "let-7a-2-3p"

                // Read FASTA file and extract sequence (2nd non-empty line or concatenated sequence lines)
                String seq = readFastaSequence( file );
                if( seq == null || seq.isEmpty() ) {
                    System.err.println( "Warning: no sequence found in " + file );
                    continue;
                }

                // Normalize sequence (optional: upper case, trim)
                seq = seq.trim().toUpperCase();

                // Put into map: sequence -> name
                String previous = map.put( seq, name );
                if( previous != null && !previous.equals( name ) ) {
                    System.err.println( "Warning: sequence collision for " + seq +
                            " (" + previous + " vs " + name + ")" );
                }
            }
        }
    }

    /**
     * Reads the FASTA file and returns the sequence.
     * Assumes:
     *   - first non-empty line starting with '>' is the header
     *   - subsequent non-empty lines are sequence lines (joined)
     */
    private static String readFastaSequence( Path file ) throws IOException {
        List<String> lines = Files.readAllLines( file, StandardCharsets.UTF_8 );

        StringBuilder sb = new StringBuilder();
        boolean inSeq = false;

        for( String line : lines ) {
            line = line.trim();
            if( line.isEmpty() ) {
                continue;
            }
            if( line.startsWith( ">" ) ) {
                // header line; start reading sequence after this
                inSeq = true;
                continue;
            }
            if( inSeq ) {
                sb.append( line );
            }
        }

        return sb.toString();
    }

    private static String stripExtension( String fileName ) {
        int dot = fileName.lastIndexOf( '.' );
        return (dot == -1) ? fileName : fileName.substring( 0, dot );
    }

    // ---------- Serialization helpers ----------

    /**
     * Serialize this index to a single file on disk.
     */
    public void saveToFile( Path file ) throws IOException {
        Files.createDirectories( file.getParent() );
        try( ObjectOutputStream oos = new ObjectOutputStream( Files.newOutputStream( file ) ) ) {
            oos.writeObject( this.seqToName );
        }
    }

    /**
     * Load a Map<String, String> (sequence -> name) from a serialized file.
     */
    @SuppressWarnings( "unchecked" )
    public static MiRBaseIndex loadFromFile( Path file ) throws IOException {
        try( ObjectInputStream ois = new ObjectInputStream( Files.newInputStream( file ) ) ) {
            Object obj = ois.readObject();
            if( !( obj instanceof Map ) ) {
                throw new IOException( "Invalid serialized miRBase index file: " + file );
            }
            return new MiRBaseIndex( ( Map<String,String>) obj );
        } catch( ClassNotFoundException e ) {
            throw new IOException( "Class not found while reading miRBase index file: " + file, e );
        }
    }

    // Convenience lookup
    public String getNameForSequence( String sequence ) {
        return seqToName.get( sequence );
    }
}
