/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import ca.iric.major.common.*;

import com.fasterxml.jackson.databind.ObjectMapper;

import java.io.IOException;
import java.io.FileNotFoundException;
import java.io.File;
import java.nio.file.*;
import java.util.List;
import java.util.LinkedList;

public final class RIMapRISCIO {

    public static final String MIRBASEDIRNAME = "MiRBase";

    static final int g1 = 0; // index of g1 in a guide
    static final int g2 = 1; // index of g2 in a guide
    static final int g8 = 7; // index of g8 in a guide
    static final int g9 = 8; // index of g9 in a guide
    static final int g12 = 11; // index of g12 in a guide
    static final int g13 = 12; // index of g13 in a guide
    static final int g17 = 16; // index of g17 in a guide
    static final int g18 = 17; // index of g18 in a guide
    static final int seedStart = g2;
    static final int seedEnd = g8;
    static final int centralStart = g9;
    static final int centralEnd = g12;
    static final int suppStart = g13;
    static final int suppEnd = g17;
    static final int tailStart = g18;

    private static final int SEED_START = 1;      // g2
    private static final int SEED_END   = 7;      // g8
    private static final int SUPP_START = 12;     // g13
    private static final int SUPP_END   = 16;     // g17
    private static final int G1         = 0;
    private static final int CENTRAL_START = 8;   // g9
    private static final int CENTRAL_END   = 11;  // g12
    
    private Path dataPath;
    private String gencodeVersion;
    private DisturbanceManager disturbanceManager;
    private Path guidesPath;
    private Path interactionsPath;
    private Path miRBasePath;

    public RIMapRISCIO( String dataPath, String gencodeVersion ) {
	this.dataPath = Paths.get( dataPath ); // includes the separator
	this.gencodeVersion = gencodeVersion;
	this.guidesPath = Paths.get( dataPath + File.separator + "guides" + File.separator );
	this.miRBasePath = Paths.get( dataPath + "MiRBase" + File.separator );
	this.interactionsPath = Paths.get( dataPath + File.separator + "gencode." + gencodeVersion + File.separator );
	// check the directory for this gencode version exists
	if( Files.exists( this.interactionsPath ) && Files.isDirectory( this.interactionsPath ) ) {
	    // Directory exists, do nothing
	} else {
	    Utils.stop( "No gencode directory matching this version: " + this.interactionsPath.toString(), 1 );
	}
	// the Path.toString() method drops the last File.separator!
	this.disturbanceManager = new DisturbanceManager( this.dataPath.toString() + File.separator + "gencode." + gencodeVersion + ".pc_transcripts.disturbance" );
    }

    // getters
    public Path getDataPath() { return this.dataPath; }
    public String getGencodeVersion() { return this.gencodeVersion; }
    public Path getGuidesPath() { return this.guidesPath; }
    public Path getInteractionsPath() { return this.interactionsPath; }
    public String getInteractionsDir() { return this.interactionsPath.toString(); }
    public Path getMiRBasePath() { return this.miRBasePath; }

    // if microRNA name, get the sequence from MiRBase/ and use a sync name, else use the guide, as file name
    public MatureMiRNA getGuide( String guide ) {
	MatureMiRNA microRNA = null;
	// check and transform guide into its sequence
	if( StringSequence.isRNA( guide.toUpperCase().replace( 'T', 'U' ) ) ) { // guide is an RNA sequence
	    guide = guide.toUpperCase().replace( 'T', 'U' ); // uppercase the sequence
	    // check is the sequence matches a microRNA in miRBase
	    MiRBaseIndex loaded;
	    String name = "";
	    try {
		loaded = MiRBaseIndex.loadFromFile( Paths.get( this.dataPath + File.separator + "miRBase-index.ser" ) ); // load the microRNA index (by mature sequence)
		name = loaded.getNameForSequence( guide );
	    } catch( Exception e ) {
		Utils.stop( "Error reading miRBase index, " + this.dataPath + File.separator + "miRBase-index.ser", 7 );
	    }
	    if( name == null ) { // we don't have a match
		// create a synthetic microRNA
		microRNA = new MatureMiRNA( new StringSequence( guide ) ); // we don't know yet the sequence
		microRNA.setFullName( SyncNameCodec.encodeSyncName( guide ) ); // since the guide was a sequence, encode a sync name
		microRNA.setMIMATCode( "" ); // no MIMAT for a sequence
		microRNA.setOrganism( "synthetic" ); // synthetic RNA
		microRNA.setName( guide ); // same as full name, which is the sequence
		return microRNA;
	    }
	    else guide = name; // assign the corresponding microRNA name to guide
	}
	// if we get here, guide is a microRNA or syncRNA name
	if( guide.contains( "sync" ) ) { // we have a syncRNA name
	    String sequence = SyncNameCodec.decodeSyncName( guide );
	    // check if the sequence is a natural microRNA (in miRBase)
	    MiRBaseIndex loaded;
	    String name = "";
	    try {
		loaded = MiRBaseIndex.loadFromFile( Paths.get( this.dataPath + File.separator + "miRBase-index.ser" ) ); // load the microRNA index (by mature sequence)
		name = loaded.getNameForSequence( sequence );
	    } catch( Exception e ) {
		Utils.stop( "Error reading miRBase index, " + this.dataPath + File.separator + "miRBase-index.ser", 7 );
	    }
	    if( name == null ) { // we don't have a match, create synthetic
		microRNA = new MatureMiRNA( new StringSequence( sequence ) );
		microRNA.setFullName( guide ); // since the guide was already a sync name
		microRNA.setMIMATCode( "" ); // no MIMAT for a sync
		microRNA.setOrganism( "synthetic" ); // synthetic RNA
		microRNA.setName( guide ); // same as full name, which is the sync name
		return microRNA;
	    }
	    else guide = name; // assign microRNA name to guide
	}
	// if we get here, we have a microRNA name in guide, so get the data from miRBase
	String subDir = guide.split( "-" )[0].toLowerCase(); // either "mir" or "let" (careful if hsa-mir-number is provided will crash)
	String fileName = this.miRBasePath + File.separator + subDir + File.separator + guide; // Path.toString() method drops the last separator
	List<MatureMiRNA> list = new LinkedList<>();
	try {
	    MiRBaseMature miRNAReader = new MiRBaseMature( fileName, list );
	    microRNA = list.getFirst();
	    guide = microRNA.getSequence().getSequence(); // get the sequence from miRBase microRNA
	} catch( IOException e ) {
	    Utils.stop( "Error reading microRNA in " + fileName, 7 );
	}
	// guide now has the correct sequence and name (which is the sequence)
	
	// String subDirL1 = guide.substring( 0, 1 ); // RIMap stores guides under four directories using the first nt in the guide
	// String subDirL2 = guide.substring( 1, 8 ); // RIMap stores guides under seed subdirectories
	// Path guideFilePath = Paths.get( this.guidesPath + File.separator + subDirL1 + File.separator + subDirL2 + File.separator + guide + ".json" ); // file name composed of full sequence + .json
	
	// if guide does not exist, create a new entry with microRNA data
	// if( !Files.exists( guideFilePath ) ) { // file does not exist, so create it
	//     // Serialize to JSON
	//     JsonGuide jsonGuide = new JsonGuide();
	//     jsonGuide.setCommonName( microRNA.getName() );
	//     jsonGuide.setId( microRNA.getMIMATCode() );
	//     jsonGuide.setSpecies( microRNA.getOrganism() );
	//     jsonGuide.setSequence( microRNA.getStringSequence() );
	//     jsonGuide.setDisturbance( this.disturbanceManager.getDisturbance( guide ) );
	//     ObjectMapper localMapper = new ObjectMapper();
	//     try {
	// 	String guideData = localMapper.writeValueAsString( jsonGuide );
	// 	this.gsm.saveGuide( guide, guideData );
	//     } catch( IOException e ) {
	// 	Utils.stop( "Problem converting or writing Json guide.", 16 );
	//     }
	// }
	
	// microRNA contains the miRBase microRNA or synthetic sequence
	return microRNA;
    }

    public Path getTranscriptFilePath( String transcript ) {
	// interactionsPath + two first letters of the transcript name + transcript name
	return Paths.get( this.interactionsPath + transcript.substring( 0, 2 ) + File.separator + transcript );
    }

    public String getTranscriptFile( String transcript ) {
	// interactionsPath + two first letters of the transcript name + transcript name; Path.toString() method drops the last separator!
	return this.interactionsPath + File.separator + transcript.substring( 0, 2 ) + File.separator + transcript;
    }

    public ProteinCodingTranscript getTranscript( String pctName ) {
	String pctAccess = getTranscriptFile( pctName );
	GencodePCTranscript target = null;
	// read transcript
	try {
	    target = new GencodePCTranscript( pctAccess );
	} catch( FileNotFoundException e ) {
	    Utils.stop( "Illegal gencode." + gencodeVersion + " transcript variant identifier:" + pctAccess, 3 );
	} catch( IOException e ) {
	    Utils.stop( "I/O problem with " + pctAccess, 6 );
	}
	return target.get( 0 );
    }

    // public Path getInteractionsFilePath( String transcript, String guide ) {
    // 	String subDirL1 = transcript.substring( 0, 2 ); // sub directory given by the two first letters of the transcript name
    // 	String subDirL2 = transcript + ".scan"; // sub directory for all interaction sets of a transcript is transcript name + extension .scan
    // 	return Paths.get( this.interactionsPath + subDirL1 + File.separator + subDirL2 + File.separator + guide ); // file name is the guide sequence
    // }

    public static Path getInteractionsFilePath(
					       Path interactionsPath,
					       String transcript,
					       String guide
					       ) {
	// First-level directory from transcript (e.g. "SI" for SIRT1-201)
	String subDirL1 = transcript.substring(0, 2);

	// Second-level directory: transcript.scan
	String subDirL2 = transcript + ".scan";

	// Extract guide segments
	String seed    = guide.substring(SEED_START, SEED_END + 1);
	String supp    = guide.substring(SUPP_START, SUPP_END + 1);
	String g1Str   = guide.substring(G1, G1 + 1);
	String central = guide.substring(CENTRAL_START, CENTRAL_END + 1);

	// Filename
	String fileName = guide + ".json";

	return interactionsPath
            .resolve(subDirL1)
            .resolve(subDirL2)
            .resolve(seed)
            .resolve(supp)
            .resolve(g1Str)
            .resolve(central)
            .resolve(fileName);
    }

    private static void writeStringToFile( String file, String content ) throws IOException {
        Files.write( Paths.get( file ), content.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING );
    }

    protected static void createAccessibilityFile( String accessibilityFile, String sequence ) {
	//System.out.println( "creating accessibility file..." );
	TiledSecondaryStructure target2D = new TiledSecondaryStructure( sequence );
	try {
	    writeStringToFile( accessibilityFile, target2D.toString() );
	    //System.out.println( "accessibility file created" );
	} catch( IOException e ) {
	    Utils.stop( "IO error writing to accessibility file: " + accessibilityFile, 18 );
	}
    }

    // return true if the results already exist; false if they were not
    //     if results were not computed, create results path and file, as well as the lock file
    //     if false is returned, it tells the calling process that results must be generated
    //         these results will need to be passed to the create method that will unlock and save them
    public boolean exist( String transcript, String guide, String transcriptSequence ) { 
	// Base transcript path (without .scan)
	String pctAccess = getTranscriptFile( transcript );

	// Transcript-level .scan directory (e.g. .../SI/SIRT1-201.scan)
	Path resultsPath = Paths.get( pctAccess + ".scan" );

	// New JSON path based on seed/supp/g1/central/guide.json
	Path resultsFile = this.getInteractionsFilePath( this.interactionsPath, transcript, guide );

	// Lock file lives next to the JSON, named guide.lock
	Path lock = resultsFile.getParent().resolve( guide + ".lock" );

	// Transcript accessibilities file stays at the top of the .scan dir
	Path accessibilityPath = resultsPath.resolve( "accessibilities" );

	try {
	    if( Files.exists( resultsPath ) && Files.isDirectory( resultsPath ) ) {
		// .scan directory exists for this transcript
		if( Files.exists( resultsFile ) || Files.exists( lock ) ) {
		    // results exist or are being processed from another request
		    while( Files.exists( lock ) ) {
			// busy-wait until lock disappears
		    }
		    return true;
		} else {
		    // directory exists but no results for this guide yet
		    try {
			// ensure seed/supp/g1/central directories exist
			Files.createDirectories( resultsFile.getParent() );
			// lock the case
			Files.createFile( lock );

			// ensure accessibility file exists
			if( !Files.exists( accessibilityPath ) ) {
			    createAccessibilityFile( accessibilityPath.toString(), transcriptSequence );
			}
		    } catch( FileAlreadyExistsException e ) {
			// lock acquired just before us → wait until it disappears
			while( Files.exists( lock ) ) {
			    // wait
			}
			return true;
		    }
		}
	    } else {
		// resultsPath (.scan) does not exist: create it, lock file, and accessibility file
		try {
		    Files.createDirectories( resultsPath );
		    Files.createDirectories( resultsFile.getParent() );
		    Files.createFile( lock );

		    // compute secondary structure and write to accessibility file
		    createAccessibilityFile( accessibilityPath.toString(), transcriptSequence );
		} catch( FileAlreadyExistsException e ) {
		    // lock acquired just before us → wait until it disappears
		    while( Files.exists( lock ) ) {
			// wait
		    }
		    return true;
		}
	    }
	} catch( IOException e ) {
	    e.printStackTrace();
	    Utils.stop( "I/O file system problem", 9 );
	}
	return false;
    }
    
    // public boolean exist( String transcript, String guide, String transcriptSequence ) { // transcript is the transcript name, guide is the guide sequence
    // 	//System.out.println( "exists( ... " );
    // 	String pctAccess = this.getTranscriptFile( transcript );
    // 	Path resultsPath = Paths.get( pctAccess + ".scan" ); // path to scan directory, if exists
    // 	Path resultsFile = Paths.get( pctAccess + ".scan" + File.separator + guide ); // path to result file, if exists
    // 	Path lock = Paths.get( pctAccess + ".scan" + File.separator + guide + ".lock" ); // lock file as a semaphore
    // 	String accessibilityFile = pctAccess + ".scan" + File.separator + "accessibilities"; // transcript accessibilities
    // 	try {
    // 	    if( Files.exists( resultsPath ) && Files.isDirectory( resultsPath ) ) { // directory exists for this transcript
    // 		if( Files.exists( resultsFile ) || Files.exists( lock ) ) { // results exist or results are being processed from another request
    // 		    // if the file is locked, wait...
    // 		    while( Files.exists( lock ) ) { }
    // 		    //System.out.println( "exists( ... ) returning true" );
    // 		    return true;
    // 		}
    // 		else { // directory exists but no results for this guide
    // 		    // lock the case
    // 		    try {
    // 			Files.createFile( lock );
    // 			if( !Files.exists( Paths.get( accessibilityFile ) ) ) { // check if secondary structure exists
    // 			    // compute secondary structure using tiled2DStructure and write to accessibility file
    // 			    createAccessibilityFile( accessibilityFile, transcriptSequence );
    // 			}
    // 		    } catch( FileAlreadyExistsException e ) {
    // 			// lock acquired just before us
    // 			// wait until it disappears => data will be available
    // 			while( Files.exists( lock ) ) { }
    // 			return true;
    // 		    }
    // 		}
    // 	    }
    // 	    else { // results path does not exit, create it as well as a lock file, and create accessibility file for this transcript
    // 		try {
    // 		    Files.createDirectories( resultsPath );
    // 		    Files.createFile( lock );
    // 		    // compute secondary structure using tiled2DStructure and write to accessibility file
    // 		    createAccessibilityFile( accessibilityFile, transcriptSequence );
    // 		} catch( FileAlreadyExistsException e ) {
    // 		    // lock acquired just before us
    // 		    // wait until it disappears => data will be available
    // 		    while( Files.exists( lock ) ) { }
    // 		    return true;
    // 		}
    // 	    }
    // 	} catch( IOException e ) {
    // 	    e.printStackTrace();
    // 	    Utils.stop( "I/O file system problem", 9 );
    // 	}
    // 	//System.out.println( "exists( ... ) returning false" );
    // 	return false;
    // }

    // public List<JsonInteraction> getInteractions( String transcript, String guide ) {
    // 	String interactionFile = this.getTranscriptFile( transcript ) + ".scan" + File.separator + guide;
    // 	List<JsonInteraction> results = null;
    // 	//System.out.println( "getting from " + interactionFile );
    // 	try {
    // 	    results = JsonInteractionIO.readJsonFile( Paths.get( interactionFile ) );
    // 	    //System.out.println( "results: " + results );
    // 	} catch( IOException e ) {
    // 	    Utils.stop( "IO error reading file: " + interactionFile, 19 );
    // 	}
    // 	return results;
    // }

    public List<JsonInteraction> getInteractions( String transcript, String guide ) {
	// New: build the full path using the seed/supp/g1/central/guide.json layout
	Path interactionPath = this.getInteractionsFilePath( this.interactionsPath, transcript, guide );

	List<JsonInteraction> results = null;
	// System.out.println("getting from " + interactionPath);

	try {
	    results = JsonInteractionIO.readJsonFile( interactionPath );
	    // System.out.println("results: " + results);
	} catch( IOException e ) {
	    Utils.stop( "IO error reading file: " + interactionPath, 19 );
	}
	return results;
    }
    
    // public void putInteractions( String transcript, String guide, List<JsonInteraction> interactions ) {
    // 	String interactionFile = this.getTranscriptFile( transcript ) + ".scan" + File.separator + guide;
    // 	Path lock = Paths.get( this.getTranscriptFile( transcript ) + ".scan" + File.separator + guide + ".lock" ); // lock file
    // 	try {
    // 	    JsonInteractionIO.writeJsonFile( Paths.get( interactionFile ), interactions );
    // 	} catch( IOException e ) {
    // 	    Utils.stop( "IO error writing file: " + interactionFile, 20 );
    // 	}
    // 	try {
    // 	    Files.delete( lock ); // delete the lock file
    // 	} catch( IOException e ) {
    // 	    Utils.stop( "IO error deleting lock file: " + lock, 21 );
    // 	}   
    // }

    public void putInteractions( String transcript, String guide, List<JsonInteraction> interactions ) {
	// New: use the same path builder (seed/supp/g1/central/guide.json)
	Path interactionPath = this.getInteractionsFilePath( this.interactionsPath, transcript, guide );

	// Lock file lives in the same directory, named by the guide (no .json)
	Path lock = interactionPath.getParent().resolve( guide + ".lock" );

	try {
	    JsonInteractionIO.writeJsonFile( interactionPath, interactions );
	} catch( IOException e ) {
	    Utils.stop( "IO error writing file: " + interactionPath, 20 );
	}

	try {
	    Files.delete( lock ); // delete the lock file
	} catch( IOException e ) {
	    Utils.stop( "IO error deleting lock file: " + lock, 21 );
	}
    }
}
