/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import ca.iric.major.common.*;

import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;

import java.io.FileNotFoundException;
import java.io.IOException;

import com.fasterxml.jackson.databind.ObjectMapper;

public class MiRScanK {

    public static void main( String[] args ) {

	// ***** JSON ObjectMapper
	ObjectMapper mapper = new ObjectMapper();

	// ***** PRESENTATION
	//	System.out.println( "MiRScanK v.1.0 June 20, 2024 - Major Lab, IRIC, Université de Montréal" );

	// ***** READ ARGUMENTS

	if( args.length != 2 ) {
	    System.out.println( "usage: MiRScan <31mer sequence> <MicroRNA name>" );
	    System.out.println( "   31mer sequence: RNA or DNA 31mer sequence" );
	    System.out.println( "   MicroRNA name: miRBase name" ); 
	    Utils.stop( "bye!" , 0 );
	}

	// ***** READ system variable for the paths
	String myVarValue = System.getenv( "MiRBooking_DATAPATH" );
	if( myVarValue != null ) {
	    //System.out.println("Using path: " + myVarValue );
	} else {
	    myVarValue = "/Users/major/Dropbox/Dev/MiRBooking/data/";
	}

	char fileSystemSeparator = '/';
	// ***** READ OS and use appropriate file system separator
	String os = System.getProperty( "os.name" ).toLowerCase();
        if (os.contains( "win" ) ) {
            //System.out.println( "running under Windows" );
	    fileSystemSeparator = '\\';
        } else if( os.contains( "mac" ) || os.contains( "nix" ) || os.contains( "nux" ) || os.contains( "aix" ) ) {
            //System.out.println( "running under MacOS/Unix/Linux" );
        } else {
            System.out.println( "Unrecognized OS: " + os);
	    Utils.stop( "Unrecognized OS", 0 );
	}

	String reporterPath = myVarValue + "Reporters";
	String miRBasePath  = myVarValue + "Homo_sapiens_mature";

	String kmer = args[0].toUpperCase().replace( 'T', 'U' ); // make it an uppercase RNA sequence
	if( kmer.length() != 31 || !StringSequence.isRNA( kmer ) ) {
            System.out.println( "Illegal 31mer: " + kmer );
	    Utils.stop( "illegal 31mer!", 0 );
	}
	    
	String miRName = args[1];

	// ***** READ miR

	String miRSubdir = miRName.split( "-" )[0].toLowerCase();
	List<MatureMiRNA> matures = new ArrayList<>();
	try {
	    MiRBaseMature miRNAReader = new MiRBaseMature( miRBasePath + fileSystemSeparator + miRSubdir + fileSystemSeparator + miRName, matures );
	} catch( FileNotFoundException e ) {
	    Utils.stop( "microRNA " + miRBasePath + fileSystemSeparator + miRSubdir + fileSystemSeparator + miRName + " not found", 0 );
	} catch( IOException e ) {
	    Utils.stop( "problem reading " + miRBasePath + fileSystemSeparator + miRSubdir + fileSystemSeparator + miRName, 0 );
	}

	MatureMiRNA miRNA = matures.get( 0 );
	String matureSequence = miRNA.getSequence().getSequence();
	String initSeed = matureSequence.substring(  1,  5 ); // seed g2-g5
	String initSupp = matureSequence.substring( 12, 16 ); // supp g13-g16

	// define target sequence
	String targetSequence = kmer;

		// get the list of seed and supp initiation sites

	List<Integer> seedInitiationSites = StringSequence.complementSites( targetSequence, initSeed, false );
	List<Integer> suppInitiationSites = StringSequence.complementSites( targetSequence, initSupp, false );
	Set<Integer> seedInitiationSitesCovered = new HashSet<>();
	Set<Integer> suppInitiationSitesCovered = new HashSet<>();

	for( Integer seedSite : seedInitiationSites )
	    for( Integer suppSite : suppInitiationSites ) {
		int g2 = seedSite + 3;
		int g9 = seedSite - 4;
		int g12 = suppSite + 4;
		int bridgeLength = g9 - g12;
		int targetLen = 31;
		if( bridgeLength > 2 && bridgeLength <= 15 ) { // if the bridge is at most 15 nts...
		    // synchronized initiation sites, make sure it fits within the target
		    //System.out.println( seedSite + " " + targetSequence.substring( seedSite, seedSite+4 ) + "... bridge(" + bridgeLength + ") ..." + targetSequence.substring( suppSite, suppSite+4 ) + " " + suppSite );
				
		    // locally fold the target by taking 3 x 31mer centered at g2 - 29, g2 + 1...
		    // TiledSecondaryStructure target2D = new TiledSecondaryStructure( targetSequence, g2 - 29 - 31, g2 + 1 + 31  );
		    Guide guide = new Guide( miRNA.getFullName(),
					     miRNA.getMIMATCode(),
					     targetSequence,
					     matureSequence,
					     g2,
					     g12,
					     StringSequence.reverseComplement( initSeed ),
					     StringSequence.reverseComplement( initSupp ),
					     bridgeLength,
					     g2 - targetLen + 2, // align with 31mer + bridgeLength
					     g2 + 1,
					     true ); // align with 31mer
		    guide.fold();
		    if( !seedInitiationSitesCovered.contains( seedSite ) && guide.bind() ) {
			System.out.println( guide ); // good seed or can at least initiate duplex formation + supp at least 4mer
			// add the detailed accessibility of the target
			//for( int nt = 31; nt <= 61; nt++ ) System.out.println( (nt-30) + "(" + target2D.getNucleotide( nt ) + "): " + target2D.getReactivity( nt ) );
			seedInitiationSitesCovered.add( seedSite );
			suppInitiationSitesCovered.add( suppSite );
			// output JSON
			try {
			    String jsonString = mapper.writeValueAsString( guide );
			    System.out.println( jsonString );
			} catch( Exception e ) {
			    e.printStackTrace();
			}
		    }
		}
	    }
	// check the sites with possible seeds only
	//System.out.println( "seed sites that bind" );
	for( Integer seedSite : seedInitiationSites ) {
	    int suppSite = seedSite - 11;
	    int g2 = seedSite + 3;
	    int g9 = seedSite - 4;
	    int g12 = suppSite + 4; // assume just the A-box                                                                                                                                                                                                         
	    int bridgeLength = g9 - g12;
	    int targetLen = 31;
	    if( !seedInitiationSitesCovered.contains( seedSite ) ) {
		if( bridgeLength > 2 && bridgeLength <= 15 ) {
		    // synchronized initiation sites, make sure it fits within the target                                                                                                                                                                           
		    //System.out.println( seedSite + " " + targetSequence.substring( seedSite, seedSite+4 ) + "... bridge(" + bridgeLength + ") ..." + targetSequence.substring( suppSite, suppSite+4 ) + " " + suppSite );                                         
		    //System.out.println( seedSite + " " + targetSequence.substring( seedSite, seedSite+4 ) + "... bridge(" + bridgeLength + ") ..." + targetSequence.substring( suppSite, suppSite+4 ) + " " + suppSite );
		    
		    // locally fold the target by taking 3 x 31mer centered at g2 - 29, g2 + 1...
		    //TiledSecondaryStructure target2D = new TiledSecondaryStructure( targetSequence, g2 - 29 - 31, g2 + 1 + 31  );
		    Guide guide = new Guide( miRNA.getFullName(),
					     miRNA.getMIMATCode(),
					     targetSequence,
					     matureSequence,
					     g2,
					     g12,
					     StringSequence.reverseComplement( initSeed ),
					     StringSequence.reverseComplement( initSupp ),
					     bridgeLength,
					     g2 - targetLen + 2, // align with 31mer + bridgeLenght
					     g2 + 1,
					     false ); // align with 31mer
		    guide.fold();
		    if( !seedInitiationSitesCovered.contains( seedSite ) && guide.bind() ) {
			System.out.println( guide ); // good seed or can at least initiate duplex formation + supp at least 4mer
			// add the detailed accessibility of the target
			//for( int nt = 31; nt <= 61; nt++ ) System.out.println( (nt-30) + "(" + target2D.getNucleotide( nt ) + "): " + target2D.getReactivity( nt ) );
			seedInitiationSitesCovered.add( seedSite );
			
			// output JSON
			try {
			    String jsonString = mapper.writeValueAsString( guide );
			    System.out.println( jsonString );
			} catch( Exception e ) {
			    e.printStackTrace();
			}
		    }
		}
	    }
	}
	// check the sites with possible supp only
	//System.out.println( "supp sites that initiate" );
	for( Integer suppSite : suppInitiationSites ) {
	    Integer seedSite = suppSite + 11;
	    int g2 = seedSite + 3;
	    int g9 = seedSite - 4;
	    int g12 = suppSite + 4; // assume just the A-box
	    int bridgeLength = g9 - g12;
	    int targetLen = 31;
	    if( bridgeLength > 2 && bridgeLength <= 15 ) {
		// synchronized initiation sites, make sure it fits within the target
		//System.out.println( seedSite + " " + targetSequence.substring( seedSite, seedSite+4 ) + "... bridge(" + bridgeLength + ") ..." + targetSequence.substring( suppSite, suppSite+4 ) + " " + suppSite );

		// locally fold the target by taking 3 x 31mer centered at g2 - 29, g2 + 1...
		//TiledSecondaryStructure target2D = new TiledSecondaryStructure( targetSequence, g2 - 29 - 31, g2 + 1 + 31  );
		Guide guide = new Guide( miRNA.getFullName(),
					 miRNA.getMIMATCode(),
					 targetSequence,
					 matureSequence,
					 g2,
					 g12,
					 StringSequence.reverseComplement( initSeed ),
					 StringSequence.reverseComplement( initSupp ),
					 bridgeLength,
					 g2 - targetLen + 2, // align with 31mer + bridgeLength
					 g2 + 1,
					 true ); // align with 31 mer
		guide.fold();
		if( !suppInitiationSitesCovered.contains( suppSite ) && guide.initiate() ) {
		    System.out.println( guide ); // print if at least seed can initiate duplex formation
		    // add the detailed accessibility of the target
		    //for( int nt = 31; nt <= 61; nt++ ) System.out.println( (nt-30) + "(" + target2D.getNucleotide( nt ) + "): " + target2D.getReactivity( nt ) );
		    suppInitiationSitesCovered.add( suppSite );
		    // output JSON
		    try {
			String jsonString = mapper.writeValueAsString( guide );
			System.out.println( jsonString );
		    } catch( Exception e ) {
			e.printStackTrace();
		    }
		}	
	    }
	}
    }
}
