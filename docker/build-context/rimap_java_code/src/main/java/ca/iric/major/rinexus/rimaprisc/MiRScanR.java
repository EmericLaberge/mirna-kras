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

public class MiRScanR {

    public static void main( String[] args ) {

	// ***** PRESENTATION
	System.out.println( "MiRScanR v.1.0 January 22, 2024 - Major Lab, IRIC, Université de Montréal" );

	// ***** READ ARGUMENTS

	System.out.println( "number of args: " + args.length );
	if( args.length != 2 ) {
	    System.out.println( "usage: MiRScanR <ReporterTranscript name> <MicroRNA name>" );
	    System.out.println( "   Transcript: Reporter Transcript name" );
	    System.out.println( "   MicroRNA name: miRBase name" ); 
	    Utils.stop( "bye!" , 0 );
	}

	String reporterPath = "/Users/major/Dropbox/Dev/MiRBooking/data/Reporters";
	String miRBasePath  = "/Users/major/Dropbox/Dev/MiRBooking/data/Homo_sapiens_mature";

	String repName = args[0];
	String miRName = args[1];

	// ***** READ reporter and miR

	List<ReporterTranscript> targets = new ArrayList<>();
	Reporter ReporterReader = new Reporter( reporterPath + "/" + repName, targets );
	
	String miRSubdir = miRName.split( "-" )[0].toLowerCase();
	List<MatureMiRNA> matures = new ArrayList<>();
	try {
	    MiRBaseMature miRNAReader = new MiRBaseMature( miRBasePath + "/" + miRSubdir + "/" + miRName, matures );
	} catch( FileNotFoundException e ) {
	    Utils.stop( "microRNA " + miRBasePath + "/" + miRSubdir + "/" + miRName + " not found", 0 );
	} catch( IOException e ) {
	    Utils.stop( "problem reading " + miRBasePath + "/" + miRSubdir + "/" + miRName, 0 );
	}

	// search all possible target sites for each mature in each target

	System.out.println( "MiRScanR will be applied to:" );
	for( int i = 0; i < matures.size(); i++ )
	    for( int j = 0; j < targets.size(); j++ )
		System.out.println( matures.get( i ).getName() + ":" + targets.get( j ).getName() );

	for( int i = 0; i < matures.size(); i++ )
	    for( int j = 0; j < targets.size(); j++ ) {
	
		// ***** PHASE 1: find all seed and supp complementary sites *****

		// seed and supp initiation sites of 4 nucleotides

		MatureMiRNA miRNA = matures.get( i );
		String matureSequence = miRNA.getSequence().getSequence();
		String initSeed = matureSequence.substring(  1,  5 ); // seed g2-g5
		String initSupp = matureSequence.substring( 12, 16 ); // supp g13-g16

		// get the protein coding transcript sequence
		ReporterTranscript target = targets.get( i );
		String targetSequence = target.getSequence().getSequence();

		// tile-fold the transcript to determine MRE accessibilities
		//TiledSecondaryStructure target2D = new TiledSecondaryStructure( targetSequence );

		// get the list of seed and supp initiation sites

		List<Integer> seedInitiationSites = StringSequence.complementSites( targetSequence, initSeed, false );
		List<Integer> suppInitiationSites = StringSequence.complementSites( targetSequence, initSupp, false );
		Set<Integer> seedInitiationSitesCovered = new HashSet<>();

		// debug initiation sites

		//System.out.println( "initSeed: " + initSeed );
		//for( Integer site : seedInitiationSites )
		//    System.out.println( site + " : " + targetSequence.substring( site, site+4 ) + " (" + target.stringRegion( site, site+3 ) + ")" );

		//System.out.println( "initSupp: " + initSupp );
		//for( Integer site : suppInitiationSites )
		//    System.out.println( site + " : " + targetSequence.substring( site, site+4 ) + " (" + target.stringRegion( site, site+3 ) + ")" );

		// check how many pairs satisfy the bridge condition
		for( Integer seedSite : seedInitiationSites )
		    for( Integer suppSite : suppInitiationSites ) {
			int g2 = seedSite + 3;
			int g9 = seedSite - 4;
			int g12 = suppSite + 4;
			int bridgeLength = g9 - g12;
			int targetLen = 31;
			if( target.inTranscript( g2 - 29, g2 + 1 ) ) {
			    if( bridgeLength > 2 && bridgeLength <= 15 ) {
				// synchronized initiation sites, make sure it fits within the target
				// locally fold the target by taking 3 x 31mer centered at g2 - 29, g2 + 1...
				//TiledSecondaryStructure target2D = new TiledSecondaryStructure( targetSequence, g2 - 29 - 31, g2 + 1 + 31  );
				Guide guide = new Guide( miRNA.getFullName(),
							 miRNA.getMIMATCode(),
							 target,
							 "reporter",
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
				    //for( int nt = 31; nt <= 61; nt++ ) System.out.println( nt + "(" + target2D.getNucleotide( nt ) + "): " + target2D.getReactivity( nt ) );
				    seedInitiationSitesCovered.add( seedSite );
				}
			    }
			}
		    }
		// check the sites with possible seeds only
		for( Integer seedSite : seedInitiationSites ) {
		    int suppSite = seedSite - 11;
		    int g2 = seedSite + 3;
		    int g9 = seedSite - 4;
		    int g12 = suppSite + 4; // assume just the A-box
		    int bridgeLength = g9 - g12;
		    int targetLen = 31;
		    if( !seedInitiationSitesCovered.contains( seedSite ) && target.inTranscript( g2 - 29, g2 + 1 ) ) {
			if( bridgeLength > 2 && bridgeLength <= 15 ) {
			    // synchronized initiation sites, make sure it fits within the target
			    // locally fold the target by taking 3 x 31mer centered at g2 - 29, g2 + 1...
			    //TiledSecondaryStructure target2D = new TiledSecondaryStructure( targetSequence, g2 - 29 - 31, g2 + 1 + 31  );
			    Guide guide = new Guide( miRNA.getFullName(),
						     miRNA.getMIMATCode(),
						     target,
						     "reporter",
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
				//for( int nt = 31; nt <= 61; nt++ ) System.out.println( nt + "(" + target2D.getNucleotide( nt ) + "): " + target2D.getReactivity( nt ) );
				seedInitiationSitesCovered.add( seedSite );
			    }
			}
		    }
		}
		// check the sites with possible supp only
		for( Integer suppSite : suppInitiationSites ) {
		    Integer seedSite = suppSite + 11;
		    int g2 = seedSite + 3;
		    int g9 = seedSite - 4;
		    int g12 = suppSite + 4; // assume just the A-box
		    int bridgeLength = g9 - g12;
		    int targetLen = 31;
		    if( target.inTranscript( g2 - 29, g2 + 1 ) ) {
			if( bridgeLength > 2 && bridgeLength <= 15 ) {
			    // synchronized initiation sites, make sure it fits within the target
			    // locally fold the target by taking 3 x 31mer centered at g2 - 29, g2 + 1...
			    //TTiledSecondaryStructure target2D = new TiledSecondaryStructure( targetSequence, g2 - 29 - 31, g2 + 1 + 31  );
			    Guide guide = new Guide( miRNA.getFullName(),
						     miRNA.getMIMATCode(),
						     target,
						     "reporter",
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
			    if( !seedInitiationSitesCovered.contains( seedSite ) && guide.bind() ) {
				System.out.println( guide ); // print if at least seed can initiate duplex formation
				// add the detailed accessibility of the target
				//for( int nt = 31; nt <= 61; nt++ ) System.out.println( nt + "(" + target2D.getNucleotide( nt ) + "): " + target2D.getReactivity( nt ) );
				seedInitiationSitesCovered.add( seedSite );
			    }	
			}
		    }
		}
	    }
    }
}
