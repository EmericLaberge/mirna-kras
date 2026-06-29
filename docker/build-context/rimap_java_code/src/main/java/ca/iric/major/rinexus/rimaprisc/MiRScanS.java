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

public class MiRScanS {

    public static void main( String[] args ) {

	// ***** JSON ObjectMapper
	ObjectMapper mapper = new ObjectMapper();

	// ***** READ ARGUMENTS

	if( args.length != 7 ) {
	    System.out.println( "usage: MiRScanS <guide name> <guide MIMAT> <guide organism> <guide sequence> <target name> <target ENST> > <target sequence>" );
	    Utils.stop( "bye!" , 0 );
	}

	String guideName = args[0];
	String guideMIMAT = args[1];
	String guideOrganism = args[2];
	String guideSequence = args[3];
	String targetName = args[4];
	String targetId = args[5];
	String targetSequence = args[6].replace( 'T', 'U' );

	// ***** PHASE 1: find all seed and supp complementary sites *****
	// seed and supp initiation sites of 4 nucleotides

	MatureMiRNA miRNA = new MatureMiRNA( guideName, guideMIMAT, guideOrganism, guideName );
	miRNA.setSequence( new StringSequence( guideSequence ) );
	String matureSequence = miRNA.getSequence().getSequence();
	String initSeed = matureSequence.substring(  1,  5 ); // seed g2-g5
	String initSupp = matureSequence.substring( 12, 16 ); // supp g13-g16

	// get the protein coding transcript sequence
	ProteinCodingTranscript target = new ProteinCodingTranscript();
	target.setSequence( new StringSequence( targetSequence ) );
	target.setName( targetName );
	target.setGeneId( targetId );

	// get the list of seed and supp initiation sites

	List<Integer> seedInitiationSites = StringSequence.complementSites( targetSequence, initSeed, false );
	List<Integer> suppInitiationSites = StringSequence.complementSites( targetSequence, initSupp, false );
	Set<Integer> seedInitiationSitesCovered = new HashSet<>();

	for( Integer seedSite : seedInitiationSites )
	    for( Integer suppSite : suppInitiationSites ) {
		int g2 = seedSite + 3;
		int g9 = seedSite - 4;
		int g12 = suppSite + 4;
		int bridgeLength = g9 - g12;
		int targetLen = 31;
		if( target.inTranscript( g2 - 29, g2 + 1 ) ) { // if the target at g2 allows to define a valid 31mer...
		    if( bridgeLength > 2 && bridgeLength <= 15 ) { // if the bridge is at most 15 nts...
				
			Guide guide = new Guide( miRNA.getFullName(),
						 miRNA.getMIMATCode(),
						 target,
						 target.getGeneId(),
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

		    Guide guide = new Guide( miRNA.getFullName(),
					     miRNA.getMIMATCode(),
					     target,
					     target.getGeneId(),
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
	for( Integer suppSite : suppInitiationSites ) {
	    Integer seedSite = suppSite + 11;
	    int g2 = seedSite + 3;
	    int g9 = seedSite - 4;
	    int g12 = suppSite + 4; // assume just the A-box
	    int bridgeLength = g9 - g12;
	    int targetLen = 31;
	    if( target.inTranscript( g2 - 29, g2 + 1 ) ) {
		if( bridgeLength > 2 && bridgeLength <= 15 ) {

		    Guide guide = new Guide( miRNA.getFullName(),
					     miRNA.getMIMATCode(),
					     target,
					     target.getGeneId(),
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
		    if( !seedInitiationSitesCovered.contains( seedSite ) && guide.initiate() ) {
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
    }
}

