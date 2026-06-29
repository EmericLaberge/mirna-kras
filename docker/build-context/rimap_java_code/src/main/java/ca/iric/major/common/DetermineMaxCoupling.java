/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.HashSet;
import java.util.Arrays;
import java.util.Collections;
//import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.IntStream;
import java.util.stream.Collectors;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.FileOutputStream;

public class DetermineMaxCoupling {

    public static List<ProteinCodingTranscript> filterLongestVariants( List<ProteinCodingTranscript> transcripts ) {
        // Group by family and keep the longest variant in each group
        Map<String,ProteinCodingTranscript> longestVariants = transcripts.stream()
            .collect(Collectors.toMap( ProteinCodingTranscript::getFamily, t -> t, (t1, t2) -> t1.getIntLength() >= t2.getIntLength() ? t1 : t2 ));

        // Return the list of longest variants
        return new ArrayList<>(longestVariants.values());
    }

    public static void main( String[] args ) {

	// ***** PRESENTATION
	System.out.println( "DetermineMaxCoupling v.1.0 November 11, 2024 - Major Lab, IRIC, Université de Montréal" );

	// ***** READ ARGUMENTS

	if( args.length != 3 ) {
	    System.out.println( "usage: DetermineMaxCoupling <Gencode version> <seed grip> <supp grip>" );
	    System.out.println( "   DetermineMaxCoupling v1.0 takes gencode version, seed and supp grip lengths as arguments, ex) v46 4 3" );
	    Utils.stop( "bye!" , 0 );
	}

	// read the environment variable and set gencodePath
	String dataPath = System.getenv( "DATA_PATH" );
	String gencodeFastaFile = dataPath + "gencode." + args[0] + ".pc_transcripts.fa";
	System.out.println( "Processing PCTs in " + gencodeFastaFile );

	int seed = Integer.parseInt( args[1] );
	int supp = Integer.parseInt( args[2] );

	// coupling

	GencodePCTranscript genPCT = null; // define a variable to process transcripts from Fasta file
	// read the PCT from gencode
	try {
	    genPCT = new GencodePCTranscript( gencodeFastaFile); // read all PCTs in gencode
	} catch( FileNotFoundException e ) {
	    e.printStackTrace();
	} catch( IOException e ) {
	    e.printStackTrace();
	}

	List<String> exclusions = Arrays.asList( "AAAA", "UUUU", "GGGG", "CCCC", "AUAUAUAU" );
	List<ProteinCodingTranscript> empty = new ArrayList<>();
	List<ProteinCodingTranscript> filteredTranscripts = filterLongestVariants( genPCT.getList() );
	int numberOfTranscripts = filteredTranscripts.size();

	System.out.println( "number of transcripts: " + numberOfTranscripts );
	System.out.println( "expected number of pairs will be below " + (long)(numberOfTranscripts * numberOfTranscripts) );

        // Now proceed with the n^2 loop using filteredTranscripts
        //ProteinCodingTranscript tX;
        //ProteinCodingTranscript tY;
        //String geneX, geneY = null;

	System.out.println( "coupling( " + seed + ", " + supp + " )..." );
		
	// int maxCoupling = Integer.MIN_VALUE;
	// int minCoupling = Integer.MAX_VALUE;
	// long sumCoupling = 0;
	// long numberOfPairs = 0;

	IntStream.range(0, numberOfTranscripts - 1).parallel().forEach(i -> {
            ProteinCodingTranscript tX = filteredTranscripts.get(i);
            String geneX = tX.getFamily();

            // Use a nested loop for the inner iteration
            for (int j = i + 1; j < numberOfTranscripts; j++) {
                ProteinCodingTranscript tY = filteredTranscripts.get(j);
                String geneY = tY.getFamily();

                // Compute coupling stats
                List<ProteinCodingTranscript> required = new ArrayList<>();
                required.add(tX);
                required.add(tY);

                GripMap gripMap = new GripMap(required, empty, empty, 3, 15, exclusions, seed, supp);
                GuideMap guideMap = new GuideMap(gripMap);
                int coupling = guideMap.getGrips().size();

                // Print results (order is not important)
                System.out.println(geneX + "\t" + geneY + "\t" + coupling);
            }
        });

	// for (int i = 0; i < numberOfTranscripts - 1; i++) {
	//     tX = filteredTranscripts.get(i);
	//     geneX = tX.getFamily();
	//     for (int j = i + 1; j < numberOfTranscripts; j++) {
	// 	tY = filteredTranscripts.get(j);
	// 	geneY = tY.getFamily();
	// 	// Compute coupling stats
	// 	List<ProteinCodingTranscript> required = new ArrayList<>();
	// 	required.add(tX);
	// 	required.add(tY);
		    
	// 	// int tmpMinCoupling = minCoupling;
	// 	// int tmpMaxCoupling = maxCoupling;
			
	// 	GripMap gripMap = new GripMap( required, empty, empty, 3, 15, exclusions, seed, supp );
	// 	GuideMap guideMap = new GuideMap(gripMap);
	// 	int coupling = guideMap.getGrips().size();
	// 	System.out.println( geneX + "\t" + geneY + "\t" + coupling );
	// 	//		System.out.println( geneX + "\t" + geneY );
  
	// 	// if (coupling > maxCoupling) maxCoupling = coupling;
	// 	// if (coupling < minCoupling) minCoupling = coupling;
	// 	// sumCoupling += coupling;
	// 	// numberOfPairs++;
		    
	// 	// if(maxCoupling > tmpMaxCoupling ) 
	// 	//     System.out.println( "maxCoupling: " + maxCoupling );
	// 	// if( minCoupling < tmpMinCoupling )
	// 	//     System.out.println( "minCoupling: " + minCoupling );
	// 	// if( maxCoupling > tmpMaxCoupling || minCoupling < tmpMinCoupling )
	// 	//     System.out.println( "number of pairs so far: " + numberOfPairs );
	//     }
	// }

	// output final values

	// System.out.println( "coupling( " + seed + ", " + supp + " ):" );
	// System.out.println( "max coupling: " + maxCoupling );
	// System.out.println( "min coupling: " + minCoupling );
	// System.out.println( "avg coupling: " + sumCoupling / numberOfPairs );

	// System.out.println( "Bye!" );
    }
}
