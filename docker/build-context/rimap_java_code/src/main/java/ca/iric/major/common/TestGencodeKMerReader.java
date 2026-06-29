/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.List;
import java.util.ArrayList;

public class TestGencodeKMerReader {

    public static void main( String[] args ) {
	// KMer kmer = new KMer();
	// // Testing the iterator version
	// GencodeKMerReader GKMreader = new GencodeKMerReader( "/Users/major/Dropbox/Dev/MiRBooking/data/gencode.v34.k31-canonical.fa" );
	// while( GKMreader.hasNext() ) {
	//     kmer = GKMreader.next();
	//     System.out.println( kmer );
	// }
	// Testing the list version
	List<KMer> myList = new ArrayList<>();
	GencodeKMer GKMreader = new GencodeKMer( "/Users/major/Dropbox/Dev/MiRTools/data/gencode.v34.k31-canonical.fa", myList );
	System.out.println( myList.size() );
    }
}
