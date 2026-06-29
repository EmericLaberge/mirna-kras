/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.List;
import java.util.ArrayList;

public class TestRNASeqReads {

    public static void main( String[] args ) {
	List<Read> myList = new ArrayList<>();
	RNASeqReads reader = new RNASeqReads( "/Users/major/Dropbox/Dev/MiRTools/data/SRR18046904.fastq", myList );
	System.out.println( myList.size() );
    }
}
