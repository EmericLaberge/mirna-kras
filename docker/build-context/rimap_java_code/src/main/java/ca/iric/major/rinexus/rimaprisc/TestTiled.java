/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import ca.iric.major.common.*;

public class TestTiled {

    public static void main( String[] args ) {
	int[] tileSizes = { 18 };
	double[] tileT = { 10 };
	TiledSecondaryStructure tss = new TiledSecondaryStructure( "UGAUGCUACGAGCUGAUUGAUGCUGAGCUUGACGUAUGCUGAUGCACUAUGCUGACUG", tileSizes, tileT, 0.0 );
	System.out.println( tss.toCSV() );
    }
}
