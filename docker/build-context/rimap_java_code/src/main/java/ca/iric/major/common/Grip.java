/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Grip links two common kmers (positioned at position5, 5' of, position3) to a transcript
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 1.0 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

public class Grip {
    private ProteinCodingTranscript pct;
    private int position5p;
    private int position3p;
    private int bridge;

    public Grip( ProteinCodingTranscript pct, int position3p, int position5p, int bridge ) {
	this.pct = pct;
	this.position3p = position3p;
	this.position5p = position5p;
	this.bridge = bridge;
    }

    public int getPosition3p() { return this.position3p; }
    public int getPosition5p() { return this.position5p; }
    public int getBridge()     { return this.bridge; }
    public ProteinCodingTranscript getPCT() { return this.pct; }

    @Override
    public String toString() {
	return "Grip( " + this.pct.getName() + ", " + this.position3p + ", " + this.position5p + ", " + this.bridge + " )";
    }
}
