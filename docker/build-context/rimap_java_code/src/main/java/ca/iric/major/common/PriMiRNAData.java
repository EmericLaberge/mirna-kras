/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

public class PriMiRNAData {
    PriMiRNA priMiRNA;
    SecondaryStructure priMiRNASS;
    public PriMiRNAData( PriMiRNA priMiRNA, SecondaryStructure priMiRNASS ) {
	this.priMiRNA = priMiRNA;
	this.priMiRNASS = priMiRNASS;
	System.out.println( this.priMiRNA.getName() + " done" );
    }
    public PriMiRNA getPriMiRNA() { return this.priMiRNA; }
    public SecondaryStructure getSecondaryStructure() { return this.priMiRNASS; }
}
