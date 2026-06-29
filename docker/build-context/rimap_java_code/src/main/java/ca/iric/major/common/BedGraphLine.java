/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

public class BedGraphLine implements Comparable<BedGraphLine>{
    private String chromosome;
    private int start;
    private int end;
    private double conservation;

    // getters
    public String getChromosome() { return this.chromosome; }
    public int getStart() { return this.start; }
    public int getEnd() { return this.end; }
    public double getConservation() { return this.conservation; }

    public BedGraphLine( String chromosome, int start, int end, double conservation ) {
	this.chromosome = chromosome;
	this.start = start;
	this.end = end;
	this.conservation = conservation;
    }

    @Override
    public int compareTo( BedGraphLine other ) {
        return Integer.compare( this.start, other.start );
    }

    @Override
    public String toString() {
	return this.chromosome + "\t" + this.start + "\t" + this.end + "\t" + this.conservation;
    }
}
