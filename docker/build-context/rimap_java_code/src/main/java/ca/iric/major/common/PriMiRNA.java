/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

public class PriMiRNA {

    String organism = "";
    String name = "";
    Sequence sequence; // Sequence holder

    // getters
    public String getOrganism() { return this.organism; }
    public String getName() { return this.name; }
    public Sequence getSequence() { return this.sequence; }
    public String getStringSequence() { return this.sequence.getSequence(); }

    // setters
    public void setOrganism( String organism ) { this.organism = organism; }
    public void setName( String name ) { this.name = name; }
    public void setSequence( Sequence sequence ) { this.sequence = sequence; }

    // constructors
    public PriMiRNA() {}
    
    public PriMiRNA( String organism, String name ) {
	this.organism = organism;
	this.name = name;
    }

    //Override
    public String toString() {
	String out = "";
	out += ">" + this.getName() + "\n";
	out += this.getSequence();
	return out;
    }
}
