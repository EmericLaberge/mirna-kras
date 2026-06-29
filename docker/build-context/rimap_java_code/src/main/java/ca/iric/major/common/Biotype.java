/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

// Biotype from Vega
public enum Biotype {
    ProteinCoding( "ProteinCoding" ),
    NonCoding( "NonCoding" ),
    Antisense( "Antisense" ),
    lincRNA( "NonCoding" ),
    miRNA( "miRNA" ),
    reporter( "reporter" ),
    rRNA( "rRNA" ),
    snRNA( "snRNA" ),
    snoRNA( "snoRNA" ),
    tRNA( "tRNA" ),
    Unclassified( "Unclassified" ),
    Unknown( "Unknown" );

    private final String name;

    // Constructor for each enum instance
    Biotype( String name ) { this.name = name; }

    // Getter to retrieve the name representation for each enum instance
    public String getName() { return this.name; } 

    // Optional: a method that returns a pretty string representation of the enum instance
    @Override
    public String toString() { return this.name; }

    // Optional: a method to return an enum instance corresponding to a given string
    // Note: This method throws an IllegalArgumentException if the string does not match any enum instance
    public static Biotype fromString( String name ) {
	for( Biotype biotype : Biotype.values() ) {
	    if( biotype.name.equalsIgnoreCase( name ) ) {
		return biotype;
	    }
	}
	throw new IllegalArgumentException( "No Biotype with name " + name + " found" );
    }
}
