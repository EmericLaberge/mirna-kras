/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

public class MatureMiRNA {

    public static final int MIRNAMINSIZE = 18; // 20 as in MiRGeneDB + possible cut of 2 nt in 3' (enzymes)
    public static final int MIRNAMAXSIZE = 26;

    String fullName = "";
    String MIMATCode = "";
    String organism = "";
    String name = "";
    Sequence sequence; // Sequence holder
    String matureSequence; // String sequence holder
    String seed;
    String offsetSeed;
    String distantSeed;
    String fullSeed;
    String supp;
    String offsetSupp;
    String distantSupp;

    // getters
    public String getFullName() { return this.fullName; }
    public String getMIMATCode() { return this.MIMATCode; }
    public String getOrganism() { return this.organism; }
    public String getName() { return this.name; }
    public Sequence getSequence() { return this.sequence; }
    public String getStringSequence() { return this.sequence.getSequence(); }
    public String getMatureSequence() { return this.matureSequence; }
    public String getSeed() { return this.seed; }
    public String getFullSeed() { return this.fullSeed; }
    public String getOffsetSeed() { return this.offsetSeed; }
    public String getDistantSeed() { return this.distantSeed; }
    public String getSupp() { return this.supp; }
    public String getOffsetSupp() { return this.offsetSupp; }
    public String getDistantSupp() { return this.distantSupp; }

    // setters
    public void setFullName( String fullName ) { this.fullName = fullName; }
    public void setMIMATCode( String MIMATCode ) { this.MIMATCode = MIMATCode; }
    public void setOrganism( String organism ) { this.organism = organism; }
    public void setName( String name ) { this.name = name; }
    public void setSequence( Sequence sequence ) {
	if( sequence.length() < MIRNAMINSIZE || sequence.length() > MIRNAMAXSIZE ) Utils.stop( "MatureMiRNA, microRNA sequence of wrong size: " + sequence.length(), 0 );
	this.sequence = sequence;
    }

    // constructors
    public MatureMiRNA( Sequence sequence ) {
	if( sequence.length() < MIRNAMINSIZE || sequence.length() > MIRNAMAXSIZE ) Utils.stop( "MatureMiRNA, microRNA sequence of wrong size: " + sequence.length(), 0 );
	this.sequence = sequence;
	this.matureSequence = this.sequence.getSequence(); /// mature full sequence
	this.seed = this.matureSequence.substring(  Guide.G2,  Guide.G5 + 1 ); // seed g2-g5
	this.offsetSeed = this.matureSequence.substring( Guide.G3, Guide.G6 + 1 ); // offset seed g3-g6
	this.distantSeed = this.matureSequence.substring( Guide.G4, Guide.G7 + 1 ); // distant seed g4-g8
	this.fullSeed = this.matureSequence.substring( Guide.G2, Guide.G8 + 1 ); // full seed g2-g8
	//this.supp = this.matureSequence.substring( 11, 14 ); // supp g12-g15
	this.supp = this.matureSequence.substring( Guide.G13, Guide.G16 + 1 ); // supp g13-g16
	this.offsetSupp = this.matureSequence.substring( Guide.G14, Guide.G17 + 1 ); // offset supp g14-g17
    }
    
    public MatureMiRNA( String fullName, String MIMATCode, String organism, String name ) {
	this.fullName = fullName;
	this.MIMATCode = MIMATCode;
	this.organism = organism;
	this.name = name;
    }

    //Override
    public String toString() {
	String out = "";
	out += ">" + this.getFullName() + " " + this.getMIMATCode() + " " + this.getOrganism() + " " + this.getName() + "\n";
	out += this.getSequence();
	return out;
    }
}
