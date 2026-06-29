/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Tuple for guides (composed of a sequence and the kmers used to build it)
 *
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2021 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.lang.IllegalArgumentException;
import java.util.Set;
import java.util.HashSet;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.fasterxml.jackson.annotation.JsonProperty;

public class Guide implements Comparable<Guide> {

    public static final int G1 = 0; // g1 index
    public static final int G2 = 1; // g2 index
    public static final int G3 = 2; // g3 index
    public static final int G4 = 3; // g4 index
    public static final int G5 = 4; // g5 index
    public static final int G6 = 5; // g6 index
    public static final int G7 = 6; // g7 index
    public static final int G8 = 7; // g8 index
    public static final int G12 = 11; // g12 index
    public static final int G13 = 12; // g13 index
    public static final int G14 = 13; // g14 index
    public static final int G15 = 14; // g15 index
    public static final int G16 = 15; // g16 index
    public static final int G17 = 16; // g17 index

    public static final double SEEDLESSKD = 1000000000;

    public static final Set<String> surrounded = Set.of( "((", "({", "{(", " (", "( " );
 
    private class Triplet {
	String s1, s2, s3;
	public String getS1() { return this.s1; }
	public String getS2() { return this.s2; }
	public String getS3() { return this.s3; }
	public Triplet( String s1, String s2, String s3 ) {
	    this.s1 = s1;
	    this.s2 = s2;
	    this.s3 = s3;
	}
    }

    public enum Seed {

	// A1-7mer (pour Bartel)
	A1SVmer( "7mer" ),           // g2-g8 + A1
	// 7mer
	SVmer( "7mer" ),                // g2-g8
	// 6mer
	A1SXmer( "6mer" ),           // g2-g7 + A1
	SXmer( "6mer" ),                // g2-g7
	OffsetSXmer( "offset-6mer" ),   // g3-g8
	// 5mer
	A1FVmer( "5mer" ),           // g2-g6 + A1
	FVmer( "5mer" ),                // g2-g6
	OffsetFVmer( "offset-5mer" ),   // g3-g7
	DistantFVmer( "distant-5mer" ), // g4-g8
	// 4mer
	A1FRmer( "4mer" ),           // g2-g5 + A1
	FRmer( "4mer" ),                // g2-g5
	OffsetFRmer( "offset-4mer" ),   // g3-g6
	DistantFRmer( "distant-4mer" ), // g4-g7
	// seedless
	SEEDLESS( "seedless" );

	private final String pprint; // a pretty print string for the type
	Seed( String pp ) { this.pprint = pp; }
	public String interactionName() { return this.pprint; }
	public double getKd() {
	    switch( this ) {
	    case A1SVmer: return 4;
	    case SVmer: return 12;
	    case A1SXmer: return 125;
	    case SXmer: return 343;
	    case OffsetSXmer: return 800;
	    case A1FVmer: return 490;
	    case FVmer: return 1504;
	    case OffsetFVmer: return 26400;
	    case DistantFVmer: return 295837;
	    case A1FRmer: return 4000;
	    case FRmer: return 12276;
	    case OffsetFRmer: return 115750;
	    case DistantFRmer: return 8455996;
	    default: return SEEDLESSKD; // seedless Kd
	    }
	}
	public int getGFirst() {
	    switch( this ) {
	    case A1SVmer: return Guide.G2;
	    case SVmer: return Guide.G2;
	    case A1SXmer: return Guide.G2;
	    case SXmer: return Guide.G2;
	    case OffsetSXmer: return Guide.G3;
	    case A1FVmer: return Guide.G2;
	    case FVmer: return Guide.G2;
	    case OffsetFVmer: return Guide.G3;
	    case DistantFVmer: return Guide.G4;
	    case A1FRmer: return Guide.G2;
	    case FRmer: return Guide.G2;
	    case OffsetFRmer: return Guide.G3;
	    case DistantFRmer: return Guide.G4;
	    }
	    return 0;
	}
	public int getGLast() {
	    switch( this ) {
	    case A1SVmer: return Guide.G8;
	    case SVmer: return Guide.G8;
	    case A1SXmer: return Guide.G7;
	    case SXmer: return Guide.G7;
	    case OffsetSXmer: return Guide.G8;
	    case A1FVmer: return Guide.G6;
	    case FVmer: return Guide.G6;
	    case OffsetFVmer: return Guide.G7;
	    case DistantFVmer: return Guide.G8;
	    case A1FRmer: return Guide.G4;
	    case FRmer: return Guide.G4;
	    case OffsetFRmer: return Guide.G6;
	    case DistantFRmer: return Guide.G7;
	    }
	    return 0;
	}
    }
	
    public enum Supplementary {
	// 5mer
	FVmer( "5mer" ),                    // g13-g17
	// 4mer
	FRmer( "4mer" ),                    // g13-g16
	OffsetFRmer( "offset-4mer" ),       // g14-g17
	SUPPLEMENTLESS( "suppless" );

	private final String pprint; // a pretty print string for the type
	Supplementary( String pp ) { this.pprint = pp; }
	public String suppTypeName() { return this.pprint; }
	public int getGFirst() {
	    switch( this ) {
	    case FVmer: return Guide.G13;
	    case FRmer: return Guide.G13;
	    case OffsetFRmer: return Guide.G14;
	    }
	    return 0;
	}
	public int getGLast() {
	    switch( this ) {
	    case FVmer: return G17;
	    case FRmer: return Guide.G16;
	    case OffsetFRmer: return G17;
	    }
	    return 0;
	}
    }

    // attributes

    private String id;                          // Guide seedPrefix/suppPrefix[PCTId.g2.g12]
    private String location;                    // Guide location "g2.g12"
    private String DBId;                        // specific database ID, such as MIMAT
    private CodingTranscript CT;                // coding transcript
    private String CTId;                        // coding transcript database name, ENST
    private String sequence;                    // sequence of the guide RNA
    private String MRE;                         // MRE sequence (CT sequence from t1 to tlast), or kmer
    private boolean A1;                         // last nt in MRE is A ? (facing g1 in seed)
    private char t1Nucleotide;                  // nucleotide at t1
    private int g1;                             // g1 position
    private int g2;                             // g2 position
    private int g5;
    private int g13;
    private int g16;
    private int g12;                            // g12 position
    private String seedPrefix;                  // seed prefix used to build it
    private String suppRegion;                  // supp region used to build it
    private int bridgeLength;                   // bridge length
    private Seed seed;                          // seed interaction type
    private String interactionType;             // interactio type as a string
    private String siteType;                    // site type as a string (Bartel's nomenclature)
    private double Kd;                          // Kd derived from interpolated RNA Bind-n-Seq data (Zamore's lab)
    private Supplementary supplementary;        // supplementary type
    private String supplementaryType;           // supplementary type as a string
    private int t1;                             // t1 position in CT
    private int tlast;                          // tlast position in CT
    private Duplex seedDuplex;                  // miniduplex from seed initiation site
    private Duplex suppDuplex;                  // miniduplex from supp initiation site
    private Duplex duplex;                      // duplex for this guide on target
    private boolean folded = false;             // folded status
    private boolean hasSupp;                    // whether we have a grip for the supp or not - affect fold()
    private boolean earlyTargetBulge = false;   // there is a bulge in the target between g1 and g2
    private boolean hasMoreThanOneSeedBulge = false;
    private boolean hasLongerThanOneNtSeedBulge = false;
    private String alignedTarget;               // aligned target in pretty-print format
    private String alignedPairs;                // aligned base pairs in the pretty-print format
    private String alignedGuide;                // aligned guide in pretty-print format
    private double phastConsScore;              // phastCons alignment score
    private double seedAccessibility;           // seed accessibility
    private double suppAccessibility;           // supp accessibility

    // getters
    @JsonProperty( "gid" )
    public String getDBId()                          { return this.DBId; }
    @JsonProperty( "gna" )
    public String getId()                            { return this.id; }
    @JsonIgnore
    public String getLocation()                      { return this.location; }
    @JsonProperty( "gseq" )
    public String getSequence()                      { return this.sequence; }
    @JsonIgnore
    public String getSeed()                          { return this.sequence.substring( 1, 8 ); }
    @JsonIgnore
    public String getSupplementaryRegion()           { return this.sequence.substring( 12, 17 ); }
    @JsonIgnore
    public CodingTranscript getCT()                  { return this.CT; }
    @JsonProperty( "tid" )
    public String getTranscriptId()                  { return this.CTId; }
    @JsonIgnore
    public String getGeneId()                        { return this.getCT().getId(); }
    @JsonProperty( "trn" )
    public String getCTName() {
	if( this.CT != null ) return this.CT.getName();
	else return "kmer";
    }
    @JsonIgnore
    public String getSequence( int start, int stop ) { return this.sequence.substring( start, stop ); }
    @JsonIgnore
    public boolean isSeedless()                      { return this.seed != Seed.SEEDLESS; }
    @JsonIgnore
    public Seed getSeedType()                        { return this.seed; }
    @JsonIgnore
    public boolean hasEarlyBulgeInTarget()           { return this.earlyTargetBulge; }
    @JsonIgnore
    public boolean initiate()                        { return this.seed.ordinal() <= Seed.DistantFRmer.ordinal(); } // initiation with 4 bps and min of 2 nts in the g2-g5,
    @JsonIgnore
    public boolean bind()                            { return this.seedBind() || this.suppBind(); } // bind by the seed or supp region
    @JsonIgnore
    public boolean seedBind()                        { return this.seed.ordinal() < Seed.SEEDLESS.ordinal(); } // bind by the seed, minimum DistantFRmer
    @JsonProperty( "isc" )
    public int getSeedScore()                        { return this.seed.ordinal(); }
    @JsonProperty( "ity" )
    public String getInteractionType()               { return this.interactionType; }
    @JsonProperty( "sty" )
    public String getSiteType()                      { return this.siteType; }
    @JsonProperty( "sen" )
    public double getSeedEnergy()                    { return this.seedDuplex.getFreeEnergy(); }
    @JsonIgnore
    public boolean isSupplementless()                { return this.supplementary != Supplementary.SUPPLEMENTLESS; }
    @JsonIgnore
    public boolean suppBind()                        { return this.supplementary.ordinal() <= Supplementary.OffsetFRmer.ordinal(); } // bind by the supp region, minimum 4 bps
    @JsonProperty( "sut" )
    public String getSuppType()                      { return this.supplementaryType; }
    @JsonProperty( "sue" )
    public double getSuppEnergy()                    { return this.suppDuplex.getFreeEnergy(); }
    @JsonIgnore
    public int getLength()                           { return this.sequence.length(); }
    @JsonIgnore
    public int getg1()                               { return this.g1; }
    @JsonIgnore
    public int getg2()                               { return this.g2; }
    @JsonIgnore
    public int getg5()                               { return this.g5; }
    @JsonIgnore
    public int getg12()                              { return this.g12; }
    @JsonIgnore
    public int getg13()                              { return this.g13; }
    @JsonIgnore
    public int getg16()                              { return this.g16; }
    @JsonIgnore
    public String getSeedPrefix()                    { return this.seedPrefix; }
    @JsonIgnore
    public String getSuppPrefix()                    { return this.suppRegion; }
    @JsonIgnore
    public String getSuppRegion()                    { return this.suppRegion; }
    @JsonIgnore
    public String getGrip()                          { return this.seedPrefix + "/" + this.suppRegion; }
    @JsonProperty( "brl" )
    public int getBridgeLength()                     { return this.bridgeLength; }
    @JsonProperty( "pst" )
    public int gett1()                               { return this.t1; }
    @JsonProperty( "pen" )
    public int gettlast()                            { return this.tlast; }
    @JsonProperty( "tseq" )
    public String getMRE()                           { return this.MRE; }
    @JsonProperty( "rgn" )
    public String getRegion() {
	if( this.CT != null ) return this.CT.printRegions( this.CT.determineRegion( this.gett1(), this.gettlast() ) );
	else return "N/A";
    }
    @JsonIgnore
    public Duplex getDuplex()                        { return this.duplex; }
    @JsonProperty( "den" )
    public Double getDeltaG()                        { return this.duplex.getFreeEnergy(); }
    @JsonProperty( "kdv" )
    public double getKd()                            { return this.Kd; }
    @JsonIgnore
    public boolean isFolded()                        { return this.folded; }
    @JsonProperty( "pcs" )
    public double getPhastConsScore()                { return this.phastConsScore; }
    @JsonProperty( "sac" )
    public double getSeedAccessibility()             { return this.seedAccessibility; }
    @JsonProperty( "suc" )
    public double getSuppAccessibility()             { return this.suppAccessibility; }

    // duplex structure":["guggcaguaacagugauagugggacaugccA","      !!|||:! |:||! |!:|| |||||","------uguugGU-CGAUu-cuGUG-ACGGu"]
    @JsonProperty( "alt" )
    public String getAlignedTarget()                 { return this.alignedTarget; }
    @JsonProperty( "alp" )
    public String getAlignedPairs()                  { return this.alignedPairs; }
    @JsonProperty( "alg" )
    public String getAlignedGuide()                  { return this.alignedGuide; }
    @JsonIgnore
    public List<String> getStructure() {
	Triplet parts = this.prettyGuidePrint();
	List<String> result = new ArrayList<>();
	result.add( parts.getS1() );
	result.add( parts.getS2() );
	result.add( parts.getS3() );
	return result;
    }

    // return the #char from the right in the target facing gi >= 0
    // -----uguugGUCGAUucu----GUGACGGu
    public int getTargetIndex( int gi ) {
	StringBuilder sb = new StringBuilder( this.alignedGuide );
        String revAlignedGuide = sb.reverse().toString();
	StringBuilder sbT = new StringBuilder( this.alignedTarget );
        String revAlignedTarget = sbT.reverse().toString();
	int i = 0;
	int dashes = 0;
	int tDashes = 0;
	while( i-dashes != gi ) {
	    if( revAlignedGuide.charAt( i ) == '-' ) dashes++;
	    if( revAlignedTarget.charAt( i ) == '-' ) tDashes++;
	    i++;
	}
	return i-tDashes;
    }

    // return "g1.g12"
    private static String getLocation( String id ) {
	//System.out.println( "getLocation( " + id + " )" );
	int brakOpen = id.indexOf( "[" );
	int brakClos = id.indexOf( "]", brakOpen + 1 );
	//System.out.println( "brakOpen: " + brakOpen + ", brakClos: " + brakClos );
	return id.substring( brakOpen + 1, brakClos );	
    }

    public static String replaceLeadingDots(String s) {
	if( s == null || s.isEmpty() ) return s;
	int i = 0;
	while( i < s.length() && s.charAt(i) == '.' ) i++;
	if( i == 0 ) return s; // no leading dots
	return "p".repeat( i ) + s.substring( i );
    }


    // setters
    public void setPhastConsScore( double score ) { this.phastConsScore = score; }
    public void setSeedAccessibility( double access ) { this.seedAccessibility = access; }
    public void setSuppAccessibility( double access ) { this.suppAccessibility = access; }


    /** Seed type is determined from this.duplex.getAntiStrandState and tlast nt in MRE
                           v
     UAUAGACACACACACAUAUGCAC  GUGCAUAUGUGUGUGUCUAUA
     ((((((((((((..(((((((((  ))))))))))))))))))))), -54.224
                               ^^^^^^^   ^^^^^^

     ((((((())))))) => 8-mer or 7-mer-m8 (depends on A1)
     ?(((((())))))? => 7-mer-A1 or 6-mer (depends on A1)
     ((((((??)))))) => 6-mer-m8

    */

    // Constructor
    public Guide( String id, String DBId, CodingTranscript CT, String CTId, String sequence, int g2, int g13, String seedPrefix, String suppRegion, int bridgeLength, int t1, int tlast, boolean hasSupp ) {
	// Utils.debug(
	// 	    "Guide( " +
	// 	    id + ", " +
	// 	    DBId + ", " +
	// 	    CT.getName() + ", " +
	// 	    CTId + ", " +
	// 	    sequence + ", " +
	// 	    "g2: " + g2 + ", " +
	// 	    "g13: " + g13 + ", " +
	// 	    seedPrefix + ", " +
	// 	    suppRegion + ", " +
	// 	    bridgeLength + ", " +
	// 	    t1 + ", " +
	// 	    tlast + ", " +
	// 	    hasSupp );
		     
	this.id = id;
	this.DBId = DBId;
	this.CT = CT;
	this.CTId = CTId;
	this.sequence = sequence;
	this.g2 = g2;
	this.g1 = this.g2 + 1;
	this.g5 = this.g2 - 3;
	this.g13 = g13;
	this.g12 = this.g13 + 1;
	this.g16 = this.g13 - 3;
	this.location = this.g12 + "." + this.g1; // "g12.g1" those ancored in RISC
	this.seedPrefix = seedPrefix;
	this.suppRegion = suppRegion;
	this.bridgeLength = bridgeLength;
	this.t1 = t1;
	this.tlast = tlast;
	this.Kd = SEEDLESSKD;
	this.hasSupp = hasSupp;
	this.MRE = this.CT.getSequence().getSequence( this.t1, this.tlast+1 );
	//Utils.debug( "completing Guide: g1: " + this.g1 + ", g5: " + this.g5 + ", g16: " + this.g16 + ", location: " + this.location + ", MRE: " + this.MRE );
    }

    public Guide( String id, String DBId, String MRE, String sequence, int g2, int g13, String seedPrefix, String suppRegion, int bridgeLength, int t1, int tlast, boolean hasSupp ) {
	this.id = id;
	this.DBId = DBId;
	this.MRE = MRE;
	this.CT = null;
	this.CTId = "N/A";
	this.sequence = sequence;
	this.g2 = g2;
	this.g1 = this.g2 + 1;
	this.g5 = this.g2 - 3;
	this.g13 = g13;
	this.g12 = this.g12 + 1;
	this.g16 = this.g13 - 3;
	this.location = this.g12 + "." + this.g1; // "g12.g1" those anchored in RISC
	this.seedPrefix = seedPrefix;
	this.suppRegion = suppRegion;
	this.bridgeLength = bridgeLength;
	this.t1 = t1;
	this.tlast = tlast;
	this.hasSupp = hasSupp;
	this.Kd = 1000;
    }

    public Guide( Guide other, String id, String DBId, String sequence ) {
	this.id = id;
	this.location = this.getLocation( this.id );
	this.DBId = DBId;
	this.sequence = sequence;
	this.CT = other.getCT();
	this.CTId = "N/A";
	this.g2 = other.getg2();
	this.g5 = this.g2 - 3;
	this.g13 = other.getg13();
	this.g12 = this.g13 + 1;
	this.g16 = this.g13 - 3;
	this.seedPrefix = other.getSeedPrefix();
	this.suppRegion = other.getSuppRegion();
	this.bridgeLength = other.getBridgeLength();
	this.t1 = other.gett1();
	this.tlast = other.gettlast();
	this.Kd = 1000;
    }

    public Guide mutate( int pos, String kMer, String id ) {
	String sequence = this.getSequence( 0, pos ) + kMer + this.getSequence( pos + kMer.length(), this.getLength() );
	Guide result = new Guide( id, this.getDBId(), this.getCT(), this.getGeneId(), sequence, this.g2, this.g12, this.getSeedPrefix(), this.getSuppRegion(), this.getBridgeLength(), this.gett1(), this.gettlast(), this.hasSupp );
	return result;
    }

    public Guide mutate( String sequence, String id ) {
	Guide result = new Guide( id, this.getDBId(), this.getCT(), this.getGeneId(), sequence, this.g2, this.g12, this.getSeedPrefix(), this.getSuppRegion(), this.getBridgeLength(), this.gett1(), this.gettlast(), this.hasSupp );
	return result;
    }

    private void adjustBridgeLength() {
	// bridge is defined by the strand between last bp in seed (starting at g8 going down to g2) and first bp in supp (from g13 to g17)
	//MiR	int oldBridgeLength = this.bridgeLength;
	this.bridgeLength = this.duplex.bpDistance( 7, 12 ); // g8, g13
	//if( oldBridgeLength != this.bridgeLength )
	//System.out.println( "new bridge length: " + oldBridgeLength + " => " + this.bridgeLength );
    }

    private static double penalizeWobble( int position ) {
	//Utils.debug( "penalizeWobble( " + position + " )" );
	switch( position ) {
	case( 1 ): // g2
	    return Math.pow( 2, 6.407 );
	case( 2 ): // g3
	    return Math.pow( 2, 7.040 );
	case( 3 ): // g4
	    return Math.pow( 2, 7.432 );
	case( 4 ): // g5
	    return Math.pow( 2, 6.451 );
	case( 5 ): // g6
	    return Math.pow( 2, 6.200 );
	case( 6 ): // g7
	    return Math.pow( 2, 4.094 );
	case( 7 ): // g8
	    return Math.pow( 2, 2.234 );
	}
	return Math.pow( 2, 6.632 ); // average
    }

    private static double penalizeDeletion( int position ) {
	//Utils.debug( "penalizeDeletion( " + position + " )" );
	switch( position ) {
	case 1: // g2
	    return Math.pow( 2, 5.004 );
	case 2: // g3
	    return Math.pow( 2, 6.038 );
	case 3: // g4
	    return Math.pow( 2, 7.737 );
	case 4: // g5
	    return Math.pow( 2, 8.202 );
	case 5: // g6
	    return Math.pow( 2, 8.067 );
	case 6: // g7
	    return Math.pow( 2, 4.748 );
	}
	return Math.pow( 2, 6.317 ); // average
    }

    private static double penalizeMismatch( int position ) {
	//Utils.debug( "penalizeMismatch( g" + (position+1) + " )" );
	switch( position ) {
	case 1: // g2
	    return Math.pow( 2, 6.266 );
	case 2: // g3
	    return Math.pow( 2, 8.323 );
	case 3: // g4
	    return Math.pow( 2, 7.644 );
	case 4: // g5
	    return Math.pow( 2, 7.486 );
	case 5: // g6
	    return Math.pow( 2, 8.050 );
	case 6: // g7
	    return Math.pow( 2, 5.660 );
	case 7: // g8
	    return Math.pow( 2, 2.998 );
	}
	return Math.pow( 2, 7.389 ); // average
    }

    private static double penalizeBulge( int position ) {
	//Utils.debug( "penalizeBulge( g" + (position+1) + ".g" + (position+2) + " )" );
	switch( position ) {
	case 1: // g2.g3
	    return Math.pow( 2, 4.992 ); // interpolated
	case 2: // g3.g4
	    return Math.pow( 2, 5.995 );
	case 3: // g4.g5
	    return Math.pow( 2, 6.998 );
	case 4: // g5.g6
	    return Math.pow( 2, 6.875 );
	case 5: // g6.g7
	    return Math.pow( 2, 4.814 );
	case 6: // g7.g8
	    return Math.pow( 2, 2.753 ); // interpolated
	}
	return Math.pow( 2, 6.170 ); // average
    }

    // k is the number of bps, gFirst is the first g that pairs in the g1-g8 region
    //    gFirst is in G2, G3, ...
    //    gLast is in G8, G7, ...
    //    gLast and gFirst are separated by at least 3
    private static Seed adjustKmer( int gFirst, int gLast, boolean A1 ) {
	Seed seed = Seed.SEEDLESS;
	int k = gLast - gFirst + 1;
	switch( k ) {
	case 7:
	    if( A1 ) seed = Seed.A1SVmer;
	    else seed = Seed.SVmer;
	    break;
	case 6:
	    if( gFirst == Guide.G2 - 1 ) {
		if( A1 ) seed = Seed.A1SXmer;
		else seed = Seed.SXmer;
	    }
	    else // gFirst == G3
		seed = Seed.OffsetSXmer;
	    break;
	case 5:
	    if( gFirst == Guide.G2 - 1 ) {
		if( A1 ) seed = Seed.A1FVmer;
		else seed = Seed.FVmer;
	    }
	    else if( gFirst == Guide.G3 - 1 ) seed = Seed.OffsetFVmer;
	    else seed = Seed.DistantFVmer;
	    break;
	case 4:
	    if( gFirst == Guide.G2 - 1 ) {
		if( A1 ) seed = Seed.A1FRmer;
		else seed = Seed.FRmer;
	    }
	    else if( gFirst == Guide.G3 - 1 ) seed = Seed.OffsetFRmer;
	    else seed = Seed.DistantFRmer;
	    break;
	}
	return seed;
    }

    // k is the number of bps, gFirst is the first g that pairs in the g1-g8 region
    //    gFirst is in G2, G3, ...
    //    gLast is in G8, G7, ...
    //    gLast and gFirst are separated by at least 3
    private void adjustKmer( int gFirst, int gLast ) {
	int k = gLast - gFirst + 1;
	switch( k ) {
	case 7:
	    if( this.A1 ) this.seed = Seed.A1SVmer;
	    else this.seed = Seed.SVmer;
	    break;
	case 6:
	    if( gFirst == Guide.G2 ) {
		if( this.A1 ) this.seed = Seed.A1SXmer;
		else this.seed = Seed.SXmer;
	    }
	    else // gFirst == G3
		this.seed = Seed.OffsetSXmer;
	    break;
	case 5:
	    if( gFirst == Guide.G2 ) {
		if( this.A1 ) this.seed = Seed.A1FVmer;
		else this.seed = Seed.FVmer;
	    }
	    else if( gFirst == Guide.G3 ) this.seed = Seed.OffsetFVmer;
	    else this.seed = Seed.DistantFVmer;
	    break;
	case 4:
	    if( gFirst == Guide.G2 ) {
		if( this.A1 ) this.seed = Seed.A1FRmer;
		else this.seed = Seed.FRmer;
	    }
	    else if( gFirst == Guide.G3 ) this.seed = Seed.OffsetFRmer;
	    else this.seed = Seed.DistantFRmer;
	    break;
	}
    }

    private static String setSiteType( Seed seed, int gFirstAdjusted, int gLastAdjusted, boolean A1, String interactionName, String postfix, boolean bulgeInGuide, boolean bulgeInTarget ) {
	// follow Bartel's nomenclature for site types: can be 5mer to 8mer
	// Target sequence: t24 (or further) to t31 (t31 if A1)
	//     8mer t24-t31, k = 7, "8mer" (initiative)
	//     7mer t25-t31, k = 6, "7mer-A1" (initiative)
	//     7mer t24-t30, k = 7, "7mer-m8" (initiative)
	//     6mer t26-t31, k = 5, "6mer-A1" (initiative)
	//     6mer t25-t30, k = 6, "6mer" (initiative)
	//     6mer t24-t29, k = 6, "6mer-m8" (offset)
	//     5mer t27-t31, k = 5, "5mer" (initiative)
	//     5mer t26-t30, k = 5, "5mer-m3.7" (offset)
	//     5mer t25-t29, k = 5, "5mer-m4.8" (distant)
	//     seedless, "seedless"
	String kmer = "seedless";
	boolean offset = seed.interactionName().contains( "offset" );
	boolean distant = seed.interactionName().contains( "distant" );
	boolean initiative = !offset && !distant;
	int k = gLastAdjusted - gFirstAdjusted + 1; // kmer base k (not considering A1)
	if( bulgeInGuide ) k--; // remove one nt in the site kmer because of a bulge in the seed
	if( bulgeInTarget ) k++; // add one nt in the site kmer because of a bulge in the target
	if( initiative ) {
	    if( k == 8 ) {
		if( A1 ) kmer = "9mer";
		else kmer = "8mer";
	    }
	    else if( k == 7 ) {
		if( A1 ) kmer = "8mer";
		else kmer = "7mer-m8";
	    }
	    else if( k == 6 ) {
		if( A1 ) kmer = "7mer-A1";
		else kmer = "6mer";
	    }
	    else if( k == 5 ) {
		if( A1 ) kmer = "6mer-A1";
		else kmer = "5mer-m2.6";
	    }
	    else kmer = "seedless";
	}
	else if( offset ) {
	    if( k == 7 ) kmer = "7mer-m8";
	    else if( k == 6 ) kmer = "6mer-m8";
	    else if( k == 5 ) kmer = "5mer-m3.7";
	    //else if( k == 5 ) kmer = "seedless";
	    else kmer = "seedless";
	}
	else if( distant ) {
	    if( k == 5 ) kmer = "5mer-m4.8";
	    //if( k == 5 ) kmer = "seedless";
	    else kmer = "seedless";
	}
	return kmer;
    }

    private void setSiteType( int gFirstAdjusted, int gLastAdjusted, boolean A1, String interactionName, String postfix, boolean bulgeInGuide, boolean bulgeInTarget ) {
	// follow Bartel's nomenclature for site types: can be 5mer to 8mer
	// Target sequence: t24 (or further) to t31 (t31 if A1)
	//     8mer t24-t31, k = 7, "8mer" (initiative)
	//     7mer t25-t31, k = 6, "7mer-A1" (initiative)
	//     7mer t24-t30, k = 7, "7mer-m8" (initiative)
	//     6mer t26-t31, k = 5, "6mer-A1" (initiative)
	//     6mer t25-t30, k = 6, "6mer" (initiative)
	//     6mer t24-t29, k = 6, "6mer-m8" (offset)
	//     5mer t27-t31, k = 5, "5mer" (initiative)
	//     5mer t26-t30, k = 5, "5mer-m3.7" (offset)
	//     5mer t25-t29, k = 5, "5mer-m4.8" (distant)
	//     seedless, "seedless"
	String kmer = "seedless";
	boolean offset = this.seed.interactionName().contains( "offset" );
	boolean distant = this.seed.interactionName().contains( "distant" );
	boolean initiative = !offset && !distant;
	int k = gLastAdjusted - gFirstAdjusted + 1; // kmer base k (not considering A1)
	if( bulgeInGuide ) k--; // remove one nt in the site kmer because of a bulge in the seed
	if( bulgeInTarget ) k++; // add one nt in the site kmer because of a bulge in the target
	if( initiative ) {
	    if( k == 8 ) {
		if( A1 ) kmer = "9mer";
		else kmer = "8mer";
	    }
	    else if( k == 7 ) {
		if( A1 ) kmer = "8mer";
		else kmer = "7mer-m8";
	    }
	    else if( k == 6 ) {
		if( A1 ) kmer = "7mer-A1";
		else kmer = "6mer";
	    }
	    else if( k == 5 ) {
		if( A1 ) kmer = "6mer-A1";
		else kmer = "5mer-m2.6";
	    }
	    else kmer = "seedless";
	}
	else if( offset ) {
	    if( k == 7 ) kmer = "7mer-m8";
	    else if( k == 6 ) kmer = "6mer-m8";
	    else if( k == 5 ) kmer = "5mer-m3.7";
	    //else if( k == 5 ) kmer = "seedless";
	    else kmer = "seedless";
	}
	else if( distant ) {
	    if( k == 5 ) kmer = "5mer-m4.8";
	    //if( k == 5 ) kmer = "seedless";
	    else kmer = "seedless";
	}
	this.siteType = kmer;
    }

    private void adjustSuppKmer( int gFirst, int gLast ) {
	int k = gLast - gFirst + 1;
	switch( k ) {
	case 5:
	this.supplementary = Supplementary.FVmer;
	break;
	case 4:
	    if( gFirst == Guide.G13 ) this.supplementary = Supplementary.FRmer;
	    else this.supplementary = Supplementary.OffsetFRmer;
	}
    }

    // determine Seed and Kd of a seed duplex that is missing g1 and t1
    public static double determineKd( Duplex duplex, boolean A1 ) {
	Seed seed = Seed.SEEDLESS; // assume seedless
	String siteType = "seedless";
	double penalty = 1.0; // assume no penalty
	String prefix = "";
	if( A1 ) prefix = "A1-";
	String postfix = "";
	String sitePostfix = "";
	double kd = SEEDLESSKD;

	// seed duplex does not contains g1 => adjust indexes
	if( duplex.getNumberOfCanonical( Guide.G2 - 1, Guide.G8 - 1 ) >= 6 ) // 7mer
	    if( A1 ) seed = Seed.A1SVmer;
	    else seed = Seed.SVmer;

	// locate the 4mer grip and determine seed type from it
	if( duplex.isCanonicallyPaired( Guide.G2 - 1, Guide.G5 - 1 ) ) { // core 4mer passed
	    if( duplex.isCanonicalGU( Guide.G8 - 1 ) ) { // 7mer
		if( A1 ) seed = Seed.A1SVmer;
		else seed = Seed.SVmer;
	    }
	    else if( duplex.isCanonicalGU( Guide.G7 - 1 ) ) { // 6mer
		if( A1 ) seed = Seed.A1SXmer;
		else seed = Seed.SXmer;
	    }
	    else if( duplex.isCanonicalGU( Guide.G6 - 1 ) ) { // 5mer
		if( A1) seed = Seed.A1FVmer;
		else seed = Seed.FVmer;
	    }
	    else {
		if( A1 ) seed = Seed.A1FRmer;
		else seed = Seed.FRmer;
	    }
	}
	else if( duplex.isCanonicallyPaired( Guide.G3 - 1, Guide.G6 - 1 ) ) { // core offset 4mer passed
	    if( duplex.isCanonicalGU( Guide.G8 - 1 ) ) // 6mer
		seed = Seed.OffsetSXmer;
	    else if( duplex.isCanonicalGU( Guide.G7 - 1 ) ) // 5mer
		seed = Seed.OffsetFVmer;
	    else seed = Seed.OffsetFRmer;
	}
	else if( duplex.isCanonicallyPaired( Guide.G4 - 1, Guide.G7 - 1 ) ) { // core distant 4mer passed
	    if( duplex.isCanonicalGU( Guide.G8 - 1 ) ) // 5mer
		seed = Seed.DistantFVmer;
	    else seed = Seed.DistantFRmer;
	}

	//Utils.debug( "duplex:\n" + duplex );
	//Utils.debug( "seed = " + seed );

	// adjust kd
	if( seed != Seed.SEEDLESS ) { // if we have a grip

	    // adjust gLast in the case there are multiple bulges or a long bulge 3' of the basic 4mer
	    int gLastAdjusted = duplex.adjustGLast( ( seed.getGFirst() - 1 ) + 3, ( seed.getGLast() - 1 ) );
	    //Utils.debug( "determineKd: gLastAdjusted = " + gLastAdjusted );
	    //Utils.debug( "iterating from " + (Guide.G2 - 1) + " to " + (Guide.G8 - 1) );

	    // identify gFirstAdjusted
	    int gFirstAdjusted = 0;
	    for( int i = Guide.G2 - 1; i <= Guide.G8 - 1; i++ )
		if( duplex.isCanonicalGU( i ) ) {
		    gFirstAdjusted = i;
		    break; // most 5' WC+GU identified, can exit the loop
		}

	    //Utils.debug( "determineKd: gFirstAdjusted = " + gFirstAdjusted );

	    // adjust seed type at the light of the adjusted gFirst and gLast
	    // gFirstAdjusted = Math.max( gFirstAdjusted, Guide.G2 );
	    seed = adjustKmer( gFirstAdjusted, gLastAdjusted, A1 );

	    boolean bulgeInGuide = false;
	    boolean bulgeInTarget = false;
	    for( int i = gFirstAdjusted; i <= gLastAdjusted; i++ ) {
		if( i < gLastAdjusted ) { // check for bulge first
		    // check for a bulge in target
		    String bulge = duplex.bulgeOfAdjacentPaired( i, i+1 );
		    //Utils.debug( "bulge of adjacent pairs( " + i + ", " + (i+1) + " ) = " + bulge );
		    if( bulge.length() > 0 ) { // there is a singel-nt bulge (no other possibilities)
			bulgeInTarget = true;
			//Utils.debug( "bulge in target!" );
			postfix += "b" + "(" + (i+2) + "." + (i+3) + ")"; // bulge
			sitePostfix += "b?" + "(" + (i+2) + "." + (i+4) + ")";
			penalty *= penalizeBulge( i+1 ); // add 1 to penalize function index, 'coz SeedDuplexes starts at g2, not g1
			//Utils.debug( "penalty now: " + penalty + ", after calling penalizeBulge( " + (i+1) + " )" );
		    }
		}
		if( duplex.isMismatch( i ) ) {
		    postfix += "x" + (i+2); // non-canonical
		    sitePostfix += "x?";
		    penalty *= penalizeMismatch( i+1 ); // add 1 to penalize function index, 'coz SeedDuplexes starts at g2, not g1
		    //Utils.debug( "penalty now: " + penalty + ", after calling penalizeMismatch( " + (i+1) + " )" );
		}
		else if( !duplex.isPaired( i ) ) {
		    bulgeInGuide = true;
		    postfix += "b" + (i+1); // unpaired (bulge in guide)
		    penalty *= penalizeDeletion( i+1 ); // add 1 to penalize function index, 'coz SeedDuplexes starts at g2, not g1
		    //Utils.debug( "penalty now: " + penalty + ", after calling penalizeDeletion( " + (i+1) + " )" );
		}
		else if( duplex.isGU( i ) ) {
		    postfix += "w" + (i+1); // wobble
		    sitePostfix += "w" + (i+1);
		    penalty *= penalizeWobble( i+1 ); // add 1 to penalize function index, 'coz SeedDuplexes starts at g2, not g1
		    //Utils.debug( "penalty now: " + penalty + ", after calling penalizeWobble( " + (i+1) + " )" );
		}
	    }

	    siteType = setSiteType( seed, gFirstAdjusted, gLastAdjusted, A1, seed.interactionName(), postfix, bulgeInGuide, bulgeInTarget );
	    if( !siteType.equals( "seedless" ) && !sitePostfix.equals( "" ) )
		if( !siteType.contains( "-" ) ) siteType += "-" + sitePostfix;
		else siteType += sitePostfix;
	    
	    kd = seed.getKd();
	    //Utils.debug( "adjusted seed: " + seed + ", kd = " + kd + ", penalty = " + penalty + ", final kd = " + kd * penalty );
	    kd *= penalty; // adjust kd

	    // add the dash at the beginning of postfix
	    if( !postfix.equals( "" ) ) postfix = "-" + postfix;
	}
	else prefix = "";

	// addup the seed type
	String interactionType = prefix + seed.interactionName() + postfix;
	//Utils.debug( "duplex =\n" + duplex + "\nreturned seed = " + seed + ", kd = " + kd + ", interaction type: " + interactionType + ", site type: " + siteType );
	return kd;
    }
    
    public void determineSeedSuppTypes() {

	this.seed = Seed.SEEDLESS; // assume seedless
	this.siteType = "seedless";
	double penalty = 1.0; // assume no penalty
	String prefix = "";
	if( this.A1 )
	    prefix = "A1-";  // init prefix based on A1
	String postfix = ""; // assume no postfix
	String sitePostfix = ""; // assume no postfix

	//Utils.debug( "determineSeedSuppTypes: this.duplex =\n" + this.duplex );

	if( !this.earlyTargetBulge ) {

	    // treat special 7mer with 1 default at g4 or g5

	    if( this.duplex.getNumberOfCanonical( Guide.G2, Guide.G8 ) >= 6 ) // 7mer
		if( this.A1 ) this.seed = Seed.A1SVmer;
		else this.seed = Seed.SVmer;

	    // locate the 4mer grip and determine seed type from it
	    if( this.duplex.isCanonicallyPaired( Guide.G2, Guide.G5 ) ) { // core 4mer passed
		if( this.duplex.isCanonicalGU( Guide.G8 ) ) { // 7mer
		    if( this.A1 ) this.seed = Seed.A1SVmer;
		    else this.seed = Seed.SVmer;
		}
		else if( this.duplex.isCanonicalGU( Guide.G7 ) ) { // 6mer
		    if( this.A1 ) this.seed = Seed.A1SXmer;
		    else this.seed = Seed.SXmer;
		}
		else if( this.duplex.isCanonicalGU( Guide.G6 ) ) { // 5mer
		    if( this.A1) this.seed = Seed.A1FVmer;
		    else this.seed = Seed.FVmer;
		}
		else {
		    if( this.A1 ) this.seed = Seed.A1FRmer;
		    else this.seed = Seed.FRmer;
		}
	    }
	    else if( this.duplex.isCanonicallyPaired( Guide.G3, Guide.G6 ) ) { // core offset 4mer passed
		if( this.duplex.isCanonicalGU( Guide.G8 ) ) // 6mer
		    this.seed = Seed.OffsetSXmer;
		else if( this.duplex.isCanonicalGU( Guide.G7 ) ) // 5mer
		    this.seed = Seed.OffsetFVmer;
		else this.seed = Seed.OffsetFRmer;
	    }
	    else if( this.duplex.isCanonicallyPaired( Guide.G4, Guide.G7 ) ) { // core distant 4mer passed
		if( this.duplex.isCanonicalGU( Guide.G8 ) ) // 5mer
		    this.seed = Seed.DistantFVmer;
		else this.seed = Seed.DistantFRmer;
	    }
	}
	    
	// complete prefix from g1 pairing case
	// if( this.duplex.getPartner( Guide.G1 ) == null ) prefix += "x-";
	// else
	//     switch( this.duplex.getBasePairType( Guide.G1 ) ) {
	//     case ')':
	// 	prefix += "m-";
	// 	break;
	//     case '}':
	// 	prefix += "w-";
	// 	break;
	//     case '>':
	// 	prefix += "x-";
	// 	break;
	//     }

	if( this.seed != Seed.SEEDLESS ) { // if we have a 4mer grip and no early bulge

	    // adjust gLast in the case there are multiple bulges or a long bulge 3' of the basic 4mer
	    int gLastAdjusted = this.duplex.adjustGLast( this.seed.getGFirst() + 3, this.seed.getGLast() );

	    // identify gFirstAdjusted
	    int gFirstAdjusted = 0;
	    for( int i = Guide.G2; i <= Guide.G8; i++ )
		if( this.duplex.isCanonicalGU( i ) ) {
		    gFirstAdjusted = i;
		    break; // most 5' WC+GU identified, can exit the loop
		}

	    // adjust seed type at the light of the adjusted gFirst and gLast
	    // gFirstAdjusted = Math.max( gFirstAdjusted, Guide.G2 );
	    this.adjustKmer( gFirstAdjusted, gLastAdjusted );
	    //Utils.debug( "determineSeedSupp: gFirstAdjusted = " + gFirstAdjusted + ", gLastAdjusted = " + gLastAdjusted );
	    //Utils.debug( "iterating from " + gFirstAdjusted + " to " + gLastAdjusted );

	    boolean bulgeInGuide = false;
	    boolean bulgeInTarget = false;
	    for( int i = gFirstAdjusted; i <= gLastAdjusted; i++ ) {
		if( i < gLastAdjusted ) { // check for bulge first
		    // check for a bulge in target
		    String bulge = this.duplex.bulgeOfAdjacentPaired( i, i+1 );
		    if( bulge.length() > 0 ) { // there is a singel-nt bulge (no other possibilities)
			bulgeInTarget = true;
			postfix += "b" + "(" + (i+1) + "." + (i+2) + ")"; // bulge
			sitePostfix += "b" + this.MRE.charAt( this.duplex.getBasePairPartner( i ) - 1 ) + "(" + (i+1) + "." + (i+3) + ")";
			penalty *= penalizeBulge( i );
		    }
		}
		if( this.duplex.isMismatch( i ) ) {
		    postfix += "x" + (i+1); // non-canonical
		    sitePostfix += "x" + this.MRE.charAt( this.duplex.getBasePairPartner( i ) ) + (i+1); 
		    penalty *= penalizeMismatch( i );
		    //Utils.debug( "penalty now: " + penalty + ", after calling penalizeMismatch( " + i + " )" );
		}
		else if( !this.duplex.isPaired( i ) ) {
		    bulgeInGuide = true;
		    postfix += "b" + (i+1); // unpaired (bulge in guide)
		    penalty *= penalizeDeletion( i );
		}
		else if( this.duplex.isGU( i ) ) {
		    postfix += "w" + (i+1); // wobble
		    sitePostfix += "w" + (i+1);
		    penalty *= penalizeWobble( i );
		}
	    }

	    this.setSiteType( gFirstAdjusted, gLastAdjusted, this.A1, this.seed.interactionName(), postfix, bulgeInGuide, bulgeInTarget );
	    if( !this.siteType.equals( "seedless" ) && !sitePostfix.equals( "" ) )
		if( !this.siteType.contains( "-" ) ) this.siteType += "-" + sitePostfix;
		else this.siteType += sitePostfix;

	    //Utils.debug( "SeedSupp adjust seed = " + this.seed );
	    this.Kd = this.seed.getKd();
	    //Utils.debug( "SeedSupp adjust kd = " + this.Kd + ", penalty = " + penalty );
	    this.Kd *= penalty; // adjust Kd
	    //Utils.debug( "SeedSupp adjust * penalty = " + this.Kd );

	    // add the dash at the beginning of postfix
	    if( !postfix.equals( "" ) ) postfix = "-" + postfix;
	}
	else prefix = "";

	// addup the seed type
	this.interactionType = prefix + this.seed.interactionName() + postfix;

	// determine supplementary case

	this.supplementary = Supplementary.SUPPLEMENTLESS;
	postfix = "";// reset postfix for supplementary region

	if( this.duplex.isCanonicallyPairedGU( Guide.G13, Guide.G16 ) ) { // core 4mer passed
	    if( this.duplex.isCanonicalGU( G17 ) ) this.supplementary = Supplementary.FVmer;
	    else this.supplementary = Supplementary.FRmer;
	}
	else if( this.duplex.isCanonicallyPairedGU( G14, G17 ) ) // core offset 4mer
	    this.supplementary = Supplementary.OffsetFRmer;

	if( this.supplementary != Supplementary.SUPPLEMENTLESS ) {

	    // adjust gLast in the case there are multiple bulges or a long bulge 3' in the basic 4mer
	    int gLastAdjusted = this.duplex.adjustGLast( this.supplementary.getGFirst() + 3, this.supplementary.getGLast() );

	    // adjust supp type at the light of the adjusted gLast
	    this.adjustSuppKmer( this.supplementary.getGFirst(), gLastAdjusted );

	    // add mismatches, wobbles and bulges
	    
	    for( int i = this.supplementary.getGFirst(); i <= gLastAdjusted; i++ ) {
		if( this.duplex.isNonCanonical( i ) && i < G17 && this.duplex.isCanonicalGU( i+1 ) ) postfix += "x" + (i+1); // non-canonical
		else if( !this.duplex.isPaired( i ) && i < G17 && this.duplex.isCanonicalGU( i+1 ) ) postfix += "b" + (i+1); // unpaired
		else if( this.duplex.isGU( i ) ) postfix += "w" + (i+1); // wobble
		// check for a bulge in target
		if( i < this.supplementary.getGLast() ) {
		    String bulge = this.duplex.bulgeOfAdjacentPaired( i, i+1 );
		    if( bulge.length() > 0 ) // there is a bulge
			postfix += "b" + "(" + (i+1) + "." + (i+2) + ")"; // bulge
		}
	    }

	    // add the dash at the beginning of postfix
	    if( !postfix.equals( "" ) ) postfix = "-" + postfix;
	}

	// addup the supplementary region type
	//Utils.debug( "SeedSupp adjust seed = " + this.seed + ", Kd = " + this.Kd );
	
	this.supplementaryType = this.supplementary.suppTypeName() + postfix;
	if( this.Kd == SEEDLESSKD )
	    this.Kd /= ( 0.25 * -this.seedDuplex.getFreeEnergy() + 0.75 * -this.suppDuplex.getFreeEnergy() );
    }

    public void fold() {
	if( !this.folded ) {
	    if( this.CT != null )
		this.MRE = this.CT.getSequence().getSequence( this.t1, this.tlast+1 );
	    // fold the bipartite guide

	    // A1 and nucleotide at t1
	    this.A1 = this.MRE.charAt( this.MRE.length() - 1 ) == 'A'; // nucleotide at tlast is A?
	    this.t1Nucleotide = this.MRE.charAt( this.MRE.length() - 1 ); // nucleotide at t1

	    // make the guide sections
	    String seed = this.sequence.substring( 0, 8 ); // g2-g8 (indexed 1 to 7; include g1 for folding)
	    String boxA = this.sequence.substring( 8, 11 ); // g9-g11 (indexed 8 to 10)
	    String supp = this.sequence.substring( 12, 17 ); // g13-g17 (indexed 12 to 16)
	    String boxD = this.sequence.substring( 17, this.getLength() ); // g18-

	    // tight seed folding
	    int bipartiteFrontier = 22; // t23 (one further than the nt facing g8)
	    String part1 = this.MRE.substring( bipartiteFrontier, this.MRE.length() ); // seed pairing potential
	    String mask1 = "p".repeat( part1.length() - 1 ) + Duplex.LOOPMASK + ")" + "q".repeat( seed.length() - 2 ); // remove g1 and t1
	    this.seedDuplex = new Duplex( part1.substring( 0, part1.length() - 1 ), seed.substring( 1, seed.length() ), mask1, 3, false, true );

	    // take the duplex that minimize kd
	    double minKd = Guide.SEEDLESSKD;
	    Duplex minKdDuplex = this.seedDuplex;
	    // if( this.seedDuplex.getConformationalSpace() == null ) {
	    // 	System.out.println( "Duplex( " +
	    // 			    part1.substring( 0, part1.length() - 1 ) + ", " +
	    // 			    seed.substring( 1, seed.length() ) + ", " +
	    // 			    mask1 + ", 2, false, true )" );
	    // 	Utils.stop( "found empty conformational space ", 0 );
	    // }
	    for( Duplex d: this.seedDuplex.getConformationalSpace() ) { // iterate all states
		if( !( d.isShifted() || d.isBulgy() ) ) {
		    double kd = Guide.determineKd( d, this.A1 );
		    //Utils.debug( "in fold, considering duplex:\n" + d + "\nkd = " + kd );
		    if( kd < minKd ) {
			minKd = kd;
			minKdDuplex = d;
		    }
		}
	    }
	    this.Kd = minKd;
	    this.seedDuplex = minKdDuplex;
	    
	    //Utils.debug( "in fold, optimizing kd give seedDuplex:\n" + this.seedDuplex + "\nkd = " + this.Kd );

	    // determine bridge's end
	    int bridgeEnd = 30;
	    if( this.hasSupp ) bridgeEnd = this.g13 - this.t1; // g13 position in 31mer
	    else {
		int maxBridgeEnd = 5;
		String halfPart = this.MRE.substring( maxBridgeEnd, bipartiteFrontier );
		String halfMask = "p".repeat( halfPart.length() ) + Duplex.LOOPMASK + ")" + "q".repeat( supp.length() - 1 );
		this.suppDuplex = new Duplex( halfPart, supp, halfMask, 2, false, true ); // fold chamber 2 (unbalanced)
		//Utils.debug( "suppDuplex:\n" + this.suppDuplex );
		bridgeEnd = this.suppDuplex.getPartner( 0 ) + maxBridgeEnd;
		//Utils.debug( "bridgeEnd: " + bridgeEnd );
	    }

	    // tight supp folding
	    //Utils.debug( "bridgeEnd: " + bridgeEnd + ", bipartiteFrontier: " + bipartiteFrontier + ", g13: " + ( this.g13 - this.t1 ) );
	    String part2 = "";
	    if( this.hasSupp ) part2 = this.MRE.substring( bridgeEnd - supp.length() + 1, bridgeEnd + 1 );
	    else part2 = this.MRE.substring( bridgeEnd - supp.length(), bridgeEnd + 1 );
	    String mask2 = "p".repeat( part2.length() ) + Duplex.LOOPMASK + "q".repeat( supp.length() );
	    this.suppDuplex = new Duplex( part2, supp.substring( 0, supp.length() ), mask2, 2, false, false );
	    
	    //Utils.debug( "suppDuplex:\n" + this.suppDuplex );

	    if( this.seedDuplex.getStrandState() == null || this.suppDuplex.getStrandState() == null ) { // one of mcff results wrong
		System.out.println( "something went wrong with folding, suppDuplex or seedDuplex is null!" );
		System.exit( 0 );
	    }

	    // make the masks for bipartite folding

	    String mask3 = "";
	    if( this.hasSupp ) mask3 = "p".repeat( bridgeEnd - 4 ); // 4 = 5 - 1, where 5 is the length of the supp region
	    else mask3 = "p".repeat( bridgeEnd - 5 ); // 5 is the length of the supp region
	    mask3 +=
		this.suppDuplex.getStrandState().replace( '.', 'p' ).replaceAll( "[{<]", "(" ) +  // keep all base pairs from tight supp folding
		"p".repeat( bipartiteFrontier - bridgeEnd - 1 ) + // allow bridge to pair with boxA
		replaceLeadingDots( this.seedDuplex.getStrandState().replaceAll( "[{<]", "(" ) ) + // keep all base pairs from tight seed folding
		"." + Duplex.LOOPMASK + "." + // . and . are for t1 and g1, respectively, and the LOOPMASK for the loop (needed for cis folding)
		this.seedDuplex.getAntiStrandState().replace( '.', 'q' ).replaceAll( "[}>]", ")" ) +
		"qqqq" + // allow boxA + g12 to pair
		this.suppDuplex.getAntiStrandState().replace( '.', 'q' ).replaceAll( "[}>]", ")" ) +
		"q".repeat( this.getLength() - 17 );  // allow boxD to pair

	    // fold full duplex using the parts' folding results as constraints
	    
	    // Global folding based on bipartite mask;
	    this.duplex = new Duplex( this.MRE, this.sequence, mask3, 2, false, false ); // uncomment for bipartite folding

	    // Analyze duplex structure, make sure that if there are bulges in the seed and supp they are at most 1 nt long
	    if( this.duplex == null ) Utils.stop( "in Guide.fold, null duplex: " + this.id + " " + this.sequence, 0 );

	    // check if bulge between g1 and g2
	    if( this.duplex.getStrandState().charAt( 29 ) == '.' ) { // bulge in the target between g1 and g2
		this.earlyTargetBulge = true;
	    }

	    // determine seed and supplementary types of the duplex
	    this.determineSeedSuppTypes();

	    // adjust the bridge length that may have changed during folding
	    this.adjustBridgeLength();

	    // generate pretty print version of the duplex
	    this.prettyGuidePrint();
	    //System.out.println( this.toString() ); // stop for debugging
	    //System.exit( 0 );

	    // indicate this duplex as folded
	    this.folded = true;
	}
    }

    public boolean hasSameLocation( Guide other ) {
	return this.g1 == other.g1 && this.CT == other.CT;
    }

    public static boolean isNotIsolatedWCPair( String dotbLeft, int i ) {
	char left = ' ';
	char righ = ' ';
	if( i - 1 >= 0 ) left = dotbLeft.charAt( i - 1 );
	if( i + 1 < dotbLeft.length() ) righ = dotbLeft.charAt( i + 1 );
	//if( i == dotbLeft.length() - 2 ) righ = ' '; // -1 is last (t1), -2 is before the last (t2)
	return left == '(' || righ == '(' || ( left == '{' && righ == '(' ) || ( righ == '{' && left == '(' );
	//return surrounded.contains( ""+left+righ );
    }

    public static boolean isNotIsolatedWobblePair( String dotbLeft, int i ) {
	char left = ' ';
	char righ = ' ';
	if( i - 1 >= 0 ) left = dotbLeft.charAt( i - 1 );
	if( i + 1 < dotbLeft.length() ) righ = dotbLeft.charAt( i + 1 );
	//if( i == dotbLeft.length() - 2 ) righ = ' '; // -1 is last (t1), -2 is before the last (t2)
	return left == '(' && righ == '(';
	//return surrounded.contains( ""+left+righ );
    }

    public boolean has4merInSeed() {
	return this.duplex.getAntiStrandState().substring( 1, 8 ).contains( "))))" ); // 4 consecutive seed nts form canonical bps
    }
    
    private Triplet prettyGuidePrint() {
	
	String target = this.duplex.getStrand();
	String guide = this.duplex.getAntiStrand();
	String dotbLeft = this.duplex.getStrandState();
	String dotbRight = this.duplex.getAntiStrandState();
	String space = " ";
	String newline = "\n";

	String alignedTarget = "";
	String bars = "";
	String alignedGuide = "";

	/*
	  target: AUAAAACACCCAGCUAGGACCAUUACUGCCA
	  guide:  UGGCAGUGUCUUAGCUGGUUGU
	  dotbLeft:  .....(((<((((((({((....{((((((.
	  dotbRight: .))))))}))})))))))>)))
                     AUAAAACACCCAGCUAGGACCAUUACUGCCA  UGGCAGUGUCUUAGCUGGUUGU
	  fold:      .....(((<((((((({((....{((((((.  .))))))}))})))))))>)))
	  result:
          auaaaacacccagcuaggaccauuacugccA
               |||:|||||||!||    !||||||
	  -----uguugGUCGAuucu----GUGACGGu

	  5'-gucccucuaccugcacccuca----cauccugccA-3' MTA2-201 [2359..2389] (CDS) Target site: 6mer-A1
	          :|:|||:||:          ||  |||||     A1-5mer(0.62) supplementless(0.13) bridge: 6 phastcons: -1.00 dG: -36.90 Kd: 490
	  3'------uguugGUCGA------uucuGU--GACGGu-5' hsa-miR-34a-5p

	  5'-gucccucuaccugcacccucacauccugccA-3' MTA2-201 [2359..2389] (CDS) Target site: 6mer-A1
	          :|:|||:||:          ||  |||||     A1-5mer(0.62) supplementless(0.13) bridge: 6 phastcons: -1.00 dG: -36.90 Kd: 490
	  3'------uguugGUCGA--uucuGU--GACGGu-5' hsa-miR-34a-5p
	 */

	// i index in left side (target)
	// j index in right side (guide)
	int j = dotbRight.length() - 1;
	int i = 0;
	// these positions in lowercase
	Integer lowerCasePos[] = { 1, 9, 10, 11, 12, 18, 19, 20, 21, 22, 23, 24, 25, 26 }; // guide positions in this list starts at 1.
	Set<Integer> lowerCasePositions = new HashSet<>( Arrays.asList( lowerCasePos ) );
	int g = guide.length();
	while( i < target.length() - 1 ) { // process all positions but t1
	    if( dotbLeft.charAt( i ) == '.' ) { // exposed nt in target
		alignedTarget += Character.toLowerCase( target.charAt( i++ ) );
		bars += space; // insert space in the bp bar
		alignedGuide += "-"; // insert - in the guide
	    }
	    else if( SecondaryStructure.bp5.contains( Character.toString( dotbLeft.charAt( i ) ) ) ) {
		if( dotbRight.charAt( j ) == '.' ) {
		    alignedTarget += "-";
		    bars += space;
		    if( lowerCasePositions.contains( g ) )
			alignedGuide += Character.toLowerCase( guide.charAt( j-- ) );
		    else alignedGuide += guide.charAt( j-- );
		    g--;
		}
		else if( dotbLeft.charAt( i ) == '(' ) {
		    alignedTarget += Character.toLowerCase( target.charAt( i++ ) );
		    //if( isNotIsolatedWCPair( dotbLeft, i-1 ) ) bars += "|"; // we have a pair
		    //else bars += " "; // no show of isolated base pairs
		    bars += "|"; // we have a pair
		    if( lowerCasePositions.contains( g ) )
			alignedGuide += Character.toLowerCase( guide.charAt( j-- ) );
		    else alignedGuide += guide.charAt( j-- );
		    g--;
		}
		else {
		    //String pair = "" + target.charAt( i ) + guide.charAt( j );
		    //if( SecondaryStructure.guBps.contains( pair ) )
		    if( dotbLeft.charAt( i ) == '{' )
			//if( isNotIsolatedWobblePair( dotbLeft, i ) ) bars += "!"; // we have a pair
			//else bars += " "; // no show of isolated base pairs
			bars += "!"; // we have a Wobble pair
		    else bars += ":"; // show non-canonical base pairs
		    //else bars += " "; // show mismatches
		    alignedTarget += Character.toLowerCase( target.charAt( i++ ) );
		    if( lowerCasePositions.contains( g ) )
			alignedGuide += Character.toLowerCase( guide.charAt( j-- ) );
		    else alignedGuide += guide.charAt( j-- );
		    g--;
		}
	    }
	}
	this.alignedTarget = alignedTarget + (this.A1 ? target.charAt( target.length() - 1 ) : Character.toLowerCase( target.charAt( target.length() - 1 ) ) ); // add t1
	this.alignedPairs = bars + " ";
	this.alignedGuide = alignedGuide + Character.toLowerCase( guide.charAt( 0 ) ); // add g1
	return( new Triplet( this.alignedTarget, this.alignedPairs, this.alignedGuide ) );
    }

    // assume left dotb of a mcff duplex mode folding and left side longer than the right side
    //    ie. MRE longer than the miR
    // ex) ...<..(((<(<(<(<{({(.....(.(((((  )))))))})}>)>)>.)>)))>
    private String barRep( String infoTop, String infoCtr, String infoBtm ) {
	Triplet structure = prettyGuidePrint();
	return
	    structure.getS1() + " " + infoTop + "\n" +
	    structure.getS2() + " " + infoCtr + "\n" +
	    structure.getS3() + " " + infoBtm;
    }

    @Override
    public int compareTo( Guide other ) {
	return this.getId().compareTo( other.getId() );
    }

    @Override
    public String toString() {
	// Returns a pretty print version of this guide, forming a duplex with target
	//System.out.println( "printing, checking folded status: " + this.folded );
	if( this.duplex == null ) return this.sequence + " duplex for this guide " + this.id + " is null\n";
	String target = this.duplex.getStrand();
	String guide = this.duplex.getAntiStrand();
	String infoTop;
	if( this.CT != null )
	    infoTop = this.CT.getName() + " [" + this.gett1() + ".." + this.gettlast() + "] (" + this.CT.printRegions( this.CT.determineRegion( this.gett1(), this.gettlast() ) ) + ")";
	else infoTop = "kmer" + " [" + this.gett1() + ".." + this.gettlast() + "]";
	String infoCtr =
	    "seed: " +
	    this.interactionType + this.seedDuplex.getFreeEnergy() +
	    " supp: " +
	    this.getSuppType() + this.suppDuplex.getFreeEnergy() +
	    " bridge: " +
	    this.getBridgeLength() +
	    " deltaG: " +
	    this.duplex.getFreeEnergy() +
	    " kcal/mol Kd = " +
	    this.getKd();
	String infoBtm = this.getId();
	String barRep = this.barRep( infoTop, infoCtr, infoBtm );
	return barRep;
    }

    // Example of a Guide duplex
    /**
       AAUAGGA/AUAGA.0.1953.1943 => Guide( PCTName: VIM-209 t1: 1939 tlast: 1960, seed prefix: AAUAGGA, supp region: AUAGA bridge length: 4, seed type: 8-mer )
       CAAGAUAGAUUUGGAAUAGGAA  UUCCUAUUCCAAUCUAUCUUG
       (((((((((((.((((((((((  ))))))))))))))))))))), -56.74

       caagauagauuuggaauaggaA VIM-209 [1939..1960]
       |||||||||| ||||||||||| 8-mer -56.74
       guuCUAUCUa-accUUAUCCUu bridge: 4
                     

     */
}
