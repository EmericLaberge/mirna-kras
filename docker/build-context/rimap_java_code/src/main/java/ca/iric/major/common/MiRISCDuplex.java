/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Arrays;

import java.lang.Runtime;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;

import java.util.LinkedList;

/** example:
//                 [bridge]
//              1         2
//   t:1234567890123456789012345678
//  5'-GGACACCCAGCUAGGACCAUUACUGCCA-3' Sirt1
//       ||| |||||||.||    |||||||
//  3'---UGUUGGUCGAUUCU----GUGACGGU-5' miR-34a
//     g:21098765432109    87654321
//         2         1
// regions:   [supp]       [ seed ]
//   seed type: 8-mer (7 bps + A1)
//   seed: g2-g8
//   bridge: bulge in target between seed and supp, up to 15 nts (Sheu-Gruttadauria et al. 2019 EMBO J)
//   supplementary region (supp) g13-g16, up to 4 bps (Sheu-Gruttadauria et al. 2019 EMBO J)
//
**/

public class MiRISCDuplex {

    // String index starts at 0
    // Duplex does not contain g1
    final static int seedLo = 0; // g2
    final static int seedHi = 6; // g8
    public final static int seedLength = seedHi - seedLo + 1;
    final static int suppLo = 10; // g12
    final static int suppHi = 15; // g17
    final static int guideLo = 0; // g2
    final static int guideHi = 19; // g21
    final static int bridgeMax = 15; // according to xtal structure

    // For MiRDesign
    //   we need 2 base paired region, split by a possible bridge
    //   in the strand: seedComp + bridge + Supp
    //   in the anti-s: seed     +          Supp
    final static int boxALo = 7; // g9
    final static int boxAHi = 9; // g11

    // Yan et al. modules
    final static int ALo = 7; // g9
    final static int AHi = 9; // g11
    final static int BLo = 10; // g12
    final static int BHi = 12; // g14
    final static int CLo = 13; // g15
    final static int CHi = 15; // g17
    final static int DLo = 16; // g18
    final static int DHi = 19; // g21

    public static int seedScore( String seed ) {
	if( seed.equals( "8-mer" ) ) return 0;
	if( seed.equals( "7-mer-m8" ) ) return 1;
	if( seed.equals( "7-mer-A1" ) ) return 2;
	if( seed.equals( "6-mer" ) ) return 3;
	if( seed.equals( "6-mer-offset" ) ) return 4;
	if( seed.equals( "1-bulge" ) ) return 5;
	if( seed.equals( "0-noncanonical" ) ) return 6;
	if( seed.equals( "unknown" ) ) return 7;
	return 8;
    }
    
    // final static Map<String,Integer> seedScore = new HashMap<>() {{
    // 	    put( "8-mer", 0 );
    // 	    put( "7-mer-m8", 1 );
    // 	    put( "7-mer-A1", 2 );
    // 	    put( "6-mer", 3 );
    // 	    put( "6-mer-offset", 4 );
    // 	    put( "1-bulge", 5 );
    // 	    put( "0-noncanonical", 6 );
    // 	    put( "unknown", 7 );
    // 	}};

    // final static Set<String> seedTypes = new HashSet<>() {{
    // 	    add( "8-mer" );
    // 	    add( "7-mer-m8" );
    // 	    add( "7-mer-A1" );
    // 	    add( "6-mer" );
    // 	    add( "6-mer-offset" );
    // 	    add( "1-bulge" );
    // 	}};

    public static final int MRELength = 31;
    public static final int cheapSeedLength = 4;

    public static int count( String s, char c ) {
	return s.length() - s.replace( Character.toString( c ), "" ).length();
    }

    // Get seed type of an RNA duplex and A1,
    //   assumes dot bracket from RNAduplex
    //   needs a minimum of 6 base pairs (6-mer or 6-mer offset))
    // ex)
    // ..(((<(((((((<(((....((((((  )))))))))>)))))))>))) -54.002
    //                              |
    //                              g2
    // Duplex from MiRISCDuplex contains g1 and A1, which must be removed
    //
    public static String determineSeedType( Duplex duplex, boolean A1 ) {
	//System.out.println( "determineSeedType( " + duplex.getState() + ", A1 = " + A1 + " )" );
	// Split at ( " " ) gives: target brackets in [0]; a space in [1]; guide brackets in [2]
	String[] splitDB = duplex.getState().split( " " );
	String targetDotb = splitDB[0].substring( 0, splitDB[0].length() - 1 ); // ex: ..(((<(((((((<(((....(((((( (length = 27) <= remove A1
	String guideDotb = splitDB[2].substring( 1, splitDB[2].length() );  // ex: )))))))))>)))))))>))) (length = 21) <= remove g1

	// tolerate one GU pair in the seed
	// GU pair considered non-canonical and indicated by <>;
	//   replace {} by () in both targetDotb and guideDotb
	int guPos = guideDotb.indexOf( SecondaryStructure.cgu3 );
	if( guPos > -1 && guPos <= seedLength ) {
	    // get partner
	    int partner = duplex.getPartner( guPos );
	    // change the > by )
	    guideDotb = guideDotb.substring( 0, guPos ) + SecondaryStructure.scanonical3 + guideDotb.substring( guPos + 1, guideDotb.length() );
	    // change the < by (
	    targetDotb = targetDotb.substring( 0, partner ) + SecondaryStructure.scanonical5 + targetDotb.substring( partner + 1,  targetDotb.length() );
	}

	// get the 7 first symbols from each Dotb (g2-g8 in guide)
     	String targetSeed = targetDotb.substring( targetDotb.length() - seedLength, targetDotb.length() ); // ex: .((((((
     	String guideSeed = guideDotb.substring( 0, seedLength ); // ex: )))))))

	int numberOfBPsInTarget = count( targetSeed, '(' );
	int numberOfBPsInSeed = count( guideSeed, ')' );
	// if there are less than 6 bps, there is no seed
	if( numberOfBPsInSeed < 6 ) {
	    //System.out.println( "0-noncanonical" );
	    return( "0-noncanonical" );
	}
	// if all pairs, treat the 7mer-m8 and 8mer types
	if( numberOfBPsInSeed == 7 && numberOfBPsInTarget == 7 )
	    if( A1 ) {
		//System.out.println( "8-mer" );
		return( "8-mer" );
	    }
	    else {
		//System.out.println( "7-mer-m8" );
		return( "7-mer-m8" );
	    }

	if( numberOfBPsInSeed == 7 && numberOfBPsInTarget == 6 ) {
	    // here, the seed is canonical, but there may be a bulge in the target
	    // check special case with 1 bulge in the target (but not in the first 4 positions in the target)
	    int bulgePos = targetSeed.indexOf( SecondaryStructure.ssingle );
	    if( bulgePos > 0 && bulgePos < 3 ) {
		//System.out.println( "1-bulge" );
		return( "1-bulge" );
	    }
	}

	// get the 6 first symbols from each Dotb (g2-g7 in guide)
     	targetSeed = targetDotb.substring( targetDotb.length() - seedLength + 1, targetDotb.length() ); // ex: .(((((
     	guideSeed = guideDotb.substring( 0, seedLength - 1 ); // ex: ))))))

	numberOfBPsInTarget = count( targetSeed, '(' );
	numberOfBPsInSeed = count( guideSeed, ')' );

	// if there are 6 bps, treat the 7mer-A1 and 6mer types
	if( numberOfBPsInTarget == 6 && numberOfBPsInSeed == 6 )
	    if( A1 ) {
		//System.out.println( "7-mer-A1" );
		return( "7-mer-A1" );
	    }
	    else {
		//System.out.println( "6-mer" );
		return( "6-mer" );
	    }

	// check Offset 6mer type
	// get the 6 offset symbols (g3-g8 in guide)
        targetSeed = targetDotb.substring( targetDotb.length() - seedLength, targetDotb.length() - 1 );
     	guideSeed = guideDotb.substring( 1, seedLength );

	numberOfBPsInTarget = count( targetSeed, '(' );
	numberOfBPsInSeed = count( guideSeed, ')' );

	// if there are 6 bps => Offset 6mer type
	if( numberOfBPsInTarget == 6 && numberOfBPsInSeed == 6 ) {
	    //System.out.println( "6-mer-offset" );
	    return( "6-mer-offset" );
	}
	
	//System.out.println( "not supposed to get here!" );
	return( "unknown" );
    }

    public void fold( int e ) {
	// lists to store dotb (states) and their energies
	LinkedList<String> states = new LinkedList<>();
	LinkedList<Double> energies = new LinkedList<>();

	// get target and guide sequences
	String thisTarget = this.getTarget();
	String thisGuide = this.getGuide();

	// ********** MC-FOLD **********
	// *****************************

	// create miRISC mask
	String mask =
	    "'" +
	    "x".repeat( thisTarget.length() - this.kmerBridge - 13 ) +
	    "(((((" +
	    "x".repeat( this.kmerBridge ) +
	    "xxx((((x..x))))xxxxxx)))))xxx" +
	    "x".repeat( thisGuide.length() - 19 ) +
	    "'";

	// fold target including last nt; and guide from g1
	String mcff = "mcff -s " + thisTarget.substring( 0, thisTarget.length() ) +
	    " -sd " +
	    thisGuide.substring( 0, this.getGuideLength() )  +
	    " -m " + mask +
	    " -e " + e +
	    " -v -alt 2> /dev/null";
	// save commandLine for further analyis (maybe)
	this.commandLine = mcff;
	//System.out.println( mcff );
	String[] commands = { "bash", "-c", mcff };
	double mfe = 0.0;
	int mfeState = -1;
	String dotb = "";
	String seed = "";
	try {
	    Process process = Runtime.getRuntime().exec( commands );
	    BufferedReader reader =
		new BufferedReader( new InputStreamReader( process.getInputStream() ) );
	    /*
	      mcff output example:
	      Explored
	      >(null)
	      GCCGGAGAGUUGGUCUCUCCCCUUCUACUGAACAGUGUCUUAGCUGGUUGU
	      mfe(-32.27481), th(-30.81095) (INFO: Duplex mode solutions are not sorted.)
	      <(.<<<<((<(<<....<...(....((((  )))))>.>>)>))>>>>)> -31.268
	      ...
	     */
	    // skip 4 lines;
	    String line;
	    line = reader.readLine();
	    line = reader.readLine();
	    line = reader.readLine();
	    line = reader.readLine();
	    int stateId = 0;
	    // build state and energy lists and save the MFE
	    while( ( line = reader.readLine() ) != null ) {
		// mfe in element[1], minus the "-" sign (3-decimal precision)
		double energy = -Double.parseDouble( line.split( "-", 2 )[1] );
		// in search of the MFE
		if( energy < mfe ) {
		    mfe = energy;
		    mfeState = stateId;
		}
		energies.add( energy );
		states.add( line );
		stateId++;
	    }
	    reader.close();
	} catch (IOException exc ) {
	    exc.printStackTrace();
	}
	// assign the MFE
	this.MFE = mfe;

	// BUG no states found
	if( states.size() == 0 ) {
	    System.out.println( "MFE is : " + mfe + " == " + this.MFE + " for " + states.size() + " dotbs" );
	    System.out.println( this.commandLine );
	    return;
	}

	// ********** SEED ANALYSIS **********
	// ***********************************

	// assign MFE state seed and state
	// assign transient state and seed (if dynamic seed)
	//    a dynamic seed is when the seed type (Bartel) is better
	//    in a transient state than the MFE, static otherwise
	int bestSeedScore = this.seedScore( "0-noncanonical" );
	int bestSeed = -1;
	Double THRESHOLD = 0.0001;
	Duplex MFEState = null;
        Duplex transientState = null;

	// iterate through the conformational states for MFE and possible
	//   transient state with better seed => dynamic duplex type
	for( int i = 0; i < states.size(); i++ ) {
	    Duplex buffer = new Duplex(
				       thisTarget.substring( 0, this.getTargetLength() - 1 ),
				       thisGuide.substring( 1, this.getGuideLength() ),
				       states.get( i ),
				       energies.get( i ) );
	    String seedType = determineSeedType( buffer, this.gettn() == 'A' );
	    //System.out.println( "seedType is " + seedType );
	
	    // case MFE
	    if( i == mfeState ) {
		this.MFESeed = seedType;
		this.MFEState = buffer.getState();
		this.MFEDuplex = buffer;
	    }
	    // case transient state with better seed score
	    else if( seedScore( seedType ) < bestSeedScore ) {
		this.transientDuplex = buffer;
	 	this.transientSeed = seedType;
		this.transientState = buffer.getState();
		this.transientDeltaG = energies.get( i ) - this.MFE;
		bestSeed = i;
	 	bestSeedScore = seedScore( seedType );
	    }
	}
	if( seedScore( this.transientSeed ) < seedScore( this.MFESeed ) )
	    this.dynamic = true;

	// no state assigned BUG
	//System.out.println( "target: " + this.target.getSequence() );
	//System.out.println( "guide : " + this.guide.getSequence() );
	//System.out.println( "command line: " + this.commandLine );

	// ********** OTHER DUPLEX FEATURES **********
	// *******************************************
    
	// Duplex contains the index and basePairs
	//     it starts at g2 as index 0
	//     therefore all indexes in Duplex have an offset of -1
	//     versus MiRISCDuplex

	// bridge
	Duplex duplex = this.getDuplex();
	//if( duplex == null )
	//    System.out.println( "duplex is null and " + this.getDynamic() );

	// determine last base paired nt in the seed
	int lastSeedBp = seedHi; // seedHi = g8 = 6; starts at g8 and go down
	while( duplex.getPartner( lastSeedBp ) == null ) lastSeedBp--;
	Integer lastSeedBpPartner = duplex.getPartner( lastSeedBp ); // partner index in target
	// this.lastSeedBasePair = duplex.getIndex( lastSeedBp ); // assign instance attribute lastSeedBasePair
	// get the target index of the partner of the next paired nt in the guide, jumping the A-box (from g12)
	//   this works since getPartner uses the base pairing index (ie. no bulge considered in the guide)
	int nextBp = suppLo; // first nt in the B-box
	while( duplex.getPartner( nextBp ) == null && nextBp < guideHi ) nextBp++;
	Integer nextPartner = duplex.getPartner( nextBp );
	// bridge length is the gap between the lastSeedPartner and the nextPartner from g12 (suppLo)
	this.bridge = lastSeedBpPartner - nextPartner - 1;
	// if( this.bridge >= 16 ) {
	//     System.out.println( duplex );
	//     System.out.println( "lastSeedBpPartner: " + lastSeedBpPartner + ", nextPartner: " + nextPartner );
	//     System.exit( 0 );
	// }

	// determine the length (number of consecutive canonical (including GU) base pairs in the supplementary region
	//   pairs involving g12-g17 (g12 = 10 = suppLo; suppHi = 15) are considered; guideHi = 19 (g21)
	this.supplementary = 0;
	this.beyondSupplementary = 0;
	this.supplementaryGU = 0;
	this.supplementaryBulges = 0;
	// consecutive bps in the supplementary and D regions (penalize bulges in strand and anti-strand)
	for( int i = suppLo; i <= guideHi; i++ )
	    if( SecondaryStructure.sCanonicalGU3.contains( String.valueOf( duplex.getBasePairType( i ) ) ) ) {
		if( i <= suppHi ) this.supplementary++;
		this.beyondSupplementary++;
		if( SecondaryStructure.sgu3.contains( String.valueOf( duplex.getBasePairType( i ) ) ) )
		    if( i <= suppHi ) this.supplementaryGU++;
	    }

	// compute the number of bulges
	Integer partnerSuppLo = duplex.getPartner( suppLo ); // suppLo = 10
	Integer partnerSuppHi = duplex.getPartner( suppHi ); // suppHi = 15
	Integer partnerDLo = duplex.getPartner( DLo ); // DLo = 16
	Integer partnerDHi = duplex.getPartner( DHi ); // DHi = 19
	int suppBulges = 0;
	int bsupBulges = 0;
	int betwBulges = 0;
	// bulges in the supplementary region
	if( partnerSuppHi != null )
	    suppBulges = ( partnerSuppLo - partnerSuppHi ) - ( suppHi - suppLo );
	// bulges in the D region
	if( partnerDHi != null )
	    bsupBulges = ( partnerDLo - partnerDHi ) - ( DHi - DLo );
	// bulges between the supplementary and D regions
	if( partnerDLo != null )
	    betwBulges = ( partnerSuppHi - partnerDLo ) - ( DLo - suppHi );
	this.supplementaryBulges = suppBulges;
	this.beyondSupplementary -= ( suppBulges + bsupBulges + betwBulges );
	
	// Yan et al. modules
	this.moduleA = 0;
	for( int i = ALo; i <= AHi; i++ )
	    if( duplex.getPartner( i ) != null && duplex.getBasePairType( i ) == SecondaryStructure.ccanonical3 ) this.moduleA++;
	this.moduleB = 0;
	for( int i = BLo; i <= BHi; i++ )
	    if( duplex.getPartner( i ) != null && duplex.getBasePairType( i ) == SecondaryStructure.ccanonical3 ) this.moduleB++;
	this.moduleC = 0;
	for( int i = CLo; i <= CHi; i++ )
	    if( duplex.getPartner( i ) != null && duplex.getBasePairType( i ) == SecondaryStructure.ccanonical3 ) this.moduleC++;
	this.moduleD = 0;
	for( int i = DLo; i <= DHi; i++ )
	    if( duplex.getPartner( i ) != null && duplex.getBasePairType( i ) == SecondaryStructure.ccanonical3 ) this.moduleD++;
    }

    public String getGuide() { return this.guide.getSequence(); }
    public int getGuideLength() { return this.guide.getLength(); }

    public String getTarget() {
	//System.out.println( "sequence: " + this.target.getSequence() + " tStart: " + tStart + " tEnd: " + tEnd );
	return this.target.getSequence().getSequence().substring( tStart, tEnd + 1 );
    }

    public int getTargetLength() { return this.getTarget().length(); }
    public char gettn() { return this.getTarget().charAt( this.getTargetLength() - 1 ); }

    public int getModuleA() { return this.moduleA; }
    public int getModuleB() { return this.moduleB; }
    public int getModuleC() { return this.moduleC; }
    public int getModuleD() { return this.moduleD; }

    public String getSeedType() {
	if( this.dynamic ) return this.transientSeed;
	else return this.MFESeed;
    }

    public double getMFE() {
	if( this.dynamic ) return this.MFE + this.transientDeltaG;
	else return this.MFE;
    }

    public String getState() {
	if( this.dynamic ) return this.transientState;
	else return this.MFEState;
    }

    public int getSeedScore() {	return seedScore( this.getSeedType() ); }
    public double getFreeEnergy() {
	if( this.dynamic ) return( this.transientDuplex.getFreeEnergy() );
	else return this.MFEDuplex.getFreeEnergy();
    }

    public String getDynamic() {
	if( this.dynamic ) return "Dynamic";
	else return "Static";
    }

    public Duplex getDuplex() {
	if( this.dynamic ) return this.transientDuplex;
	else return this.MFEDuplex;
    }

    public String getRegion() {	return this.region; }
    public int getSupplementary() { return this.supplementary; }
    public int getSupplementaryBulges() { return this.supplementaryBulges; }
    public int getSupplementaryGU() { return this.supplementaryGU; }
    public int getBeyondSupplementary() { return this.beyondSupplementary; }
    public int getBridge() { return this.bridge; }

    public static double pairingScore = 10.0d;

    public static double miScore( int A, int B, int C ) {
	double score = 0.0;
	double percentA = A / 3d;
	double percentB = B / 3d;
	double percentC = C / 3d;
	if( percentB > 0 )
	    score += percentB * pairingScore;
	else return score;
	if( percentC > 0 )
	    score += percentC * pairingScore;
	else return score;
	if( percentA > 0 )
	    score += percentA * pairingScore;
	return score;
    }

    public double cheapScore() {
	int numberBps = this.moduleA + this.moduleB + this.moduleC;
	int score = seedScore( this.getSeedType() );
	if( score == 0 ) numberBps += 8;
	else if( score > 1 && score < 3 ) numberBps += 7;
	else if( score > 2 && score < 6 ) numberBps += 6;
	return score * 3.0;
    }
		       
    public double miScore() {
	if( this.getSeedScore() <= 2 ) // perfect seed matching
	    return miScore( this.moduleA, this.moduleB, this.moduleC );
	else return 0;
    }	    	    

    // attributes
    protected Guide guide;
    protected ProteinCodingTranscript target;
    protected int tStart; // t1
    protected int tEnd;   // tn
    //protected int extension; // extension for target folding
    protected String commandLine;
    protected boolean dynamic = false;
    protected String MFESeed = "unknown";
    protected Double MFE = 0.0;
    //protected Double targetMFE = 0.0;
    protected int t8 = -1;
    //protected int seedAccessibility = 0;
    protected String MFEState = null;
    //protected String targetState = null;
    protected Duplex MFEDuplex = null;
    protected String transientSeed = "unknown";
    protected Double transientDeltaG = 0.0;
    protected String transientState = null;
    protected Duplex transientDuplex = null;
    protected String region;
    protected int lastSeedBasePair;
    protected int kmerBridge;
    protected int bridge;
    protected int supplementary;
    protected int beyondSupplementary;
    protected int supplementaryGU;
    protected int supplementaryBulges;
    protected int moduleA;
    protected int moduleB;
    protected int moduleC;
    protected int moduleD;

    /** example:
     //                 [bridge]
     //              1         2
     //   t:1234567890123456789012345678
     //  5'-GGACACCCAGCUAGGACCAUUACUGCCA-3' Sirt1 [tStart-tEnd]
     //       ||| |||||||.||    |||||||
     //  3'---UGUUGGUCGAUUCU----GUGACGGU-5' miR-34a
     //     g:21098765432109    87654321
     //         2         1
     // regions:    [su]        [ seed ]
     //   bridge in target up to 15 nts
     //   supplementary region (su) 12-17
     //   beyond supplementary (bsu) 18-21
     //   seed type: 8-mer (7 bps + A1)
     //   tEnd - tStart = MRE length
     //   tEnd is the position fronting g1,
     //       it determines the A1 feature of Bartel 2009.
     //
     **/


    // create an instance of MiRISC duplex with miR and RNA, MRE defined by
    //     tStart and tEnd.
    public MiRISCDuplex( Guide guide, ProteinCodingTranscript target, int tStart, int tEnd, int bridgeLen ) {
	this.guide = guide;
	this.target = target;
	this.tStart = tStart;
	this.tEnd = tEnd;
	this.kmerBridge = bridgeLen;
	//this.extension = extension;
	
	// fold init -e value is 4
	int initE = 4;
	// fold until dot-brackets are generated
	while( this.MFE == 0.0 ) fold( initE++ );
    }

    // assume left dotb of a mcff duplex mode folding and left side longer than the right side
    //    ie. MRE longer than the miR
    // ex) ...<..(((<(<(<(<{({(.....(.(((((  )))))))})}>)>)>.)>)))>
    public static String barRep( String target, String guide, String dotbLeft, String dotbRight, String space, String newline ) {
	String alignedTarget = "";
	String bars = "";
	String alignedGuide = "";

	/*
	  Left+Right: ......(..(....(<((((((((<(((((((  )))))))>))..))))))>)))
	  target: AUUCUAUUACAAUGAAUUUUGCAGUUUUGCAC
	  guide: GUGCAAAUCUAUGCAAAACUGA
	  result:
          AUUCUAUUACAAUGAAUUUUGC--AGUUUUGCAC
                |  |    |!||||||  ||!|||||||   
          ------A--G----UCAAAACGUAUCUAAACGUG
	 */

	// i index in left side (target)
	// j index in right side (guide)
	int j = dotbRight.length() - 1;
	int i = 0;
	// these positions in lowercase
	Integer lowerCasePos[] = { 8, 9, 10, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
	Set<Integer> lowerCasePositions = new HashSet<>( Arrays.asList( lowerCasePos ) );
	int g = guide.length();
	while( i < target.length() ) {
	    if( dotbLeft.charAt( i ) == '.' ) {
		alignedTarget += Character.toLowerCase( target.charAt( i++ ) );
		bars += space;
		alignedGuide += "-";
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
		    bars += "|";
		    if( lowerCasePositions.contains( g ) )
			alignedGuide += Character.toLowerCase( guide.charAt( j-- ) );
		    else alignedGuide += guide.charAt( j-- );
		    g--;
		}
		else {
		    String pair = "" + target.charAt( i ) + guide.charAt( j );
		    if( StringSequence.guBps.contains( pair ) )
			bars += "!";
		    else bars += ":";
		    alignedTarget += Character.toLowerCase( target.charAt( i++ ) );
		    if( lowerCasePositions.contains( g ) )
			alignedGuide += Character.toLowerCase( guide.charAt( j-- ) );
		    else alignedGuide += guide.charAt( j-- );
		    g--;
		}
	    }
	}
	return alignedTarget + newline + bars + newline + alignedGuide;
    }

    // reverse String
    public static String reverse( String str ) {
	String rev = "";
	for( int i = str.length() - 1; i >= 0; i-- )
	    rev = rev + str.charAt(i);
	return rev;
    }

    public String dataPrint() {
	String output = "";
	if( this.dynamic )
	    output += "dynamic, " + this.getSeedType() + ", " + (this.MFE + this.transientDeltaG) + ", ";
	else
	    output += "static, " + this.getSeedType() + ", " + this.MFE + ", ";
	output += "A: " + this.moduleA + ", B: " + this.moduleB + ", C: " + this.moduleC + ", D: " + this.moduleD + ", miscore: " + String.format("%.2f", this.miScore() ) + ", bridge: " + this.bridge + ", suppBps: " + this.supplementary + ", extBps: " + this.getBeyondSupplementary();
	output += ", suppGU: " + this.supplementaryGU + ", bulgeNtsSupp: " + this.supplementaryBulges + ", bulgesTarget: " + this.getDuplex().getNumberOfBulgesInStrand() + ", bulgesGuide: " + this.getDuplex().getNumberOfBulgesInAntiStrand() + "\n";
	//output += this.target.getKd() + "\n"; // + "\t" + this.target.getKcat() + "\n";
	return output;
    }

    /*** example
5'-UCUCAAUGACAACUUUGUCAAGCUCAUUUCC-3' [991-1021] null NM_001289726.1 overlapping
          :||:|||::|   |::|||||||     Static 7-mer-m8 deltaG:-39.687 A:1 B:1 C:3 D:2 miscore:16.67 bridge:3 suppBps:4 extBps:6 suppGU:1 bulgeNtsSupp:0 bulgesTarget:1 bulgesGuide:0
3'--------UUGGUGAUGC---UUUAGUAAAGG-5' shNM_139305.2:1425
     */

    public String numbersPrint( String newline, String space ) {
	/*
          Transform:
	  miR-19a-3p 5'-UGUGCAAAUCUAUGCAAAACUGA-3'
	  PTEN 3252-5'-AUUCUAUUACAAUGAAUUUUGCAGUUUUGCAC-3'-3283
	  type: static 8-mer
	  A1: true
	  ground state: ......(..(....(<((((((((<(((((((  )))))))>))..))))))>)))
	  seed: 8-mer
	  MFE: -46.011
	  into:
	  "5'-cagcaagaauaaaaaagaaauccauaucuuaaagaA-3'
                            ||| ||||||:||:|||||||
	   3'---------------uuc-uUAGGUUuacAAUUUCUU-5'"
	 */
	String output = "";
	String guide = this.getGuide().substring( 1, this.getGuideLength() );
	String target = this.getTarget().substring( 0, this.getTargetLength() - 1 );
	String state = this.getState();
	String dotbLeft = state.split( " ", 3 )[0]; // left dotb
	String dotbRight = state.split( " ", 3 )[2]; // right dotb
	String bars = barRep( target, guide, dotbLeft, dotbRight, space, newline );
	// data
        //output += this.getDynamic() + space + this.getSeedType() + "\t" + this.getMFE() + "\t" ;
        //output += this.bridge + "\t" + this.supplementary + "\t";
	//output += this.getBeyondSupplementary() + "\t" + this.supplementaryGU + "\t";
	//output += this.supplementaryBulges + "\t" + this.getDuplex().getNumberOfBulgesInStrand() + "\t" + this.getDuplex().getNumberOfBulgesInAntiStrand() + "\t";
	// bars has 3 components separated by '\n': target + bars + guide
	String[] bar = bars.split( newline );
	output +=
	    "5'-" +
	    bar[0] +
	    this.gettn() +
	    "-3'" + space +
	    this.target.getName() + space + "[" + this.tStart + "-" + this.tEnd + "]" + space + this.target.stringRegion( this.tStart, this.tEnd ) +
	    newline;
	output +=
	    space + space + space +
	    bar[1] +
	    space + space + space + space + space +
	    this.getDynamic() + space +
	    this.getSeedType() + space +
	    "deltaG:" + this.getMFE() + " A:" + this.moduleA + " B:" + this.moduleB + " C:" + this.moduleC + " D:" + this.moduleD + " miscore:" + String.format("%.2f", this.miScore() ) + newline;
	output +=
	    "3'-" +
	    bar[2] +
	    this.guide.getSequence().charAt( 0 ) +
	    "-5'" + space +
            "guide" + " bridge:" + this.bridge + " suppBps:" + this.supplementary + " extBps:" + this.getBeyondSupplementary() +
	    " suppGU:" + this.supplementaryGU + " bulgeNtsSupp:" + this.supplementaryBulges + " bulgesTarget:" + this.getDuplex().getNumberOfBulgesInStrand() + " bulgesGuide:" + this.getDuplex().getNumberOfBulgesInAntiStrand();
	return output;
    }
    
    public String prettyPrint( String newline, String space ) {
	/*
          Transform:
	  miR-19a-3p 5'-UGUGCAAAUCUAUGCAAAACUGA-3'
	  PTEN 3252-5'-AUUCUAUUACAAUGAAUUUUGCAGUUUUGCAC-3'-3283
	  type: static 8-mer
	  A1: true
	  ground state: ......(..(....(<((((((((<(((((((  )))))))>))..))))))>)))
	  seed: 8-mer
	  MFE: -46.011
	  into:
	  5'-AUUCUAUUACAAUGAAUUUUGCAGUUUUGCACA-3' PTEN[3252-3283]
              |!|!|!!|| |      ||!|!  |||||||     static 8-mer (-46.011 kcal/mol)
          3'--AGUCAAAAC-G------UAUCU--AAACGUGU-5' hsa-miR-19a-3p
	 */
	String output = "";
	String guide = this.getGuide().substring( 1, this.getGuideLength() );
	String target = this.getTarget().substring( 0, this.getTargetLength() - 1 );
	String state = this.getState();
	String dotbLeft = state.split( " ", 3 )[0]; // left dotb
	String dotbRight = state.split( " ", 3 )[2]; // right dotb
	String bars = barRep( target, guide, dotbLeft, dotbRight, space, newline );
	// bars has 3 components separated by '\n': target + bars + guide
	String[] bar = bars.split( newline );
	//output += this.commandLine + "\n";
	output += "5'-" + bar[0] + this.gettn() + "-3'" + space + this.target.getName() + space + "[" + this.tStart + "-" + this.tEnd + "]" + space + this.target.stringRegion( this.tStart, this.tEnd ) + newline;
	output += space + space + space + bar[1] + space + space + space + space + space + this.getDynamic() + space + this.getSeedType() + space + "(" + this.getMFE() + space + "kcal/mol)" + space;
	//	output += "A:" + this.moduleA + space + "B:" + this.moduleB + space + "C:" + this.moduleC + space + "D:" + this.moduleD + space + "bridge:" + this.bridge + space + "supp:" + this.supplementary + space + "suppGU:" + this.supplementaryGU + space + "suppBulge:" + this.supplementaryBulges + space + "bsup:" + this.getBeyondSupplementary() + space + "btarget:" + this.getDuplex().getNumberOfBulgesInStrand() + space + "bguide:" + this.getDuplex().getNumberOfBulgesInAntiStrand() + space + newline;
	output += "bridge:" + this.bridge + space + "supp-bps:" + this.supplementary + space + "supp-GU:" + this.supplementaryGU + space + "bulges-supp:" + this.supplementaryBulges + space + "ext-supp-bps:" + this.getBeyondSupplementary() + space + "bulge-target:" + this.getDuplex().getNumberOfBulgesInStrand() + space + "bulges-guide:" + this.getDuplex().getNumberOfBulgesInAntiStrand() + space + newline;
	output += "3'-" + bar[2] + this.guide.getSequence().charAt( 0 ) + "-5'" + newline;
	//output += this.getDuplex() + "\n";
	//output += "lastSeedBasePair: " + this.lastSeedBasePair + "\n";
	//output += "bridge: " + this.bridge + "\n";
	//output += "supplementary: " + this.supplementary + "\n";
	return output;
    }

    // Override
    public String toString() {
	//return prettyPrint( "<br>", "&nbsp;" );
	return numbersPrint( "\n", " " );
	//return dataPrint();
    }
}
