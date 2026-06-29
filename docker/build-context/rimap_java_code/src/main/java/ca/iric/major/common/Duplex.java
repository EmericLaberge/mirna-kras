/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
* Duplex is a class for generic RNA duplexes formed by a strand and antistrand
*
* @license     MIT
* @author      Francois Major, Université de Montréal
* @version     %I%, %G%
* @since       1.0
*/

import java.lang.StringBuffer;

import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Arrays;
import java.util.Collections;

import java.lang.Runtime;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.io.File;
import java.io.Serializable;


/** -------------------------------------------
      Duplex
        is a class to handle RNA duplex conformational states, which
        has the following attributes:
        Double freeEnergy; // DeltaG energy
        String state = ""; // dot-bracket of the duplex, initially empty
        String strand; // sequence of the strand
	String antiStrand; // sequence of the anti strand
	Map<Integer,Integer> index; // index for the strand
	Map<Integer,Integer> basePairs; // base pairs
*/

public class Duplex implements Serializable {

    private static final long serialVersionUID = 1L;

    protected static final String LOOP =     "AAAAGGAAACAAAA";
    protected static final String LOOPMASK = "....((..))....";
    // private ProcessBuilder processBuilder = new ProcessBuilder();
    private static String os = "";
    private static String mcftablesPath = "";
    private static String mcfPath = "";

    private void setOS() {
	// ***** READ OS and use appropriate file system separator
	Duplex.os = System.getProperty( "os.name" ).toLowerCase();
	    // setup the processBuilder
	    // Map<String,String> environment = this.processBuilder.environment();
	    // Duplex.mcftablesPath = "C:" + File.separator + "Users" + File.separator + "shayy" + File.separator + "OneDrive" + File.separator + "Bureau" + File.separator + "work" + File.separator + "tables" + File.separator;
	    // Duplex.mcfPath = "C:" + File.separator + "Users" + File.separator + "shayy" + File.separator + "OneDrive" + File.separator + "Bureau" + File.separator + "work" + File.separator;
	    // environment.put( "MCFTABLES", Duplex.mcftablesPath );
	    //environment.put( "PATH", Duplex.mcfPath );
	    //}
	    //else if( Duplex.os.contains( "mac" ) || Duplex.os.contains( "nix" ) || Duplex.os.contains( "nux" ) || Duplex.os.contains( "aix" ) ) {
	    //} else Util.stop( "Duplex does not recognize the OS", 0 );
    }

    protected Double freeEnergy;
    protected String state = "";
    protected boolean maximizeDoubleHelix = false;
    protected boolean keepAll = false;
    protected String strand;
    protected String strandState;
    protected String antiStrand;
    protected String antiStrandState;
    protected String mask;
    protected boolean unbalancedMask;
    protected int e;
    protected Map<Integer,Integer> index;
    protected Map<Integer,Integer> basePairs;
    protected Map<Integer,Character> basePairTypes;
    protected List<String> states;
    protected List<Double> energies;
    protected List<Integer> numberOfBPS;
    protected List<Duplex> duplexes = null;

    // getters

    public Double  getFreeEnergy()                 { return this.freeEnergy; }
    public String  getState()                      { return this.state; }
    public String  getStrandState()                { return strandState; }
    public String  getAntiStrandState()            { return antiStrandState; }
    public String  getStrand()                     { return this.strand; }
    public char    getStrand( int i )              { return this.strand.charAt( i ); }
    public String  strandSubstring( int i, int j ) { return this.strand.substring( i, j ); }
    public String  getAntiStrand()                 { return this.antiStrand; }
    public Integer getIndex( int i )               { return this.index.get( i ); }
    public String  getMask()                       { return this.mask; }
    public Integer getE()                          { return this.e; }
    public Map<Integer,Integer> getIndex()         { return this.index; }
    public Map<Integer,Integer> getBasePairs()     { return this.basePairs; }
    public List<String> getStates()                { return this.states; }
    public List<Double> getEnergies()              { return this.energies; }
    public List<Duplex> getConformationalSpace()   { return this.duplexes; }

    // buildIndex creates the bps map (this.basePairs), antistrand index (this.index),
    //   replaces the <> by {} for the GU bps and remove energy at the end of this.state
    //
    //  assumes mcff dotb line with energy at the end in this.state
    //  side effect: modifies the dotb using ( ) canonical (AU, CG)
    //                                       { } GU
    //                                       < > non-canonical (all others)
    //  side effect: modifies this.state
    //  side effect: assigns this.strandState and this.antistrandState
    //  strand, antistrand, and state example:
    //
    //          [    strand   ]LOOP[ antistrand ]
    //  mcff -s AUAAAACACCCAGCUGAAAUCUUAGCUGGUUGU -m 'ppppppppppppppp....qqqqqqqqqqqqqq' -t 2
    //	        .....((((((((((........)))))))))) -24.384 ()
    //
    public void buildIndex() {
	String strandDotb = this.state.substring( 0, this.strand.length() );
	String antiStrandDotb = this.state.substring( this.strand.length() + LOOPMASK.length(), this.strand.length() + LOOPMASK.length() + this.antiStrand.length() );
	String initialDotb = strandDotb + "  " + antiStrandDotb;

	// initialize the strandState
	this.strandState = "";
	// for( int k = strandDotb.length() - 1; k >= 0; k-- ) // insert the terminal dots
	//     if( strandDotb.charAt( k ) == '.' ) this.strandState += Utilities.ssingle;
	//     else break;
	// System.out.println( "init value for this.strandState: " + this.strandState + "; added " + this.strandState.length() + " dots" );

	// initialize the antiStrandState
	this.antiStrandState = "";
	Stack<Integer> s = new Stack<>(); // list of parentheses
	this.index = new HashMap<>(); // antiStrand index { antiStrand_i:strand_i }
	this.basePairs = new HashMap<>(); // result
	this.basePairTypes = new HashMap<>(); // result; bp types: ), >, or }

	// if( strandDotb.length() != this.strand.length() ) {
	//     System.out.println( "state : " + this.state );
	//     System.out.println( "strand: " + this.strand );
	//     System.out.println( "stndDB: " + strandDotb );
	//     System.out.println( "antiSd: " + this.antiStrand );
	//     System.out.println( "antiDB: " + antiStrandDotb );
	//     System.out.println( "mask  : " + this.mask );
	//     System.out.println( "state : " + this.state );
	//     Util.stop( "strand and strand dotb don't match: " + strandDotb + "[" + strandDotb.length() + "] strand: " + this.strand + "[" + this.strand.length() + "]", 0 );
	// }
	
	// iterate on each symbol in dotb, keep its index (j)
	// suck the strand 5' opening bps
	int j = 0;
	for( int indexC = 0; indexC < strandDotb.length(); indexC++ ) {
	    char c = strandDotb.charAt( indexC );
	    // push 5' bp parenthesis
	    if( SecondaryStructure.parent5.indexOf( c ) != -1 ) s.push( j );
	    else
		// only the '.' is allowed
		if( c != SecondaryStructure.csingle ) Utils.stop( "Illegal symbol in dot bracket [" + strandDotb + "]", 12 );
	    j += 1;
	}
	// suck the antiStrand 3' closing bps; build the antiStrand index,
	// and the new dot brackets: <> replaced by {} for GUs
	int previousI = j; int g = 0; int i = strandDotb.length();
	j = 0;
	for( int indexC = 0; indexC < antiStrandDotb.length(); indexC++ ) {
	    char c = antiStrandDotb.charAt( indexC );
	    if( SecondaryStructure.parent3.indexOf( c ) != -1 ) {
		// if no opening parenthesis => error
		if( s.isEmpty() ) Utils.stop( "Unbalanced dot bracket [" + initialDotb + "]", 13 );
		else {
		    // get 5' partner and store the bp; c is the symbol in antiStrand
		    // j is the index of nt in antiStrand; i = s.pop() the index of the partner in strand
		    // antiStrand.charAt(j) is the base in the antiStrand
		    // strand.charAt(i) is the base of the strand partner; strandDotb.charAt(i) is the strand symbol
		    i = s.pop();
		    // deal with the single stranded symbol
		    for( int k = 0; k < previousI - i - 1; k++ )
			this.strandState += SecondaryStructure.ssingle; // ssingle is the String "."
		    this.basePairs.put( j, i ); // put 5':3' pair
		    //this.basePairs.put( i, j ); // put 3':5' pair
		    this.index.put( g, j );
		    // if GUPair, change symbols in dotb
		    if( StringSequence.guBps.contains( "" + this.strand.charAt( i ) + this.antiStrand.charAt( j ) ) ) {
			this.strandState += SecondaryStructure.sgu5;
			this.antiStrandState += SecondaryStructure.sgu3;
			this.basePairTypes.put( j, SecondaryStructure.cgu3 );
		    }
		    else {
			this.strandState += strandDotb.charAt( i );
			this.antiStrandState += c;
			this.basePairTypes.put( j, c );
		    }
		}
	    }
	    else {
		if( c != SecondaryStructure.csingle ) Utils.stop( "Illegal symbol in dot bracket [" + initialDotb + "]", 12 );
		this.antiStrandState += SecondaryStructure.ssingle;
		g -= 1;
	    }
	    previousI = i;
	    g += 1;
	    j += 1;
	}
	// add the loop symbols at the beginning
	for( int k = 0; k < previousI; k++ ) this.strandState += SecondaryStructure.ssingle;
	// reverse the strandState
	this.strandState = new StringBuffer( this.strandState ).reverse().toString();
	// if stack is empty, it balances => correct 2D structure
	if( s.isEmpty() ) {
	    if( initialDotb.length() != ( this.strandState.length() + 2 + this.antiStrandState.length() ) )
		Utils.stop( "Something went wrong with the building of the new dotb " + initialDotb + " " + this.strandState + "  " + this.antiStrandState, 14 );
	}
	else Utils.stop( "Unbalanced dot bracket [" + initialDotb + "]", 13 );
	this.state = this.strandState + "  " + this.antiStrandState;
    }

    // i is the nt in the antistrand/guide
    public Integer getPartner( int i ) {
    	return this.basePairs.get( i );
    }

    public Integer getBasePairPartner( int i ) { return this.basePairs.get( i ); }

    public Character getBasePairType( int i ) {
	return this.basePairTypes.get( i );
    }

    public boolean isShifted() { // when the first base pair is diagonal. ex) .((((((.  ).))))), -18.012
	String strandState = this.getStrandState();
	String antiStrandState = this.getAntiStrandState();
	int asFirst = 0;
	int sLast = strandState.length() - 1;
	int asFirstBP = -1;
	int sFirstBP = -1;
	for( int i = asFirst; i < antiStrandState.length(); i++ )
	    if( antiStrandState.charAt( i ) == ')' ) { asFirstBP = i; break; }
	for( int i = sLast; i >= 0; i-- )
	    if( strandState.charAt( i ) == '(' ) { sFirstBP = i; break; }
	if( sLast - sFirstBP > asFirstBP ) return true;
	return false;
    }

    public boolean isBulgy() { // when there is a >1-nt bulge in the strand
	return this.getStrandState().contains( "(.." );
    }

    public boolean isCanonical( int i ) {
	Character bp = '?';
	try {
	    bp = this.basePairTypes.get( i );
	} catch( NullPointerException e ) {
	    System.out.println( "isCanonical: basePairTypes map is null, i = " + i );
	    System.exit( 1 );
	}
	return bp == null ? false : bp == SecondaryStructure.ccanonical3;
    }

    public boolean isNonCanonical( int i ) {
	Character bp = this.basePairTypes.get( i );
	return bp == null ? false : bp == SecondaryStructure.cnonCanonical3;
    }

    public boolean isGU( int i ) {
	Character bp = this.basePairTypes.get( i );
	return bp == null ? false : bp == SecondaryStructure.cgu3;
    }

    public boolean isCanonicalGU( int i ) {
	//System.out.print( "isCanonicalGU( " + i + " ) = " );
	Character bp = this.basePairTypes.get( i );
	//System.out.println( ""+bp );
	return bp == null ? false : bp == SecondaryStructure.ccanonical3 || bp == SecondaryStructure.cgu3;
    }

    public boolean isPaired( int i ) {
	//Utils.debug( "isPaired( " + i + " ) = " );
	Character bp = this.basePairTypes.get( i );
	//Utils.debug( "" + bp );
	return bp == null ? false : true;
    }

    public int getNumberOfCanonical( int start, int end ) {
	int count = 0;
	for( int i = start; i <= end; i++ )
	    if( isCanonical( i ) ) count++;
	return count;
    }

    public boolean isPaired( int start, int end ) {
	Integer startPartner = this.getBasePairPartner( start );
	Integer endPartner = this.getBasePairPartner( end );
	//System.out.println( start+":"+startPartner + " " + end+":"+endPartner );
	if( startPartner == null || endPartner == null ) return false;
	if( startPartner - endPartner == end - start )
	    for( int i = start; i <= end; i++ ) {
		if( !isPaired( i ) ) return false;
	    }
	else return false;
	return true;
    }

    public boolean isCanonicallyPaired( int start, int end ) {
	Integer startPartner = this.getBasePairPartner( start );
	Integer endPartner = this.getBasePairPartner( end );
	//System.out.println( start+":"+startPartner + " " + end+":"+endPartner );
	if( startPartner == null || endPartner == null ) return false;
	if( startPartner - endPartner == end - start )
	    for( int i = start; i <= end; i++ ) {
		if( !isCanonical( i ) ) return false;
	    }
	else return false;
	return true;
    }

    public boolean isCanonicallyPairedGU( int start, int end ) {
	Integer startPartner = this.getBasePairPartner( start );
	Integer endPartner = this.getBasePairPartner( end );
	//System.out.println( start+":"+startPartner + " " + end+":"+endPartner );
	if( startPartner == null || endPartner == null ) return false;
	if( startPartner - endPartner == end - start )
	    for( int i = start; i <= end; i++ ) {
		if( !isCanonicalGU( i ) ) return false;
	    }
	else return false;
	return true;
    }

    public boolean isCanonicallyPairedWithMax1GU( int start, int end ) {
	Integer startPartner = this.getBasePairPartner( start );
	Integer endPartner = this.getBasePairPartner( end );
	int countGU = 0;
	//System.out.println( start+":"+startPartner + " " + end+":"+endPartner );
	if( startPartner == null || endPartner == null ) return false;
	if( startPartner - endPartner == end - start )
	    for( int i = start; i <= end; i++ ) {
		if( !isCanonicalGU( i ) ) return false;
		if( isGU( i ) ) countGU++;
	    }
	else return false;
	if( countGU > 1 ) return false;
	return true;
    }

    // return the index of last g before first bulge of 1 nt if another or longer bulge, gy otherwise
    public int adjustGLast( int gx, int gy ) {
	int gLast = gx;
	boolean firstBulge = true;
	for( int i = gx; i < gy; i++ ) {
	    String res = this.bulgeOfAdjacentPaired( gLast, gLast + 1 );
	    if( res.length() == 0 ) gLast++;
	    else
		if( res.length() == 1 && firstBulge ) {
		    gLast++;
		    firstBulge = false;
		}
		else break;
	}
	for( int i = gLast; i >= gx; i--)
	    if( this.isCanonicalGU( i ) ) return i; 
	return gy;
    }

    // assume gx == gy - 1
    public String bulgeOfAdjacentPaired( int gx, int gy ) {
	if( gy - gx > 1 ) return ""; // they are not adjacent in the guide
	Integer gxPartner = this.getBasePairPartner( gx );
	Integer gyPartner = this.getBasePairPartner( gy );
	if( gxPartner == null || gyPartner == null ) return ""; // one of the two is not paired
	if( gxPartner - gyPartner >= 2 )
	    return this.strand.substring( gyPartner+1, gxPartner );
	return "";
    }

    public boolean isCanonicallyPairedWithMax1GU1Bulge( int start, int end ) {
	Integer startPartner = this.getBasePairPartner( start );
	Integer endPartner = this.getBasePairPartner( end );
	int countGU = 0;
	//System.out.println( start+":"+startPartner + " " + end+":"+endPartner );
	if( startPartner == null || endPartner == null ) return false;
	if( startPartner - endPartner - 1 == end - start )
	    for( int i = start; i <= end; i++ ) {
		if( !isCanonicalGU( i ) ) return false;
		if( isGU( i ) ) countGU++;
	    }
	else return false;
	if( countGU > 1 ) return false;
	return true;
    }

    public int getNumberOfBulgesInStrand() {
	int n = 0; // initially 0 bulge
	// find first paired nt in the duplex
	int firstPairedNt = 0;
	for( int i = 0; i < this.basePairs.size(); i++ )
	    if( this.basePairs.containsKey( i ) ) break;
	    else firstPairedNt = i + 1;
	int bpi = firstPairedNt;
	int bpj = this.basePairs.get( bpi );
	for( int i = firstPairedNt + 1; i < this.antiStrand.length(); i++ )
	    if( this.basePairs.containsKey( i ) ) {
		int j = this.basePairs.get( i );
		if( Math.abs( j - bpj ) > 1 ) n++;
		bpj = j;
		bpi = i;
	    }
	return n;
    }

    public int getNumberOfBulgesInAntiStrand() {
	int n = 0; // initially 0 bulge
	// find first paired nt in the duplex
	int firstPairedNt = 0;
	for( int i = 0; i < this.basePairs.size(); i++ )
	    if( this.basePairs.containsKey( i ) ) break;
	    else firstPairedNt = i + 1;
	int bpi = firstPairedNt;
	for( int i = firstPairedNt + 1; i < this.antiStrand.length(); i++ )
	    if( this.basePairs.containsKey( i ) ) {
		if( Math.abs( i - bpi ) > 1 ) n++;
		bpi = i;
	    }
	return n;
    }

    public boolean isMismatch( int i ) {
	Integer partner = this.getPartner( i );
	if( partner != null )
	    return this.basePairTypes.get( i ) == SecondaryStructure.cnonCanonical3;
	// no partner
	Integer prevBp = this.getPartner( i-1 );
	Integer succBp = this.getPartner( i+1 );
	if( prevBp != null && succBp != null ) // both have partners
	    return Math.abs( prevBp - succBp ) == 2; // distance of 2 indicates the middle matches i and is not paired
	return false;
    }

    public int bpDistance( int first, int second ) { // compute strand length between first and second bp
	int firstPartner = -1;
	int secondPartner = -1;
	for( int i = first; i > 0; i-- ) {
	    Integer partner = this.getPartner( i );
	    if( partner != null && isPaired( i ) ) {
		firstPartner = partner;
		break;
	    }
	}
	//System.out.println( firstPartner );
	for( int i = second; i < 17; i++ ) {
	    Integer partner = this.getPartner( i );
	    if( partner != null && isPaired( i ) ) {
		secondPartner = partner;
		break;
	    }
	}
	//System.out.println( secondPartner );
	return firstPartner - secondPartner - 1;
    }

    // state assumed to come from an mcff state
    //   <dotb>
    public Duplex( String strand, String antiStrand, String state, Double deltaG ) {
	this.setOS();
	this.strand = strand;
	this.antiStrand = antiStrand;
	this.state = state;
	this.freeEnergy = deltaG;
	buildIndex();
    }

    // this version selects the stats with longest stretch of canonical base pairs
    public Duplex( String strand, String antiStrand, String mask, int e, boolean maximizeDoubleHelix, boolean unbalancedMask ) {
	this.keepAll = true;
	this.setOS();
	this.strand = strand;
	this.antiStrand = antiStrand;
	if( Duplex.os.contains( "win" ) ) this.mask = "\"" + mask + "\"";
	else this.mask = "'" + mask + "'";
	this.e = e;
	this.maximizeDoubleHelix = maximizeDoubleHelix;
	String mcff = "";
	this.unbalancedMask = unbalancedMask;
	String mcffMaskOption = " -m ";
	if( this.unbalancedMask ) mcffMaskOption = " -um ";
	if( this.fold() == 0 ) { // set this.state and this.freeEnergy; return 0 if okay, -1 otherwise
	    this.buildIndex();
	    //System.out.println( "bps: " + this.basePairs );
	}
	else {
	    if( !this.mask.equals( "" ) ) {
		mcff = "mcff -s " + this.strand + Duplex.LOOP + this.antiStrand +
		mcffMaskOption + this.mask + " -t " + this.e;
	    }
	    else
		mcff = "mcff -s " + this.strand + Duplex.LOOP + this.antiStrand +
		" -t " + this.e;
	    //Utils.debug( mcff );
	    Utils.stop( "fold() is not okay", 15 );
	}
    }

    // this version keeps the list of states (optimal and suboptimals)
    public Duplex( String strand, String antiStrand, String mask, int e, boolean unbalancedMask ) {
	this.keepAll = true;
	this.setOS();
	this.strand = strand;
	this.antiStrand = antiStrand;
	if( Duplex.os.contains( "win" ) ) this.mask = "\"" + mask + "\"";
	else this.mask = "'" + mask + "'";
	this.e = e;
	this.maximizeDoubleHelix = maximizeDoubleHelix;
	String mcff = "";
	this.unbalancedMask = unbalancedMask;
	String mcffMaskOption = " -m ";
	if( this.unbalancedMask ) mcffMaskOption = " -um ";
	this.states = new LinkedList<>();
	this.energies = new LinkedList<>();
	this.numberOfBPS = new ArrayList<>();
	if( this.fold() == 0 ) { // set this.state and this.freeEnergy; return 0 if okay, -1 otherwise
	    this.buildIndex();
	    //System.out.println( "bps: " + this.basePairs );
	}
	else {
	    if( !this.mask.equals( "" ) ) {
		mcff = "mcff -s " + this.strand + Duplex.LOOP + this.antiStrand +
		mcffMaskOption + this.mask + " -t " + this.e;
	    }
	    else
		mcff = "mcff -s " + this.strand + Duplex.LOOP + this.antiStrand +
		" -t " + this.e;
	    //Utils.debug( mcff );
	    Utils.stop( "fold() is not okay", 15 );
	}
    }

    private static int countConsecutiveLeftParens( String s ) {
        int mid = s.length() / 2;
        int count = 0;

        for( int i = 0; i < mid; i++ ) {
            if( s.charAt(i) == '(' ) {
                count++;
            } else {
                count = 0; // reset count on interruption
            }
        }
        return count;
    }

    private static int countLeftParens( String s ) {
        int mid = s.length() / 2;
        int count = 0;

        for( int i = 0; i < mid; i++ )
            if( s.charAt(i) == '(' )
                count++;
        return count;
    }

    private void createDuplexes() {
	this.duplexes = new ArrayList<>();
	for( int i = 0; i < this.states.size(); i++ )
	    this.duplexes.add( new Duplex( this.strand, this.antiStrand, this.states.get( i ), this.energies.get( i ) ) );
    }

    // Forms the strand:antistrand duplex using mcff
    public int fold() {
	// lists to store dotb (states) and their energies
	List<String> states = new ArrayList<>();
	List<Double> energies = new ArrayList<>();
	List<Integer> numberOfBPS = new ArrayList<>();

	// ********** MC-FOLD **********
	// *****************************

	String mcff = "";
	String commandLine = "";
	int theEvalue = this.e;
	double mfe = 0.0;
	String mask = "";
	String mcffMaskOption = " -m ";
	if( this.unbalancedMask ) mcffMaskOption = " -um ";
	this.freeEnergy = 0.0;
	while( this.freeEnergy == 0.0 && theEvalue < 11 ) {
	    if( !this.mask.equals( "" ) ) {
		mcff = "mcff -s " + this.strand + Duplex.LOOP + this.antiStrand +
		    mcffMaskOption + this.mask +
		    " -t " + theEvalue + " -alt";
	    }
	    else {
	    mcff = "mcff -s " + this.strand + Duplex.LOOP + this.antiStrand +
		" -t " + theEvalue + " -alt";
	    }
	    // save commandLine for further analyis (maybe)
	    commandLine = mcff;

	    String[] commands = new String[3];
	    if( Duplex.os.contains( "win" ) ) {
		commands[0] = "cmd.exe";
		commands[1] = "/c";
		commands[2] = mcff;
	    }
	    else {
		commands[0] = "bash";
		commands[1] = "-c";
		commands[2] = mcff;
	    }

	    try {
		Process process = Runtime.getRuntime().exec( commands );

		// Create a thread-safe list to store the output lines
		List<String> outputLines = Collections.synchronizedList( new ArrayList<>() );

		Thread inputThread = new Thread(() -> {
			try( BufferedReader reader = new BufferedReader( new InputStreamReader( process.getInputStream() ) ) ) {
			    String line;
			    while( ( line = reader.readLine() ) != null ) {
				outputLines.add( line ); // Add each line to the shared list
			    }
			} catch( IOException e ) {
			    e.printStackTrace();
			}
		});
		
		inputThread.start();
		inputThread.join(); // Wait for the thread to finish

		// Consume the error stream
		Thread errorThread = new Thread( () -> {
			try( BufferedReader error = new BufferedReader( new InputStreamReader( process.getErrorStream() ) ) ) {
			    error.lines().forEach( System.err::println );
			} catch (IOException e) {
			    e.printStackTrace();
			}
		});

		errorThread.start();
		errorThread.join();

		// Wait for the process to complete
		int exitCode = process.waitFor();
		process.destroy();
		
		/*
		  mcff output example:
		  ....(<(((((....)))))>) -19.773 () 
		  ..(..(<((((....))))>)) -16.308 () 
		  ....<(<((((....))))>)> -15.482 () 
		  .<<(...((((....)))))>> -15.290 () 
		  ...
		*/

		//Utils.debug( mcff );
		int stateId = 0;
		// build state and energy lists and save the MFE; // skip 4 lines;
		for( int lineNo = 0; lineNo < outputLines.size(); lineNo++ ) {
		    //line = outputLines.get( lineNo );
		    // mfe in element[1], minus the "-" sign (3-decimal precision)
		    //System.out.println( outputLines.get( lineNo ) );
		    Double energy;
		    String[] parts = outputLines.get( lineNo ).split( "[+-]", 2 );
		    String forEnergy = parts[1].split( " ", 2 )[0]; // get rid of the shape

		    // Determine the sign
		    if( outputLines.get( lineNo ).contains( "-" ) ) {
			energy = -Double.parseDouble( forEnergy ); // Negate the value for "-"
		    } else {
			energy = Double.parseDouble( forEnergy );  // Keep the value positive for "+"
		    }
		    energies.add( energy );
		    numberOfBPS.add( countLeftParens( parts[0].substring( 0, parts[0].length() - 1 ) ) );
		    states.add( outputLines.get( lineNo ) );
		    stateId++;
		    if( !this.maximizeDoubleHelix && !this.keepAll ) break; // if we're not searching the space, we have the MFE in the first state
		}

	    } catch( IOException exc ) {
		exc.printStackTrace();
	    } catch( InterruptedException exc ) {
		exc.printStackTrace();
	    }

	    if( states.size() == 0 ) theEvalue += this.e; // no states found, increase e and try again
	    else { // assign freeEnergy and state
		//Utils.debug( "fold() -> number of states: " + states.size() );
		if( !this.maximizeDoubleHelix ) { // return MFE
		    this.freeEnergy = energies.get( 0 ); // MFE is the one of the first state
		    this.state = states.get( 0 );
		    if( this.keepAll ) {
			this.states = states;
			this.energies = energies;
			this.numberOfBPS = numberOfBPS;
			createDuplexes(); // create Duplex list
		    }
		}
		else { // return the structure with longer perfect helix
		    int bestDoubleHelixIndex = 0;
		    int bestDoubleHelixSoFar = 0;
		    //System.out.println( "searching best helix..." );
		    for( int i = 0; i < numberOfBPS.size(); i++ ) {
			//System.out.println( "state: " + states.get( i ) + ", #bps: " + numberOfBPS.get( i ) );
			if( numberOfBPS.get( i ) > bestDoubleHelixSoFar ) { // update best helix
			    //System.out.println( "new state: " + states.get( i ) + ", #bps: " + numberOfBPS.get( i ) );
			    bestDoubleHelixIndex = i;
			    bestDoubleHelixSoFar = numberOfBPS.get( i );
			}
		    }
		    //if( bestDoubleHelixIndex != 0 )
			//System.out.println( "index of best helix = " + bestDoubleHelixIndex );
		    this.freeEnergy = energies.get( bestDoubleHelixIndex );
		    this.state = states.get( bestDoubleHelixIndex );
		    if( this.keepAll ) {
			this.states = states;
			this.energies = energies;
			this.numberOfBPS = numberOfBPS;
			createDuplexes(); // create Duplex list
		    }
		}
	    }
	}
	if( this.freeEnergy == 0.0 ) return -1;
	return 0;
    }
    
    // Override
    public String toString() {
	String out = "";
	// index
	// out += "Index: ";
	// // Iterating over the index keys
	// out += "{ "
	// for( Integer key : this.index.keySet() )
	//     out += key + ":" + this.index.get( key ) + " ";
	// out += "}\n";
	// // basePairs
	// out += "Bps: ";
	// // Iterating over the basePairs keys
	// out += "{ ";
	// for( Integer key : this.basePairs.keySet() )
	//     out += key + ":" + this.basePairs.get( key ) + " ";
	// out += "}\n";
	// // Iterating over the basePairTypes keys
	// out += "{ ";
	// for( Integer key : this.basePairTypes.keySet() )
	//     out += key + ":" + this.basePairTypes.get( key ) + " ";
	// out += "}\n";
	out += strand + "  " + antiStrand + "\n";
	out += this.state + ", " + this.freeEnergy;
	return out;
    }
}
