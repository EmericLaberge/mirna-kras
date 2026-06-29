/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
* SecondaryStructure is a class for generic RNA 2D structures formed by a strand
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
import java.util.LinkedHashSet;
import java.util.Stack;
import java.util.List;
import java.util.LinkedList;
import java.util.ArrayList;

import java.lang.Runtime;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;


/** -------------------------------------------
      SecondaryStructure
        is a class to handle RNA secondary structure conformational states, which
        has the following attributes:
        Double freeEnergy; // DeltaG energy
        String state = ""; // dot-bracket of the strand, initially empty
        String strand; // sequence of the strand
	double[] canonicalBasePairProbability; // base pair probabilities
	double[] reactivity; // reactivities
	Map<Integer,Integer> basePairs // base pairs
*/

public class SecondaryStructureBkp {

    // Dot bracket symbols
    public final static String parent5 = "(<{";
    public final static String parent3 = ")>}";
    public final static char csingle = '.';
    public final static String ssingle = ".";
    public final static String regxsingle = "\\.";
    public final static char cgu5 = '{';
    public final static char cgu3 = '}';
    public final static String sgu5 = "{";
    public final static String sgu3 = "}";
    public final static char ccanonical5 = '(';
    public final static char ccanonical3 = ')';
    public final static String scanonical5 = "(";
    public final static String scanonical3 = ")";
    public final static char cnonCanonical5 = '<';
    public final static char cnonCanonical3 = '>';
    public final static String snonCanonical5 = "<";
    public final static String snonCanonical3 = ">";
    public final static String sCanonicalGU5 = "({";
    public final static String sCanonicalGU3 = ")}";
    public final static String bp5 = "(<{";
    public final static String bp3 = ")>}";

    public static int countBPs( String sense, String antisense ) { // assume sense and antisense same size
	String antiantisense = Utils.reverse( antisense ); // reverse the antisense, as it is given 5'->3'
	int numBPs = 0;
	for( int i = 0; i < sense.length(); i++ )
	    if( StringSequence.canonicalGUBps.contains( Character.toString( sense.charAt( i ) ) + Character.toString( antiantisense.charAt( i ) ) ) ) numBPs++;
	return numBPs;
    }

    protected Double freeEnergy = 0.0; // of the MFE
    protected List<Double> energies; // states' energies
    protected String state = ""; // of the MFE
    protected List<String> states; // for all states
    protected String strand; // strand to fold
    protected String mask = ""; // optional mask
    protected double e; // initial -e parameter for mcff
    protected String name = ""; // name when abstract shape is used
    protected String abstractShape; // abstract shape to be considered
    protected List<String> shapes; // states' shapes
    protected boolean checkShape = false;
    protected int ultimateNumberOfStates = 0; // number in the entire conformational space
    protected int numberOfStates = 0; // number of states in the folding result
    protected int[] canonicalBasePairCount; // counts of canonical base pairs in states
    protected int[] noncanonicalBasePairCount; // counts of non-canonical base pairs in states
    protected int[] dotCount; // counts of unpaired/singlestrandedness in states
    protected double[] canonicalBasePairProbability; // probabilities of forming canonical base pairs
    protected double[] noncanonicalBasePairProbability; // probabilities of forming non-canonical base pairs
    protected double[] reactivity; // reactivities
    protected String flexibilityMask = ""; // flexibility vector in Nr format
    protected boolean folded = false; // boolean for parallell synchronization
    // protected Map<Integer,Integer> basePairs; // for future needs
    // protected Map<Integer,Character> basePairTypes; // for future needs

    // getters

    public Double  getFreeEnergy()                   { return this.freeEnergy; }
    public Double  getFreeEnergy( int i )            { return this.energies.get( i ); }
    public String  getState()                        { return this.state; }
    public List<String> getStates()                  { return this.states; }
    public String  getState( int i )                 { return this.states.get( i ); }
    public String  getStrand()                       { return this.strand; }
    public String  getMask()                         { return this.mask; }
    public double  getE()                            { return this.e; }
    public String  getAbstractShape()                { return this.abstractShape; }
    public String  getAbstractShape( int i )         { return this.shapes.get( i ); }
    public int     getCanonicalCount( int i )        { return this.canonicalBasePairCount[i]; }
    public int     getReactivityCount( int i )       { return this.dotCount[i]; }
    public double  getCanonicalProbability( int i )  { return this.canonicalBasePairProbability[i]; }
    public double  getReactivity( int i )            { return this.reactivity[i]; }
    public int     getUltimateNumberOfStates()       { return this.ultimateNumberOfStates; }
    public int     getNumberOfStates()               { return this.numberOfStates; }
    public double[] getReactivities()                { return this.reactivity; }
    public double[] getCBPProbabilities()            { return this.canonicalBasePairProbability; }
    public double[] getNCBPProbabilities()           { return this.noncanonicalBasePairProbability; }
    public double[] getBPProbabilities() {
	double[] bpP = new double[this.strand.length()];
	for( int i = 0; i < this.strand.length(); i++ )
	    bpP[i] = this.canonicalBasePairProbability[i] + this.noncanonicalBasePairProbability[i];
	return bpP;
    }
    public boolean isFolded()                        { return this.folded; }
    public String  getName()                         { return this.name; }
    //public Map<Integer,Integer> getBasePairs() { return this.basePairs; }

    // i is the nt index in the strand
    //public Integer getPartner( int i ) {
    //	return this.basePairs.get( i );
    //}

    //public Character getBasePairType( int i ) {
    //return this.basePairTypes.get( i );
    //}

    public static String computeAbstractShape( String dotBracket ) {
	String shape = dotBracket.replace( ".", "" ); // remove all dots
	int n = shape.length();
	// find base pairs in the filtered shape
	Stack<Integer> stackOfIndices = new Stack<>();
	Map<Integer,Integer> bps = new HashMap<>();
	for( int j = 0; j < shape.length(); j++ )
	    if( shape.charAt( j ) == '(' ) stackOfIndices.push( j );
	    else if( shape.charAt( j ) == ')' ) bps.put( stackOfIndices.pop(), j );
	// built the shape by removing the superfluous parentheses
	for( int i = 0; i < n-1; i++ )
	    if( bps.keySet().contains( i ) && shape.charAt( i ) != ' ' && bps.keySet().contains( i+1 ) && shape.charAt( i+1 ) != ' ' && bps.get( i ) == (bps.get( i+1 ) + 1) )
		// remove parentheses at i and bps.get( i )
		shape = shape.substring( 0, i ) + " " + shape.substring( i+1, bps.get( i ) ) + " " + shape.substring( bps.get( i )+1, shape.length() );
	return shape.replace( " ", "" );
    }

    public SecondaryStructureBkp( String strand, double e ) { // strand to fold, initial -e parameter
	//System.out.println( "SecondaryStructure( " + strand + ", " + e + " )" );
	this.strand = strand;
	this.mask = "";
	this.e = e;
	this.canonicalBasePairProbability = new double[this.strand.length()];
	this.noncanonicalBasePairProbability = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }
    
    public SecondaryStructureBkp( String strand, String mask, double e ) { // strand to fold, mask, initial -e parameter
	this.strand = strand;
	if( !this.mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.canonicalBasePairProbability = new double[this.strand.length()];
	this.noncanonicalBasePairProbability = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructureBkp( String strand, String mask, double e, String shape ) { // strand to fold, mask, initial -e parameter, and abstract shape
	this.strand = strand;
	if( !this.mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.abstractShape = shape;
	this.checkShape = true;
	this.canonicalBasePairProbability = new double[this.strand.length()];
	this.noncanonicalBasePairProbability = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructureBkp( String strand, double e, String abstractShape ) { // strand to fold, initial -e parameter and abstractShape
	this.strand = strand;
	this.mask = "";
	this.e = e;
	this.abstractShape = abstractShape;
	this.checkShape = true;
	this.canonicalBasePairProbability = new double[this.strand.length()];
	this.noncanonicalBasePairProbability = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructureBkp( String strand, double e, String mask, String abstractShape ) { // strand to fold, initial -e parameter, mask, and abstractShape
	this.strand = strand;
	if( !mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.abstractShape = abstractShape;
	this.checkShape = true;
	this.canonicalBasePairProbability = new double[this.strand.length()];
	this.noncanonicalBasePairProbability = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructureBkp( String strand, String mask, double e, String abstractShape, String highFlx, String mediumFlx, String lowFlx ) { // strand to fold, initial -e parameter, mask, abstractShape, and flexibility vectors
	this.strand = strand;
	if( !mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.abstractShape = abstractShape;
	if( !this.abstractShape.isEmpty() ) this.checkShape = true;
	this.canonicalBasePairProbability = new double[this.strand.length()];
	this.noncanonicalBasePairProbability = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.flexibilityMask = " " + Utils.convertToNrFormat( highFlx, " 8" );
	this.flexibilityMask += Utils.convertToNrFormat( mediumFlx, " 4" );
	this.flexibilityMask += Utils.convertToNrFormat( lowFlx, " 1" );
	
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    // public SecondaryStructure( String strand, double e, String abstractShape, String name ) { // strand to fold, initial -e parameter and keep abstractShape only
    // 	this.strand = strand;
    // 	this.mask = "";
    // 	this.e = e;
    // 	this.name = name;
    // 	this.abstractShape = abstractShape;
    // 	this.checkShape = true;
    // 	this.canonicalBasePairProbability = new double[this.strand.length()];
    // 	this.noncanonicalBasePairProbability = new double[this.strand.length()];
    // 	this.reactivity = new double[this.strand.length()];
    // 	this.fold(); // set this.state and this.freeEnergy
    // 	this.folded = true;
    // 	//this.buildBasePairs();
    // }

    // Forms the 2D structure of this.strand using mcff with provided -e initial parameter
    public void fold() {
	//System.out.println( "fold " + this.strand + " shape => " + this.abstractShape );
	// lists to store dotb (states) and their energies)
	Map<String,Map<String,Double>> shapedStates = new HashMap<>();
	Set<String> states = new LinkedHashSet<>(); // removes duplicates
	List<String> shapes = new LinkedList<>();
	List<Double> energies = new LinkedList<>();
	String mfeState = "";
	this.canonicalBasePairCount = new int[this.strand.length()]; // assume initialized to 0
	this.noncanonicalBasePairCount = new int[this.strand.length()]; // assume initialized to 0
	this.dotCount = new int[this.strand.length()]; // assume initialized to 0

	// ********** MC-FOLD **********
	// *****************************

	String mcff = "";
	String commandLine = "";
	double theEvalue = this.e;
	double mfe = 0.0;
	String mfeShape = "";
	int numberWithAbstractShape = 0; // number of states with the requested shape
	while( this.freeEnergy == 0.0 && theEvalue < 20 ) {
	    this.freeEnergy = 0.0;
	    numberWithAbstractShape = 0;
	    mfe = 0.0;
	    this.ultimateNumberOfStates = 0;
	    energies.clear();
	    states.clear();
	    shapes.clear();
	    // System.out.println( "fold( " + theEvalue + " )" );
	    // build the mcff command line using mask and e received as constructor's arguments
	    if( !this.mask.isEmpty() )
		mcff = "mcff -s " + this.strand + " -ns " + " -m " + this.mask + " -t " + theEvalue + this.flexibilityMask;
	    else
		mcff = "mcff -s " + this.strand + " -ns " + " -t " + theEvalue + this.flexibilityMask;

	    // save commandLine for further analyis (maybe)
	    commandLine = mcff;
	    System.out.println( mcff );
	    String[] commands = { "bash", "-c", mcff };

	    try {
		Process process = Runtime.getRuntime().exec( commands );
		BufferedReader reader =
		    new BufferedReader( new InputStreamReader( process.getInputStream() ) );
		/*
		  mcff output example: NOTE, using the -ns option does not generate the abstract shapes
		  (((((((((((((.((((..))))(((..)))))))))))))).))(((((....))))) -60.958
		  (((((((((((((.((((..))))((....))))))))))))).))(((((....))))) -60.246
		  (((((((((((((.((((..))))(((..)))))))))))))).))((((((..)))))) -60.221
		  (((((((((((((.((((..))))((....))))))))))))).))((((((..)))))) -59.510
		  ((((((((((((((.(((..))))(((..)))))))))))))).))((((((..)))))) -59.166
		  ((((((((((((.(((((..))))(((..)))))))))))))).))((((((..)))))) -59.155
		  ...
		*/
		// skip 3 lines;
		String line;
		//line = reader.readLine();
		//line = reader.readLine();
		//line = reader.readLine();
		int stateId = 0;
		// build state and energy lists and save the MFE
		while( ( line = reader.readLine() ) != null ) {
		    // mfe in element[1] (3-decimal precision)
		    String[] splitLine = line.split( " ", 3 );
		    double energy = Double.parseDouble( splitLine[1] );
		    String state = splitLine[0];
		    //String shape = splitLine[2].strip(); // when shape in the output
		    String shape = computeAbstractShape( state );
		    if( energy < mfe ) {
			mfe = energy;
			mfeState = state;
			mfeShape = shape;
		    }
		    this.ultimateNumberOfStates++;
		    // collect shapedStates
		    shapedStates.computeIfAbsent( shape, k -> new HashMap<>() ).put( state, energy ); 
		    if( this.checkShape && this.abstractShape.equals( shape ) ) {
			numberWithAbstractShape++; // will equal numberOfStatesFound when !this.checkShape
			//System.out.println( state + " " + shape + " " + energy );
			//System.out.println( shape + " equals " + this.abstractShape );
			// updtate counts
			for( int i = 0; i < state.length(); i++ ) {
			    char current = state.charAt( i );
			    if( current == ccanonical5 ||
				current == ccanonical3 ) this.canonicalBasePairCount[i]++;
			    else if( current == csingle ) this.dotCount[i]++;
			    else if( current == cnonCanonical5 ||
				     current == cnonCanonical3 ) this.noncanonicalBasePairCount[i]++;
			}
			if( !states.contains( state ) ) {
			    energies.add( energy );
			    states.add( state );
			    shapes.add( shape );
			}			
		    }
		    else if( !this.checkShape )
			if( !states.contains( state ) ) {
			    energies.add( energy );
			    states.add( state );
			    shapes.add( shape );
			}
		}
		reader.close();
	    } catch (IOException exc ) {
		exc.printStackTrace();
	    }

	    // no states found
	    if( this.checkShape && numberWithAbstractShape == 0 ) {
		theEvalue += 1; // 1.0 increment of the -e parameter
		//System.out.println( "increasing -e value to " + theEvalue + " for " + this.name );
	    }
	    else { // assign freeEnergy and state
		this.freeEnergy = mfe;
		this.state = mfeState;
		this.abstractShape = mfeShape;
		this.states = new ArrayList<String>( states );
		this.energies = new ArrayList<Double>( energies );
		this.shapes = new ArrayList<String>( shapes );
	    }
	} // end while( this.freeEnergy == 0.0 && theEvalue < 20 ) {
	
	if( this.freeEnergy == 0.0 ) Utils.stop( "Cannot fold or no conformation found: " + commandLine, 0 );
	// assign canonical bp probabilities
	this.numberOfStates = states.size();
	for( int i = 0; i < this.strand.length(); i++ ) {
	    this.canonicalBasePairProbability[i] = (double)this.canonicalBasePairCount[i] / states.size();
	    this.noncanonicalBasePairProbability[i] = (double)this.noncanonicalBasePairCount[i] / states.size();
	    this.reactivity[i] = (double)this.dotCount[i] / states.size();
	}
	this.e = theEvalue;
    }
    
    // Override
    public String toString() {
	String out = "";
	out += strand + "\n";
	out += this.state + " " + this.freeEnergy + " kcal/mol\n";
	//for( int i = 0; i < this.strand.length(); i++ )
	//    out += this.canonicalBasePairProbability[i] + " ";
	return out;
    }
}
