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
import java.util.Collections;

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
	double[] canonicalPairingProbabilities; // base pair probabilities
	double[] reactivity; // reactivities
	Map<BasePair,List<Integer>> basePairs // base pairs
*/

public class SecondaryStructure {

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

    public final static double R = 0.0019872;
    public final static double T = 310.15;

    public static int countBPs( String sense, String antisense ) { // assume sense and antisense same size
	String antiantisense = Utils.reverse( antisense ); // reverse the antisense, as it is given 5'->3'
	int numBPs = 0;
	for( int i = 0; i < sense.length(); i++ )
	    if( StringSequence.canonicalGUBps.contains( Character.toString( sense.charAt( i ) ) + Character.toString( antiantisense.charAt( i ) ) ) ) numBPs++;
	return numBPs;
    }

    protected Double freeEnergy = 0.0; // of the MFE
    protected List<Double> energies; // states' energies (mcff values)
    protected List<Double> adjustedEnergies; // adjusted states' energies
    protected String state = ""; // of the MFE
    protected String shape = ""; // of the MFE
    protected Set<String> states; // for all states
    protected String strand; // strand to fold
    protected String mask = ""; // optional mask
    protected double e; // initial -e parameter for mcff
    protected String name = ""; // name when abstract shape is used
    protected String abstractShape; // abstract shape to be considered
    protected List<String> shapes; // states' shapes
    protected boolean checkShape = false;
    protected int ultimateNumberOfStates = 0; // number in the entire conformational space
    protected int numberOfStates = 0; // number of states in the folding result
    protected int[] canonicalPairingCounts; // counts of canonical base pairs in states
    protected int[] noncanonicalPairingCounts; // counts of non-canonical base pairs in states
    protected int[] dotCounts; // counts of unpaired/singlestrandedness in states
    protected Map<BasePair,Double> basePairs = new HashMap<>(); // base pairs
    protected Map<BasePair,Integer> basePairCounts = new HashMap<>(); // counting the number of times we add energy to a base pair
    protected Map<Loop,Double> loops = new HashMap<>(); // loops
    protected Map<Loop,Integer> loopCounts = new HashMap<>(); // loop counts
    protected Map<String,Double> abstractShapes = new HashMap<>(); // abstract shapes
    protected Map<String,Map<String,Double>> shapedStates = new HashMap<>(); // map of states per shape
    protected Map<String,Set<BasePair>> statesBasePairs = new HashMap<>(); // list of base pairs per dotb
    protected Map<String,Set<Loop>> statesLoops = new HashMap<>(); // list of loops per dotb
    protected Map<Loop,List<Loop>> pseudoknots = new HashMap<>(); // map of loop interactions
    protected double Z; // partition function for the conformational space
    protected double userT = T; // temperature
    protected double RT; // Boltzmann constant * temperature
    protected double X; // to avoid overflow using very small RT
    
    protected double[] canonicalPairingProbabilities; // probabilities of forming canonical base pairs
    protected double[] noncanonicalPairingProbabilities; // probabilities of forming non-canonical base pairs
    protected double[] reactivity; // reactivities
    protected String flexibilityMask = ""; // flexibility vector in Nr format
    protected boolean folded = false; // boolean for parallell synchronization
    // protected Map<Integer,Integer> basePairs; // for future needs
    // protected Map<Integer,Character> basePairTypes; // for future needs

    // getters

    public double  getFreeEnergy()                         { return this.freeEnergy; }
    public double  getFreeEnergy( int i )                  { return this.energies.get( i ); }
    public double  getX()                                  { return this.X; }
    public String  getState()                              { return this.state; }
    public double  getStateProbability( int i )            { return this.adjustedEnergies.get( i ) / this.Z; }
    public String  getShape()                              { return this.shape; }
    public String  getShape( int i )                       { return this.shapes.get( i ); }
    public Set<String> getStates()                         { return this.states; }
    public double  getPartitionFunction()                  { return this.Z; }
    public double  getZ()                                  { return this.Z; }
    //public String  getState( int i )                 { return this.states.get( i ); }
    public String  getStrand()                             { return this.strand; }
    public String  getMask()                               { return this.mask; }
    public double  getE()                                  { return this.e; }
    public String  getAbstractShape()                      { return this.abstractShape; }
    public int     getCanonicalCount( int i )              { return this.canonicalPairingCounts[i]; }
    public int     getReactivityCount( int i )             { return this.dotCounts[i]; }
    public List<BasePair> getBasePairs()                    {
	List<BasePair> sortedBasePairs = new ArrayList<>( this.basePairs.keySet() );
	Collections.sort( sortedBasePairs );
	return sortedBasePairs;
    }
    public List<Loop> getLoops()                           {
	List<Loop> sortedLoops = new ArrayList<>( this.loops.keySet() );
	Collections.sort( sortedLoops );
	return sortedLoops;
    }
    public Map<Loop,List<Loop>> getPseudoknots()           { return this.pseudoknots; }
    public Set<String> getAbstractShapes()                 { return this.abstractShapes.keySet(); }
    public int     getNumberOfBasePairs()                  { return this.basePairs.keySet().size(); }
    public int     getNumberOfLoops()                      { return this.loops.keySet().size(); }
    public int     getNumberOfAbstractShapes()             { return this.abstractShapes.keySet().size(); }
    public double  getBasePairProbability( BasePair bp )   { return this.basePairs.get( bp ) / this.Z; }
    public int     getBasePairCount( BasePair bp )         { return this.basePairCounts.get( bp ); }
    public double  getLoopProbability( Loop l )            { return ( this.loops.get(l) == null ) ? 0.0 : this.loops.get(l) / this.Z; }
    public double  getAbstractShapeProbability( String s ) { return this.abstractShapes.get( s ) / this.Z; }
    public double  getBasePairEnergy( BasePair bp )        { return this.basePairs.get( bp ); }
    public double  getLoopEnergy( Loop l )                 { return this.loops.get( l ); }
    public int     getLoopCount( Loop l )                  { return ( this.loopCounts.get(l) == null ) ? 0 : this.loopCounts.get(l); }
    public double  getCanonicalPairingProbability( int i ) { return this.canonicalPairingProbabilities[i]; }
    public double  getReactivity( int i )                  { return this.reactivity[i]; }
    public int     getUltimateNumberOfStates()             { return this.ultimateNumberOfStates; }
    public int     getNumberOfStates()                     { return this.numberOfStates; }
    public double[] getReactivities()                      { return this.reactivity; }
    public double[] getCPProbabilities()                   { return this.canonicalPairingProbabilities; }
    public double[] getNCPProbabilities()                  { return this.noncanonicalPairingProbabilities; }
    public double[] getPProbabilities() {
	double[] bpP = new double[this.strand.length()];
	for( int i = 0; i < this.strand.length(); i++ )
	    bpP[i] = this.canonicalPairingProbabilities[i] + this.noncanonicalPairingProbabilities[i];
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


    private static final Set<String> validPairs = Set.of( "AU", "UA", "CG", "GC", "GU", "UG" );

    public static Map<Loop,List<Loop>> computePK( String sequence, List<Loop> loops ) {
        Map<Loop,List<Loop>> pkMap = new HashMap<>();

        for( int i = 0; i < loops.size(); i++ ) {
            Loop loop1 = loops.get( i );
            int start1 = loop1.position;
            int end1 = loop1.position + loop1.length;

            for( int j = i + 1; j < loops.size(); j++ ) {
                Loop loop2 = loops.get( j );
                int start2 = loop2.position;
                int end2 = loop2.position + loop2.length;

                // Get substrings
                String s1 = sequence.substring( start1, end1 );
                String s2 = sequence.substring( start2, end2 );

                // Check for possible interaction: antiparallel or parallel
                boolean interacts = canFormPK( s1, s2 ) || canFormPK( s1, new StringBuilder( s2 ).reverse().toString() );
                if( interacts ) {
                    pkMap.computeIfAbsent( loop1, k -> new ArrayList<>()).add( loop2 );
                }
            }
        }
        return pkMap;
    }

    // Checks if at least two valid base pairs exist between s1 and s2
    private static boolean canFormPK( String s1, String s2 ) {
        int count = 0;
        int len = Math.min( s1.length(), s2.length() );

        for( int i = 0; i < len; i++ ) {
            String pair = "" + s1.charAt(i) + s2.charAt(i);
            if( validPairs.contains( pair ) ) {
                count++;
                if( count >= 2 ) return true;
            }
        }
        return false;
    }

    public static String computeMotifs(
				       String dotBracket,
				       Map<String,Set<BasePair>> statesBasePairs,
				       Map<String,Set<Loop>> statesLoops ) {
	// Step 1: Remove dots and keep track of original indices
	StringBuilder filtered = new StringBuilder();
	List<Integer> originalIndices = new ArrayList<>();
	List<Integer> dotIndices = new ArrayList<>();
	for (int i = 0; i < dotBracket.length(); i++) {
	    char c = dotBracket.charAt( i );
	    if( c != '.' ) { // collect parentheses and their original indices
		filtered.append( c );
		originalIndices.add( i );
	    }
	    else // collect dot indices
		dotIndices.add( i );
	}

	// Step 2A: Find base pairs
	Stack<Integer> stack = new Stack<>();
	Map<Integer, Integer> bps = new HashMap<>();
	
	for( int i = 0; i < filtered.length(); i++ ) {
	    if( filtered.charAt(i) == '(' ) {
		stack.push( i );
	    } else if( filtered.charAt( i ) == ')' ) {
		int open = stack.pop();
		bps.put( open, i );

		// Add base pairs using original indices; assume basePairs != null; discipline the users
		int i_orig = originalIndices.get( open );
		int j_orig = originalIndices.get( i );
		statesBasePairs.computeIfAbsent( dotBracket, k -> new HashSet<>()).add( new BasePair( i_orig, j_orig ) );
	    }
	}

	// Step 2B: Find loops
        for( int i = 0; i < dotIndices.size(); i++ ) {
	    int start = dotIndices.get( i );
	    int length = 1;
            while( i + 1 < dotIndices.size() && dotIndices.get( i + 1 ) == dotIndices.get( i ) + 1 ) {
		i++;
                length++;
	    }
	    statesLoops.computeIfAbsent( dotBracket, k -> new HashSet<>()).add( new Loop( start, length ) );
	}

	// Step 3: Mark parentheses to remove
	boolean[] toRemove = new boolean[filtered.length()];
	for (int i = 0; i < filtered.length() - 1; i++) {
	    Integer j1 = bps.get(i);
	    Integer j2 = bps.get(i + 1);
	    if (j1 != null && j2 != null && j1 == j2 + 1) {
		toRemove[i] = true;
		toRemove[j1] = true;
	    }
	}

	// Step 4: Build final shape
	StringBuilder result = new StringBuilder();
	for (int i = 0; i < filtered.length(); i++) {
	    if (!toRemove[i]) {
		result.append(filtered.charAt(i));
	    }
	}

	return result.toString();
    }

    // public static String computeAbstractShape2( String dotBracket ) {
    // 	String shape = dotBracket.replace( ".", "" ); // remove all dots
    // 	int n = shape.length();
    // 	// find base pairs in the filtered shape
    // 	Stack<Integer> stackOfIndices = new Stack<>();
    // 	Map<Integer,Integer> bps = new HashMap<>();
    // 	for( int j = 0; j < shape.length(); j++ )
    // 	    if( shape.charAt( j ) == '(' ) stackOfIndices.push( j );
    // 	    else if( shape.charAt( j ) == ')' ) bps.put( stackOfIndices.pop(), j );
    // 	// built the shape by removing the superfluous parentheses
    // 	for( int i = 0; i < n-1; i++ )
    // 	    if( bps.keySet().contains( i ) && shape.charAt( i ) != ' ' && bps.keySet().contains( i+1 ) && shape.charAt( i+1 ) != ' ' && bps.get( i ) == (bps.get( i+1 ) + 1) )
    // 		// remove parentheses at i and bps.get( i )
    // 		shape = shape.substring( 0, i ) + " " + shape.substring( i+1, bps.get( i ) ) + " " + shape.substring( bps.get( i )+1, shape.length() );
    // 	return shape.replace( " ", "" );
    // }

    public SecondaryStructure( String strand, double e ) { // strand to fold, initial -e parameter
	//System.out.println( "SecondaryStructure( " + strand + ", " + e + " )" );
	this.strand = strand;
	this.mask = "";
	this.e = e;
	this.abstractShape = "";
	this.checkShape = false;
	this.canonicalPairingProbabilities = new double[this.strand.length()];
	this.noncanonicalPairingProbabilities = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.states = new LinkedHashSet<String>();
	this.energies = new ArrayList<Double>();
	this.adjustedEnergies = new ArrayList<Double>();
	this.shapes = new ArrayList<String>();
	this.Z = 0.0; // initialize partition function
	this.RT = this.userT * R;
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }
    
    public SecondaryStructure( String strand, String mask, double e ) { // strand to fold, mask, initial -e parameter
	this.strand = strand;
	this.mask = mask;
	if( !this.mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.abstractShape = "";
	this.checkShape = false;
	this.canonicalPairingProbabilities = new double[this.strand.length()];
	this.noncanonicalPairingProbabilities = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.states = new LinkedHashSet<String>();
	this.energies = new ArrayList<Double>();
	this.adjustedEnergies = new ArrayList<Double>();
	this.shapes = new ArrayList<String>();
	this.Z = 0.0; // initialize partition function
	this.RT = this.userT * R;
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructure( String strand, String mask, double e, String shape ) { // strand to fold, mask, initial -e parameter, and abstract shape
	this.strand = strand;
	if( !this.mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.abstractShape = shape;
	this.checkShape = true;
	this.canonicalPairingProbabilities = new double[this.strand.length()];
	this.noncanonicalPairingProbabilities = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructure( String strand, double e, String abstractShape ) { // strand to fold, initial -e parameter and abstractShape
	this.strand = strand;
	this.mask = "";
	this.e = e;
	this.abstractShape = abstractShape;
	this.checkShape = true;
	this.canonicalPairingProbabilities = new double[this.strand.length()];
	this.noncanonicalPairingProbabilities = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructure( String strand, double e, String mask, String abstractShape ) { // strand to fold, initial -e parameter, mask, and abstractShape
	this.strand = strand;
	if( !mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.abstractShape = abstractShape;
	this.checkShape = true;
	this.canonicalPairingProbabilities = new double[this.strand.length()];
	this.noncanonicalPairingProbabilities = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.fold(); // set this.state and this.freeEnergy
	this.folded = true;
	//this.buildBasePairs();
    }

    public SecondaryStructure( String strand, String mask, double e, String abstractShape, String highFlx, String mediumFlx, String lowFlx ) { // strand to fold, initial -e parameter, mask, abstractShape, and flexibility vectors
	this.strand = strand;
	if( !mask.isEmpty() ) this.mask = "'" + mask + "'";
	this.e = e;
	this.abstractShape = abstractShape;
	if( !this.abstractShape.isEmpty() ) this.checkShape = true;
	this.canonicalPairingProbabilities = new double[this.strand.length()];
	this.noncanonicalPairingProbabilities = new double[this.strand.length()];
	this.reactivity = new double[this.strand.length()];
	this.flexibilityMask = " " + Utils.convertToNrFormat( highFlx, " 8" );
	this.flexibilityMask += Utils.convertToNrFormat( mediumFlx, " 4" );
	this.flexibilityMask += Utils.convertToNrFormat( lowFlx, " 1" );
	this.states = new LinkedHashSet<String>();
	this.energies = new ArrayList<Double>();
	this.adjustedEnergies = new ArrayList<Double>();
	this.shapes = new ArrayList<String>();
	this.Z = 0.0; // initialize partition function
	this.RT = this.userT * R;
	this.fold();
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
    // 	this.canonicalBasePairingProbabilities = new double[this.strand.length()];
    // 	this.noncanonicalBasePairingProbabilities = new double[this.strand.length()];
    // 	this.reactivity = new double[this.strand.length()];
    // 	this.fold(); // set this.state and this.freeEnergy
    // 	this.folded = true;
    // 	//this.buildBasePairs();
    // }

    // Forms the 2D structure of this.strand using mcff with provided -e initial parameter
    //    remove duplicated states (yes, this happens with mc-flashfold)
    public void fold() {
	//System.out.println( "fold " + this.strand + " shape => " + this.abstractShape );
	// lists to store dotb (states) and their energies)
	String mfeState = "";
	this.canonicalPairingCounts = new int[this.strand.length()]; // assume initialized to 0
	this.noncanonicalPairingCounts = new int[this.strand.length()]; // assume initialized to 0
	this.dotCounts = new int[this.strand.length()]; // assume initialized to 0

	// ********** MC-FOLD **********
	// *****************************

	String mcff = "";
	String commandLine = "";
	double theEvalue = this.e;
	double mfe = 0.0;
	String mfeShape = "";
	this.freeEnergy = 0.0;
	int numberWithAbstractShape = 0; // number of states with the requested shape
	while( this.freeEnergy == 0.0 && theEvalue < 20 ) {
	    this.ultimateNumberOfStates = 0;
	    //System.out.println( "fold( " + theEvalue + " )" );
	    // build the mcff command line using mask and e received as constructor's arguments
	    if( !this.mask.isEmpty() )
		mcff = "mcff -s " + this.strand + " -ns " + " -m " + this.mask + " -t " + theEvalue + this.flexibilityMask;
	    else
		mcff = "mcff -s " + this.strand + " -ns " + " -t " + theEvalue + this.flexibilityMask;

	    // save commandLine for further analyis (maybe)
	    commandLine = mcff;
	    //System.out.println( mcff );
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

		String line;
		double energy;
		String state;
		String shape;
		// skip 3 lines in some mcff modes
		//line = reader.readLine();
		//line = reader.readLine();
		//line = reader.readLine();

		// read dotbs and energies, save in shapedStates
		while( ( line = reader.readLine() ) != null ) {
		    String[] splitLine = line.split( " ", 2 ); // dotb + energy
		    state = splitLine[0];
		    this.ultimateNumberOfStates++;
		    energy = Double.parseDouble( splitLine[1] );
		    // compute shape and add it to statesBasePairs and statesLoops
		    shape = computeMotifs( state, this.statesBasePairs, this.statesLoops );
		    if( energy < mfe ) { // adjust global mfe (we use mcff -ns)
			mfe = energy;
			mfeState = state;
			mfeShape = shape;
		    }
		    this.shapedStates.computeIfAbsent( shape, k -> new HashMap<>() ).put( state, energy );
		}
		reader.close();
	    } catch (IOException exc ) {
		exc.printStackTrace();
	    }

	    if( this.checkShape ) {
		Map<String, Double> states = this.shapedStates.get( this.abstractShape );
		numberWithAbstractShape = ( states != null ) ? states.size() : 0;
	    }
	    else numberWithAbstractShape = this.ultimateNumberOfStates; // include duplicates

	    if( numberWithAbstractShape == 0 ) {
		theEvalue += 1; // 1.0 increment of the -e parameter
		//System.out.println( "increasing -e value to " + theEvalue + " for " + this.name );
		this.shapedStates.clear(); // clear previously accumulated shapedStates
	    }
	    else { // assign freeEnergy and state
		this.freeEnergy = mfe;
		this.state = mfeState;
		this.shape = mfeShape;
		this.X = -this.freeEnergy / this.RT; // determine adjustment's factor
		this.Z = 0.0; // initial partition function

		// extract the states from shapedStates
		Set<Map.Entry<String,Double>> statesToConsider = new HashSet<>();
		if( this.checkShape ) {
		    Map<String,Double> shapeMap = this.shapedStates.get( this.abstractShape );
		    if( shapeMap != null ) {
			statesToConsider.addAll( shapeMap.entrySet() );
			for( int i = 0; i < shapeMap.size(); i++ ) this.shapes.add( this.abstractShape );
		    }
		} else {
		    // flatten all states from all shapes
		    for( Map.Entry<String,Map<String, Double>> shapeEntry : this.shapedStates.entrySet() ) {
			String shape = shapeEntry.getKey();
			Map<String,Double> subMap = shapeEntry.getValue();
			statesToConsider.addAll( subMap.entrySet() );
			for( int i = 0; i < subMap.size(); i++ ) this.shapes.add( shape ); // check this, maybe not necessary to loop the subMap?
		    }
		}

		// construct the solution set
		int stateIndex = 0;
		for( Map.Entry<String, Double> innerEntry : statesToConsider ) {
		    String dotb = innerEntry.getKey();
		    Double energy = innerEntry.getValue();
		    this.states.add( dotb ); // add state
		    this.energies.add( energy ); // add mcff energy
		    energy = Math.exp( -energy / this.RT - this.X ); // adjust for partition function
		    this.adjustedEnergies.add( energy );
		    this.Z += energy; // adjusted accumulate partition function
		    // compute the basePairs energies
		    for( BasePair bp: this.statesBasePairs.get( dotb ) ) {
			this.basePairs.merge( bp, energy, Double::sum );
			this.basePairCounts.merge( bp, 1, Integer::sum );
		    }
		    // compute the loops energies
		    for( Loop l: this.statesLoops.get( dotb ) ) {
			this.loops.merge( l, energy, Double::sum );
			this.loopCounts.merge( l, 1, Integer::sum );
		    }
		    // compute the abstract shapes energies
		    this.abstractShapes.merge( this.shapes.get( stateIndex++ ), energy, Double::sum );
		}

		// sort the base pairs and loops
		

		// stats of the pairing and total energy
		for( String state: this.states ) {
		    for( int i = 0; i < state.length(); i++ ) {
			char current = state.charAt( i );
			if( current == ccanonical5 ||
			    current == ccanonical3 ) this.canonicalPairingCounts[i]++;
			else if( current == csingle ) this.dotCounts[i]++;
			else if( current == cnonCanonical5 ||
				 current == cnonCanonical3 ) this.noncanonicalPairingCounts[i]++;
		    }
		}
	    }
	} // end while( this.freeEnergy == 0.0 && theEvalue < 20 ) {

	// compute pairing and base pair probabilities
	if( this.freeEnergy == 0.0 ) Utils.stop( "Cannot fold or no conformation found: " + commandLine, 0 );
	// assign canonical bp probabilities
	this.numberOfStates = this.states.size();
	for( int i = 0; i < this.strand.length(); i++ ) {
	    this.canonicalPairingProbabilities[i] = (double)this.canonicalPairingCounts[i] / states.size();
	    this.noncanonicalPairingProbabilities[i] = (double)this.noncanonicalPairingCounts[i] / states.size();
	    this.reactivity[i] = (double)this.dotCounts[i] / states.size();
	}

	// compute pseudoknots
	this.pseudoknots = computePK( this.strand, this.getLoops() );
	this.e = theEvalue;
    }

    // Override
    public String toString() {
	String out = "";
	out += strand + "\n";
	out += this.state + " " + this.freeEnergy + " kcal/mol\n";
	//for( int i = 0; i < this.strand.length(); i++ )
	//    out += this.canonicalPairingProbabilities[i] + " ";
	return out;
    }
}
