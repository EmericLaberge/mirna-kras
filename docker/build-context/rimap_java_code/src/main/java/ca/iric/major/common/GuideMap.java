/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * GuideMap
 *
 * Provides a guide map <name,sequence>, where name is build from transcript names and position of the kmers in map A.
 *          a guide map <name,target sites>, where target sites is a list of list of positions of kmers in map A.
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 1.0 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.LinkedHashSet;

import java.util.function.Predicate;

import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;
import java.util.Iterator;

import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;                                                                                                      
import java.util.concurrent.atomic.AtomicInteger;
import java.lang.InterruptedException;

/**
 * 
 */
public class GuideMap {

    // add the antisense generators to each grip
    //    *** size is not use as the antisense generators is only for 21-nt for now
    private static Map<String,AntisenseGenerator> designAntisensePool2( GripMap gripMap, int size, List<String> exclusions, double gcMin, double gcMax ) {
	Map<String,AntisenseGenerator> pool = new LinkedHashMap<>(); // pool container
	Predicate<String> filter = guide ->
	    guide.charAt( guide.length() - 1 ) != 'U' &&
	    !containsExclusions( guide, exclusions ) &&
	    containsRightGCPercentage( guide, gcMin, gcMax );
	for( String grip : gripMap.getRequiredGripMap().keySet() ) { // iterate through the required transcripts' grips
	    String reverseKMerA = StringSequence.reverseComplement( gripMap.get3pKey( grip ) );
	    String reverseKMerB = StringSequence.reverseComplement( gripMap.get5pKey( grip ) );
	    pool.put( grip, new AntisenseGenerator( reverseKMerA, reverseKMerB, filter ) );
	}
	return pool;
    }

    private static Map<String, AntisenseGenerator> designAntisensePool( GripMap gripMap, int size, List<String> exclusions, double gcMin, double gcMax) {

	Predicate<String> filter = guide ->
	    //guide.charAt( guide.length() - 1 ) != 'U' &&
	    !containsExclusions( guide, exclusions ) &&
	    containsRightGCPercentage( guide, gcMin, gcMax );

	return gripMap.getRequiredGripMap().keySet().parallelStream()
	    .collect( Collectors.toConcurrentMap(
						grip -> grip,
						grip -> {
						    String reverseCompKMerA = StringSequence.reverseComplement( gripMap.get3pKey( grip ) );
						    String reverseCompKMerB = StringSequence.reverseComplement( gripMap.get5pKey( grip ) );
						    return new AntisenseGenerator( reverseCompKMerA, reverseCompKMerB, filter );
						}
						));
    }

    // for each grip, filter the antisense that affect each required transcripts
    private void filterAntisensePool2() {
	Predicate<Guide> filter = guide ->
	    guide.bind() &&
	    guide.getKd() < this.maxKd;

	int gripNumber = 0;
	for( String grip : this.gripMap.getRequiredGripMap().keySet() ) { // for each grip in the required transcripts
	    //System.out.println( ++gripNumber + "/" + this.gripMap.getRequiredGripMap().keySet().size() + " grips done" );
	    //int nbAdded = 0;
	    //int nbRejected = 0;
	    System.out.println( grip + ":" );
            String kmerA = gripMap.get3pKey( grip );
            String kmerB = gripMap.get5pKey( grip );
	    Set<ProteinCodingTranscript> targetedTranscripts = new HashSet<>(); // keep track of the targeted transcripts
	    AntisenseGenerator ag = this.antisensePool.get( grip ); // get the antisense generator
	    int nbAntisense = 0;
	    while( ag.hasNext() ) { // create the set of antisense that bind and satisfy the filter for each required transcript
		String antisense = ag.next(); // test the next antisense
		nbAntisense++;
		for( Grip g : this.gripMap.getRequiredGripMap().get( grip ) ) { // for each grip location across the transcripts
		    ProteinCodingTranscript pct = g.getPCT(); // get the pct of the grip
		    if( !targetedTranscripts.contains( pct ) ) {
			int posA = g.getPosition3p();
			int posB = g.getPosition5p();
			int bridgeLen = g.getBridge();
			int tlast = posA + kmerA.length(); // target's position facing g1 (last)
			int t1 = tlast - ( this.guideSize + bridgeLen - 4 ); // target's position facing glast; bridge contains A-Box
			t1 = tlast - 30; // fix t1 to make an MRE of 31 nts                             
			if( pct.inTranscript( t1, tlast ) ) { // check the 31mer is in the transcript
			    int g2 = posA + kmerA.length() - 1;
			    int g12 = posB + kmerB.length();
			    Guide guide = new Guide( "", "design", pct, "", antisense, g2, g12, kmerA, kmerB, bridgeLen, t1, tlast, true );
			    guide.fold();
			    if( filter.test( guide ) ) {
				//nbAdded++;
				targetedTranscripts.add( pct );
				this.designs.computeIfAbsent( grip, k -> new HashSet<>()).add( antisense );
				//System.out.println( antisense );
				if( targetedTranscripts.size() == this.targets.size() ) break; // exit loop as all required transcript works for this antisense
			    }
			    else { // does not bind or KD way too high
				//nbRejected++;
				//ag.addSeedToAvoid( antisense.substring( 1, 8 ) ); // drop that seed
			    }
			}
		    }
		}
	    } // try next antisense
	    //System.out.println( nbAntisense + " tried, " + nbAdded + " added, and " + nbRejected + " were rejected " );
	} // go to next grip
	
    }

    private void filterAntisensePool3() {

	Predicate<Guide> filter = guide ->
	    guide.getKd() <= this.maxKd && // respect designer max Kd
	    guide.bind() && // make sure the guide binds ot the transcript
	    !guide.hasEarlyBulgeInTarget(); // make sure there is no bulge in target

	this.gripMap.getRequiredGripMap().keySet().parallelStream().forEach( grip -> { // for all grips common to targeted transcripts
		//System.out.println(grip + ":");
		String kmerA = gripMap.get3pKey( grip );
		String kmerB = gripMap.get5pKey( grip );
		Set<ProteinCodingTranscript> targetedTranscripts = new HashSet<>();
		AntisenseGenerator ag = this.antisensePool.get( grip );

		while( ag.hasNext() ) { // for all antisense related to this grip
		    String antisense = ag.next();
		    for( Grip g : this.gripMap.getRequiredGripMap().get( grip ) ) { // for all transcripts and position pairs (Grip instances)
			ProteinCodingTranscript pct = g.getPCT(); // get the transcript of the Grip
			if( !targetedTranscripts.contains( pct ) ) {
			    int posA = g.getPosition3p(); // get position in 3'
			    int posB = g.getPosition5p(); // get postions in 5'
			    int bridgeLen = g.getBridge();
			    int tlast = posA + kmerA.length();
			    int t1 = tlast - 30; // fixed 31mer
			    if( pct.inTranscript( t1, tlast ) ) {
				int g2 = posA + kmerA.length() - 1;
				int g12 = posB + kmerB.length();
				Guide guide = new Guide( "", "design", pct, "", antisense, g2, g12, kmerA, kmerB, bridgeLen, t1, tlast, true );
				guide.fold();
				if( filter.test( guide ) ) {
				    synchronized( targetedTranscripts ) {
					targetedTranscripts.add( pct );
				    }
				    if( targetedTranscripts.size() == this.targets.size() ) { // this antisense binds all targeted transcripts
					//System.out.println(antisense); // print the antisense
					synchronized( this.designs ) { // add this antisense to the designs
					    this.designs.computeIfAbsent( grip, k -> new HashSet<>()).add( antisense );
					}
					break; // so end the loop, this antisense is kept in the designs
				    }
					//else { // does not bind or produce low KD on one of the transcript
					//ag.addSeedToAvoid( guide.getSeed() );
				    //}
				}
			    }
			}
		    } // end for( Grip g: ...)
		    // all Grips tested for this antisense that does not bind to all trageted transcripts, so avoid this seed in the future search
		    //ag.addSeedToAvoid( antisense.substring( 1, 8 ) ); // drop that seed for this generator
		}
	    });
    }

    private void filterAntisensePool() {

	// Per-site filter on a folded Guide
	Predicate<Guide> filter = guide ->
	    guide.getKd() <= this.maxKd;      // respect designer max Kd
	    // guide.bind() &&                     // guide actually binds the transcript
	    // !guide.hasEarlyBulgeInTarget();     // no early bulge in target

	Map<String, Set<Grip>> requiredGripMap = this.gripMap.getRequiredGripMap();

	// --------- Limit global parallelism to at most 8 threads ----------
	int maxThreads = 128;
	int parallelism = Math.min(maxThreads, Runtime.getRuntime().availableProcessors());
	ExecutorService executor = Executors.newFixedThreadPool(parallelism);

	List<Callable<Void>> tasks = new ArrayList<>();

	for (String gripKey : requiredGripMap.keySet()) {
	    tasks.add(() -> {
		    processGripKey(gripKey, requiredGripMap, filter);
		    return null;
		});
	}

	try {
	    executor.invokeAll(tasks);
	} catch (InterruptedException e) {
	    Thread.currentThread().interrupt();
	    Utils.stop("filterAntisensePool interrupted", 900);
	} finally {
	    executor.shutdown();
	    try {
		executor.awaitTermination(7, TimeUnit.DAYS);
	    } catch (InterruptedException e) {
		Thread.currentThread().interrupt();
	    }
	}
    }

    /**
     * Process all antisense sequences for a single gripKey, sequentially.
     * Called in parallel for different gripKeys by filterAntisensePool().
     */
    private void processGripKey(String gripKey,
				Map<String, Set<Grip>> requiredGripMap,
				Predicate<Guide> filter) {

	// 3p and 5p kmers for this gripKey
	String kmerA = this.gripMap.get3pKey(gripKey);
	String kmerB = this.gripMap.get5pKey(gripKey);

	// Antisense generator for this gripKey (must not be shared across gripKeys)
	AntisenseGenerator ag = this.antisensePool.get(gripKey);
	if (ag == null) {
	    return; // nothing to do
	}

	// --------- Precompute per-site geometry ONCE per gripKey ----------

	// Helper class to store all site-level constants
	class Site {
	    final ProteinCodingTranscript pct;
	    final int t1;
	    final int tlast;
	    final int g2;
	    final int g12;
	    final int bridgeLen;

	    Site(ProteinCodingTranscript pct, int t1, int tlast, int g2, int g12, int bridgeLen) {
		this.pct = pct;
		this.t1 = t1;
		this.tlast = tlast;
		this.g2 = g2;
		this.g12 = g12;
		this.bridgeLen = bridgeLen;
	    }
	}
	
	// Group grips by transcript and precompute the geometry per site
	Map<ProteinCodingTranscript, List<Site>> sitesByTranscript = new LinkedHashMap<>();

	for (Grip g : requiredGripMap.get(gripKey)) {
	    ProteinCodingTranscript pct = g.getPCT();

	    int posA = g.getPosition3p();
	    int posB = g.getPosition5p();
	    int bridgeLen = g.getBridge();

	    int tlast = posA + kmerA.length();
	    int t1    = tlast - 30;         // fixed 31mer
	    int g2    = posA + kmerA.length() - 1;
	    int g12   = posB + kmerB.length();

	    sitesByTranscript
		.computeIfAbsent(pct, p -> new ArrayList<>())
		.add(new Site(pct, t1, tlast, g2, g12, bridgeLen));
	}

	// Local result set for this gripKey (deterministic order within gripKey)
	Set<String> antisensesForGrip = new LinkedHashSet<>();

	// --------- Iterate lazily over antisense sequences for this gripKey ---------
	while (ag.hasNext()) {
	    String antisense = ag.next();

	    boolean okForAllTranscripts = true;

	    // For each required transcript, we just need ONE valid site
	    for (Map.Entry<ProteinCodingTranscript, List<Site>> entry : sitesByTranscript.entrySet()) {
		ProteinCodingTranscript pct = entry.getKey();
		List<Site> sitesForPct = entry.getValue();

		boolean okForThisTranscript = false;

		for (Site s : sitesForPct) {
		    // Check transcript bounds once for precomputed geometry
		    if (!pct.inTranscript(s.t1, s.tlast)) {
			continue;
		    }

		    Guide guide = new Guide(
					    "",           // name
					    "design",     // source
					    pct,          // transcript
					    "",           // comment / label
					    antisense,    // antisense sequence
					    s.g2,
					    s.g12,
					    kmerA,
					    kmerB,
					    s.bridgeLen,
					    s.t1,
					    s.tlast,
					    true          // design mode
					    );

		    // Expensive part
		    guide.fold();

		    if (filter.test(guide)) {
			okForThisTranscript = true;
			break; // one good site for this transcript is enough
		    }
		}

		if (!okForThisTranscript) {
		    okForAllTranscripts = false;

		    // >>> This is where your "throw the grip" idea lives <<< 
		    // If you are *certain* that fail for this transcript
		    // depends only on the seed/supp geometry (same for
		    // all antisense of this gripKey), you can uncomment:
		    //
		    // Utils.debug("All sites failed for transcript " + pct.getName()
		    //             + " with grip " + gripKey + ", abandoning gripKey.");
		    // return;  // drop the whole gripKey early
		    
		    // Otherwise, we just discard this *antisense* and try the next one.
		    break;
		}
	    }

	    // This antisense has at least one valid site in *every* required transcript
	    if (okForAllTranscripts) {
		//Utils.debug( antisense + " works on all transcripts using " + gripKey );
		antisensesForGrip.add(antisense);
	    }
	} // end while antisense sequence

	// Merge local results into global designs map once per gripKey
	if (!antisensesForGrip.isEmpty()) {
	    synchronized (this.designs) {
		this.designs
		    .computeIfAbsent(gripKey, k -> new LinkedHashSet<>())
		    .addAll(antisensesForGrip);
	    }
	}
    }

    // Attributes
    KMerMap kmerMap3p;
    KMerMap kmerMap5p;
    GripMap gripMap;
    Map<String,AntisenseGenerator> antisensePool;
    Map<String,Set<String>> designs;
    Set<ProteinCodingTranscript> targets;
    Set<ProteinCodingTranscript> excludedSubset;
    Map<String,Guide> guides = new LinkedHashMap<>(); // a map for <guide name,guide info> of the guides
    Map<String,Set<Guide>> grips = new LinkedHashMap<>(); // a map to keep the guides under their kmer prefix grips: <"seedPrefix/suppPrefix",Set<guide>>
    Map<String,Set<Guide>> sequences = new LinkedHashMap<>(); // a map to keep the guide sequences linked to their set of guides: <guide sequence,set of guides>
    int distance = 0; // distance constraint between sites in maps A and B
    int guideSize;
    double maxKd = 1000000;
    List<String> exclusions;
    boolean siRNAInAllTargets = false;
    static int guideReferenceNumber = 0; // global guide reference number to avoid guide key duplicates
    private double gcPercentMin = 0.3;
    private double gcPercentMax = 0.64;
    private int guideAdded = 0;
    private int guideRemoved = 0;

    // Constructor
    // Build guides from grips of the required and optional transcripts (gripMap is not empty)
    public GuideMap( GripMap gripMap, int guideSize, double gcPercentMin, double gcPercentMax, double maxKd ) {
	//System.out.println( "gripMap:\n" + gripMap );
	this.gripMap = gripMap;
	this.guideSize = guideSize;
	this.maxKd = maxKd;
	this.gcPercentMin = gcPercentMin;
	this.gcPercentMax = gcPercentMax;
	this.distance = gripMap.getMaxDistance();
	this.exclusions = this.gripMap.getExclusions();
	this.targets = new LinkedHashSet<>( this.gripMap.getRequired() );
	this.designs = new LinkedHashMap<>();
	System.out.println( "building antisensePool..." );
	this.antisensePool = designAntisensePool( this.gripMap, this.guideSize, this.exclusions, this.gcPercentMin, this.gcPercentMax );
	// check the number of sequences for each grip
	// for( String grip: this.antisensePool.keySet() ) {
	//     int numberOfAntisense = 0;
	//     while( this.antisensePool.get( grip ).hasNext() ) numberOfAntisense++;
	//     System.out.println( grip + " has " + numberOfAntisense + " antisenses" );
	// }
	System.out.println( "filtering antisensePool... " );
	this.filterAntisensePool(); // we have the grips for which we have at least one antisense that binds all required transcripts with Kd <= maxKd
	System.out.println( "done." );
	
	// here, we have generators for all grips, to be tried in each transcript where the grip is present
	
	// this.extendSeeds( gripMap.getK3() ); // extend the seed region
	// this.extendG12(); // extend G12
	// this.unfoldedGuides();
	// this.extendSupp( gripMap.getK5() ); // extend the supplementary region
	// this.unfoldedGuides();
	// this.crossHybridize(); // cross hybridize on all transcripts: required, optional, and excluded
	// this.unfoldedGuides(); 
	// this.foldAndFilter();  // fold the remaining duplexes
	// this.unfoldedGuides(); 
    }

    // here one for grip analysis only
    public GuideMap( GripMap gripMap ) {
	//System.out.println( "gripMap:\n" + gripMap );
	this.guideSize = 21; // default for coupling assessment
	this.distance = gripMap.getMaxDistance();
	this.exclusions = gripMap.getExclusions();
	this.targets = new LinkedHashSet<>( gripMap.getRequired() );
	// design siRNA for all grips
	// 1. merge the grips: required+optional+excluded
	Map<String,Set<Grip>> allGrips = new LinkedHashMap<>( gripMap.getRequiredGripMap() );
	if( allGrips.size() > 0 ) { // common grips for the required
	    // adding the guides from optional
	    for( Map.Entry<String,Set<Grip>> entry : gripMap.getOptionalGripMap().entrySet() ) {
		Set<Grip> fromAllGrips = allGrips.get( entry.getKey() ); // grips are garanteed to be in the required grip set, but they can be empty
		fromAllGrips.addAll( entry.getValue() );
		allGrips.put( entry.getKey(), fromAllGrips );
	    }
	    for( Map.Entry<String,Set<Grip>> entry : gripMap.getExcludedGripMap().entrySet() ) {
		Set<Grip> fromAllGrips = allGrips.get( entry.getKey() ); // grips are garanteed to be in the required grip set, but they can be empty
		fromAllGrips.addAll( entry.getValue() );
		allGrips.put( entry.getKey(), fromAllGrips );
	    }
	    // 2. generate the sequences for all grips and MREs
	    //    this will generate sequences (guides/duplexes) in all transcripts: required+optional+excluded
	    for( String grip : allGrips.keySet() ) {
		String kmerA = gripMap.get3pKey( grip );
		String kmerB = gripMap.get5pKey( grip );
		for( Grip g : allGrips.get( grip ) ) {
		    int posA = g.getPosition3p();
		    int posB = g.getPosition5p();
		    int bridgeLen = g.getBridge();
		    int tlast = posA + kmerA.length(); // target's position facing g1 (last)
		    int t1 = tlast - ( this.guideSize + bridgeLen - 4 ); // target's position facing glast; bridge contains A-Box
		    String sirna = g.getPCT().getAntisense( tlast, t1, bridgeLen, gripMap.getRegion(), this.guideSize );
		    //System.out.println( "MRE length: " + (tlast - t1 + 1) );
		    if( sirna != null && checkSequence( sirna ) && ( tlast - t1 ) >= 30 ) {
			String newGrip = kmerA + "/" + kmerB;
			String guideId = newGrip + "[" + g.getPCT().getName() + "." + posA + "." + posB + "]." + guideReferenceNumber++;
			int g2 = posA + kmerA.length() - 1;
			int g12 = posB + kmerB.length();
			Guide newGuide = new Guide( guideId, "design", g.getPCT(), "", sirna, g2, g12, kmerA, kmerB, bridgeLen, t1, tlast, true );
			this.addGuide( newGuide );
		    }
		}
	    }
	}
    }

    
    // ceate a guide map for a set of targets (> 1)
    public GuideMap( KMerMap kmerMap3p, KMerMap kmerMap5p, List<ProteinCodingTranscript> required, List<ProteinCodingTranscript> optional, List<ProteinCodingTranscript> excluded, int distance, int guideSize, List<String> exclusions, double gcPercentMin, double gcPercentMax ) {
	guideReferenceNumber = 0;
	this.kmerMap3p = kmerMap3p;
	this.kmerMap5p = kmerMap5p;
	this.targets = new LinkedHashSet<>( required );
	if( this.targets.size() > 1 ) {
	    this.excludedSubset = new LinkedHashSet<>( excluded );
	    this.distance = distance;
	    this.guideSize = guideSize;
	    this.exclusions = exclusions;
	    this.gcPercentMin = gcPercentMin;
	    this.gcPercentMax = gcPercentMax;
	    this.buildMapAB(); // A -> seed; B -> supp
	    this.extendSeeds( this.kmerMap3p.getK() ); // extend the seed region
	    this.extendG12(); // extend G12
	    this.unfoldedGuides();
	    this.extendSupp( kmerMap5p.getK() ); // extend the supplementary region
	    this.unfoldedGuides();

	    this.crossHybridize(); // cross hybridize all design sequences
	    this.unfoldedGuides();

	    this.addOptionalTranscripts( optional );
	    this.unfoldedGuides();
	    this.foldAndFilter();
	    this.unfoldedGuides();
	}
	else
	    Utils.stop( "GuideMap for 1 target or less is not possible", 0 );	    
    }

    // utilities

    private void unfoldedGuides() {
	int count = 0;
	for( Guide g : this.guideSet() )
	    if( !g.isFolded() ) count++;
	System.out.println( "=> " + this.size() + " duplexes, " + count + " are unfolded" );
	count = 0;
	int countDuplexes = 0;
	for( String s : this.getSequences() )
	    for( Guide g : this.sequences.get( s ) ) {
		countDuplexes++;
		if( !g.isFolded() ) 
		    count++;
	    }
	System.out.println( "=> " + this.numberOfSequences() + " sequences (" + countDuplexes + " duplexes), " + count + " are unfolded" );
	count = 0;
	countDuplexes = 0;
	for( String grip : this.getGrips() )
	    for( Guide g : this.guideSetFromGrip( grip ) ) {
		countDuplexes++;
		if( !g.isFolded() ) 
		    count++;
	    }
	System.out.println( "=> " + this.grips.size() + " grips (" + countDuplexes + " duplexes), " + count + " are unfolded" );
    }
    
    private void addGuide( Guide g ) {
	this.guides.put( g.getId(), g ); // add to the guide map
	this.addToGrips( g ); // add to the grip map
	this.addToSequences( g ); // add to the sequence map
	this.guideAdded++;
	//this.checkNumbersOfGuides();
	//System.out.println( "added guide: " + g );
    }

    private void removeGuide( Guide g ) {
	this.guideRemoved++;
	this.removeFromGrips( g );
	this.removeFromSequences( g );
	this.guides.remove( g.getId() );
    }

    private void clear() {
        this.guides = new LinkedHashMap<>();
	this.grips = new LinkedHashMap<>();
	this.sequences = new LinkedHashMap<>();
    }

    // public methods
    
    public boolean hasSiRNAInAllTargets() { return this.siRNAInAllTargets; }
    public int numberOfSequences() { return this.sequences.size(); } // number of guide sequences
    public int size() { return this.guides.size(); } // number of duplexes
    
    public Set<Guide> guideSetFromSequence( String sequence ) { return this.sequences.get( sequence ); }
    
    public Set<Map.Entry<String,Guide>> entrySet() { return this.guides.entrySet(); }
    public Set<String> keySet() { return this.guides.keySet(); }
    public Guide get( String key ) { return this.guides.get( key ); }
    public Set<Guide> guideSet() { return new LinkedHashSet<Guide>( this.guides.values() ); }
    public Set<String> getSequences() { return this.sequences.keySet(); } // number of guide sequences
    public Set<String> getGrips() { return this.grips.keySet(); }
    public Set<Guide> guideSetFromGrip( String grip ) { return this.grips.get( grip ); }
    public Set<ProteinCodingTranscript> getTargetSubset() { return this.targets; }
    public int getNumberOfTargets() { return this.targets.size(); }
    public int getGuideAdded() { return this.guideAdded; }
    public int getGuideRemoved() { return this.guideRemoved; }
    public GripMap getGripMap() { return this.gripMap; }
    public Map<String,Set<String>> getDesigns() { return this.designs; }
    public double averageGuidesPerGrip() {
	double avg = 0;
	for( String grip : this.getGrips() ) // for each grip
	    avg += this.grips.get( grip ).size();
	return avg / this.getGrips().size();
    }

    // utilities

    // get seed from grip
    private String getSeedFromGrip( String grip ) {
	return grip.substring( 0, this.kmerMap3p.getK() );
    }

    private String getSuppFromGrip( String grip ) {
	return grip.substring( this.kmerMap5p.getK() + 1, grip.length() );
    }

    // remove the grips that are not in all transcripts
    //   to be called after initial siRNA map (buildMapAB)
    private void removeIncompleteGrips() {
	System.out.print( " removing incomplete grips... " );
	for( String grip : this.getGrips() ) {
	    // build the sets of transcripts and guides for this grip
	    Set<Guide> thisGripGuides = new HashSet<>();
	    Set<CodingTranscript> thisGripTranscripts = new LinkedHashSet<>();
	    for( Guide g : this.guideSetFromGrip( grip ) ) {
		thisGripGuides.add( g );
		thisGripTranscripts.add( g.getCT() );
	    }
	    if( !thisGripTranscripts.containsAll( this.targets ) ) { // if the grip does not exist in all targets
		for( Guide g : thisGripGuides ) this.removeGuide( g ); // remove the guides with this grip
	    }
	}
	// remove the grips that have no guides
	Set<String> gripsToBeRemoved = new HashSet<>();
	for( String grip : this.getGrips() )
	    if( this.guideSetFromGrip( grip ).size() == 0 )
		gripsToBeRemoved.add( grip );
	// remove the grips to be removed; this is needed to avoid concurrent exception
	this.grips.keySet().removeAll( gripsToBeRemoved );
    }

    // add a Guide element to the grips set
    public void addToGrips( Guide g ) {
	//System.out.println( "addToGrips( " + g.getGrip() + " )" );
	// Check if the grip already exists in the map
	String grip = g.getGrip();
	if( !this.grips.containsKey( grip ) ) {
	    // create a new Set
	    this.grips.put( grip, new HashSet<>() );
	}
	// add guide
	this.grips.get( grip ).add( g );
    }

    public void removeFromGrips( Guide g ) {
	String grip = g.getGrip();
	this.grips.get( grip ).remove( g );
    }

    public void removeFromSequences( Guide g ) {
	String sequence = g.getSequence();
	if( this.sequences.containsKey( sequence ) ) { // perhaps the sequence was remove previously
	    this.sequences.get( sequence ).remove( g ); // Guide g must be in matches (as it was added previously)
	    if( this.sequences.get( sequence ).size() == 0 )
		this.sequences.remove( sequence ); // remove associated sequence if it binds nowhere
	}
    }

    // add a Guide element to the grips set
    public void addToSequences( Guide g ) {
	//System.out.println( "addToMatches( " + g.getSequence() + " )" );
	// Check if the grip already exists in the map
	String sequence = g.getSequence();
	if( !this.sequences.containsKey( sequence ) ) {
	    // create a new Set
	    this.sequences.put( sequence, new LinkedHashSet<>() );
	}
	// add guide
	this.sequences.get( sequence ).add( g );
    }

    private void removeIncompleteSequences4() {
	System.out.print(" removing incomplete sequences... ");

	// Thread-safe set to collect valid sequences
	Set<String> validSequences = ConcurrentHashMap.newKeySet();

	// Process guides in parallel
	this.guideSet().parallelStream()
	    .filter(g -> this.targets.contains(g.getCT())) // Consider only guides in targets
	    .peek(Guide::fold) // Fold the duplex
	    .filter(this::checkCondition) // Apply the checkCondition
	    .forEach(g -> validSequences.add(g.getSequence()));

	System.out.println("Valid sequences found: " + validSequences.size());

	// Identify sequences to be removed
	Set<String> allSequences = new HashSet<>(this.getSequences());
	allSequences.removeAll(validSequences); // Sequences without valid guides

	// Collect guides associated with invalid sequences
	List<Guide> guidesToBeRemoved = this.guideSet().stream()
	    .filter(g -> allSequences.contains(g.getSequence()))
	    .collect(Collectors.toList());

	System.out.println("Guides to be removed: " + guidesToBeRemoved.size());

	// Remove invalid guides
	synchronized (this) {
	    guidesToBeRemoved.forEach(this::removeGuide);
	}
    }

    private void removeIncompleteSequences() {
	System.out.print(" removing incomplete sequences... ");
    
	// Thread-safe map to organize guides by ProteinCodingTranscript
	Map<ProteinCodingTranscript, Set<Guide>> guidesByPCT = new ConcurrentHashMap<>();
	this.guideSet().parallelStream()
	    .filter(g -> this.targets.contains(g.getCT())) // add only guides in targetSet
	    .forEach(g -> guidesByPCT
		     .computeIfAbsent((ProteinCodingTranscript) g.getCT(), k -> ConcurrentHashMap.newKeySet())
		     .add(g)
		     );

	System.out.println( "PCTs considered: " + guidesByPCT.keySet() );

	// Thread-safe queue to hold guides marked for removal
	ConcurrentLinkedQueue<Guide> guidesToBeRemoved = new ConcurrentLinkedQueue<>();

	// Process each sequence in parallel
	this.getSequences().parallelStream().forEach(sequence -> {
		//System.out.println( sequence );
		Set<CodingTranscript> thisSequenceTranscripts = ConcurrentHashMap.newKeySet();
		AtomicBoolean sequenceInvalid = new AtomicBoolean(false);

		// Check each PCT for the current sequence
		guidesByPCT.keySet().forEach(pct -> {
			if (sequenceInvalid.get()) return;  // Skip if already marked for removal

			boolean foundAtLeastOneMRE = guidesByPCT.get(pct).stream()
			    .filter(g -> g.getSequence().equals(sequence))
			    .peek(Guide::fold)  // Fold the duplex
			    .anyMatch(g -> checkCondition(g) && thisSequenceTranscripts.add(pct));

			if (!foundAtLeastOneMRE) {
			    // Mark all guides with this sequence for removal and set sequenceInvalid
			    guidesToBeRemoved.addAll(this.guideSetFromSequence(sequence));
			    sequenceInvalid.set(true);  // Mark sequence as invalid
			}
		    });
	    });

	// Single synchronized removal step
	synchronized (this) {
	    guidesToBeRemoved.forEach(this::removeGuide);
	}
    }

    // remove the sequences that do not bind all targets
    //   to be called after cross hybridization (assume duplexes might not be folded
    private void removeIncompleteSequences2() {
	System.out.print( " removing incomplete sequences... " );
	Set<Guide> guidesToBeRemoved = new HashSet<>();
	// organize the guides per transcripts to optimize the procedure
	Map<ProteinCodingTranscript,Set<Guide>> guidesByPCT= new HashMap<>();
	for( Guide g : this.guideSet() ) guidesByPCT.computeIfAbsent( (ProteinCodingTranscript)(g.getCT()), k -> new HashSet<>() ).add( g );
	
	// for all sequences, check if it can bind to at least one MRE on each this.targets
	for( String sequence : this.getSequences() ) {
	    Set<CodingTranscript> thisSequenceTranscripts = new HashSet<>();
	    for( ProteinCodingTranscript pct : guidesByPCT.keySet() ) {
		Iterator<Guide> gIterator = guidesByPCT.get( pct ).iterator();
		boolean foundAtLeastOneMRE = false;
		while( gIterator.hasNext() && !foundAtLeastOneMRE && !thisSequenceTranscripts.equals( this.targets ) ) {
		    Guide g = gIterator.next();
		    if( g.getSequence().equals( sequence ) ) { // only for the guides of the running sequence
			// fold the duplex and check binding
			g.fold();
			if( checkCondition( g ) ) {
			    thisSequenceTranscripts.add( g.getCT() );
			    foundAtLeastOneMRE = true;
			}
		    }
		}
		if( !foundAtLeastOneMRE ) {
		    guidesToBeRemoved.addAll( this.guideSetFromSequence( sequence  ) ); // mark the sequence to be removed
		    break; // one pct has no MRE, this sequence will be removed, go to the next sequence
		}
	    }
	}
	// remove the guides under the sequences to be removed
	for( Guide g : guidesToBeRemoved ) this.removeGuide( g );
    }
    
    public Set<Guide> bindsExcluded( Set<ProteinCodingTranscript> excludedSubset ) {
	Set<Guide> binders = new HashSet<>();
	for( Guide g : this.guideSet() )
	    if( bindsTo( g, excludedSubset ) ) binders.add( g );
	return binders;
    }

    private boolean bindsTo( Guide g, Set<ProteinCodingTranscript> excluded ) {
	String seedPrefix = g.getSeedPrefix(); // get seed prefix of the guide
	String suppPrefix = g.getSuppPrefix(); // get supp prefix of the guide
	// create guides where the seed and supp appear
	for( ProteinCodingTranscript excludedTranscript : excluded ) {
	    String targetRegion = excludedTranscript.getTargetableSequence(); // get the sequence to possibly bind to
	    // get positions of seed in exluded transcript
	    List<Integer> positionsSeed = new ArrayList<>();
	    List<Integer> positionsSupp = new ArrayList<>();
	    int i = 0;
	    while( i < targetRegion.length() ) {
		int index = targetRegion.indexOf( seedPrefix );
		if( index > -1 ) {
		    positionsSeed.add( index );
		    index++;
		}
		else break;
	    }
	    i = 0;
	    while( i < targetRegion.length() ) {
		int index = targetRegion.indexOf( suppPrefix );
		if( index > -1 ) {
		    positionsSupp.add( index );
		    index++;
		}
		else break;
	    }
	    for( Integer posSeed : positionsSeed ) { // for each position of the seed in excluded transcript
		for( Integer posSupp : positionsSupp ) { // for each position of the supp in excluded transcript
		    int bridgeLen = ( posSeed - ( 7 - seedPrefix.length() ) ) - ( posSupp + suppPrefix.length() ) - 1; // pos(facing g8) - pos(facing g12) - 1
		    if( bridgeLen >= 4 && bridgeLen <= this.distance )  { // B site within distance; bridgeLength: A-Box size of 3 + bulge
			// create new guide 
			int tlast = posSeed + 4; // target's position facing g1 (last)
			int t1 = tlast - ( this.guideSize + bridgeLen - 4 ); // target's position facing glast; bridge contains A-Box
			String guideRNA = g.getSequence();
			int g2 = posSeed + seedPrefix.length() - 1;
			int g12 = posSupp + suppPrefix.length();
			Guide tmpGuide = new Guide( "tmp", "design", excludedTranscript, "", guideRNA, g2, g12, seedPrefix, suppPrefix, bridgeLen, t1, tlast, true );
			tmpGuide.fold();
			if( checkCondition( tmpGuide ) ) return true;
		    }
		}
	    }
	}
	return false; // if the guide was not found to fold any of the excluded transcript
    }

    // check that the number of guides in guides equals the number of guides in matches
    public boolean checkNumbersOfGuides() {
	int numberOfGuidesUnderSequences = 0;
	for( String seq : this.getSequences() )
	    numberOfGuidesUnderSequences += this.guideSetFromSequence( seq ).size();
	if( numberOfGuidesUnderSequences != this.guides.size() )
	    System.out.println( "imbalance " + numberOfGuidesUnderSequences + " duplexes under sequences, and " + this.guides.size() + " duplexes in guides" );
	return numberOfGuidesUnderSequences == this.guides.size();
    }

    
    private static boolean containsExclusions( String guide, List<String> exclusions ) {
	for( String exclusion : exclusions )
	    if( guide.contains( exclusion ) ) return true;
	return false;
    }

    // accept between 30 and 64%
    private static boolean containsRightGCPercentage( String sirna, double gcMin, double gcMax ) {
        // Count G and C nucleotides
        int gcCount = 0;
        for( char nucleotide : sirna.toCharArray() ) {
            if( nucleotide == 'G' || nucleotide == 'C' ) {
                gcCount++;
            }
        }
        // Calculate the percentage of GC content
	double gcPercentage = gcCount / (double) sirna.length();
	// Return true if GC content is between gcPercentMin and gcPercentMax, inclusive
        return  gcPercentage >= gcMin && gcPercentage <= gcMax;
    }

    // add guides for an optional transcript
    public void addOptionalTranscripts( List<ProteinCodingTranscript> pcts ) {
	System.out.print( "adding optional transcripts (step 1) ... " );
	long startTime = System.currentTimeMillis(); // to measure time
	// requires 2 steps: 1) create sequences from optional targets in the context of the current grips
	Set<String> newSiRNAs = new HashSet<>(); // keep the new siRNA ids
	for( ProteinCodingTranscript pct : pcts ) { // for all optional targets
	    System.out.println( pct.getName() + "..." );
	    for( String grip : this.getGrips() ) { // for current grips
		System.out.println( "grip: " + grip );
		String kmerA = this.getSeedFromGrip( grip ); // get kmerA
		String kmerB = this.getSuppFromGrip( grip ); // get kmerB
		int lenKmerA = kmerA.length();
		int lenKmerB = kmerB.length();
		List<Integer> positionsA = this.kmerMap3p.getPositions( kmerA, pct ); // get positions of kmerA in pct
		List<Integer> positionsB = this.kmerMap5p.getPositions( kmerB, pct ); // get positions of kmerB in pct
		System.out.println( "kmerA: " + kmerA + ", " + positionsA );
		System.out.println( "kmerB: " + kmerB + ", " + positionsB );
		if( !positionsA.isEmpty() && !positionsB.isEmpty() ) { // if kmerA and kmerB were found in this target, try every kmerA/kmerB combinations
		    for( Integer posA : positionsA ) { // for each position in kmer map A
			for( Integer posB : positionsB ) { // for each position in kmer map B
			    int bridgeLen = ( posA - ( 7 - lenKmerA ) ) - ( posB + lenKmerB ) - 1; // pos(facing g8) - pos(facing g12)
			    if( bridgeLen >= 4 && bridgeLen <= this.distance )  { // B site within distance; bridgeLength: A-Box size of 3 + bulge
				// design siRNA in pct at posA and posB
				int tlast = posA + lenKmerA; // target's position facing g1 (last)
				int t1 = tlast - ( this.guideSize + bridgeLen - 4 ); // target's position facing glast; bridge contains A-Box
				String sirna = pct.getAntisense( tlast, t1, bridgeLen, this.kmerMap3p.getRegion(), this.guideSize ); // region of both kmerMaps are the same
				Set<String> extendedSiRNAs = this.generateAllExtensions( sirna );
				for( String si : extendedSiRNAs )
				    if( checkSequence( si ) ) {
					if( !this.sequences.containsKey( si ) ) { // if not already existing
					    int g2 = posA + lenKmerA - 1;
					    int g12 = posB + lenKmerB;
					    String newGuideId = grip + "[" + pct.getName() + "." + posA + "." + posB + "]." + guideReferenceNumber++;
					    Guide newGuide = new Guide( newGuideId, "design", pct, "", si, g2, g12, kmerA, kmerB, bridgeLen, t1, tlast, true );
					    this.addGuide( newGuide );
					    newSiRNAs.add( newGuideId );
					    System.out.println( "added: " + newGuideId + " seq: " + si);
					}
				    }
			    }
			}
		    } // all posB
		} // all posA
	    } // end if kmerA and kmerB found in pct
	// here, we have new siRNAs from the optional targets
	// try to cross hybridize them with the others
	}

	System.out.println( "added " + newSiRNAs.size() + " siRNA from optional targets" );
	int count = 0;
    	Set<Guide> newGuides = new LinkedHashSet<>(); // to collect the new duplexes
	for( String newSiId : newSiRNAs ) { // for each added siRNA
	    Guide g = this.get( newSiId ); // guide to hybridize
	    String grip = g.getGrip(); // grip for this new siRNA
	    Set<String> locations = new HashSet<>(); // keep track of locations
	    for( Guide other : this.guideSetFromGrip( grip ) ) { // try to bind the newSi to all other duplexes (with same grip)
		if( g != other ) { // don't duplicate same guide
		    if( !locations.contains( other.getLocation() ) ) { // don't bind to the same location twice
			locations.add( other.getLocation() );
			String id = g.getId() + "[" + other.getLocation() + "]" + guideReferenceNumber++;
			Guide newGuide = new Guide( other, id, "miRDesign", g.getSequence() );
			newGuides.add( newGuide ); // add new guide to the set of new guides
			count++;
		    }
		}
	    }
	}
	System.out.println( "after cross-hybridization, adding " + newGuides.size() + " guides (count=" + count + ")" );
	for( Guide g : newGuides )
	    this.addGuide( g ); // add new guides to the data structure (maps)
	this.foldAndFilter(); // filter unbinders
	this.removeIncompleteSequences(); // remove sequences that don't bind all required targets
	this.unfoldedGuides();
	System.out.println( "done in " + ( ( System.currentTimeMillis() - startTime ) / 1000 ) + " seconds" );

	// // step 2
	// long startTime = System.currentTimeMillis(); // to measure time spent in buildMapAB
	// Set<String> sequences = new HashSet<>(); // treat each sequence one time
	// int nbAdded = 0;
	// for( Guide guide : this.guideSet() ) // for all duplexes
	//     for( ProteinCodingTranscript pct : pcts ) { // see if we can bind it to any optional transcript
	// 	String kmerA = guide.getSeedPrefix(); // get kmerA (seed starts at g2, we want g2, index starts at 0 so 1)
	// 	String kmerB = guide.getSuppRegion(); // get kmerB (supp starts at g12, we want g13, index starts at 0 so 12)
	// 	int lenKmerA = kmerA.length();
	// 	int lenKmerB = kmerB.length();
	// 	List<Integer> positionsA = this.kmerMap.getSeedPositions( kmerA, pct ); // get positions of kmerA in pct
	// 	List<Integer> positionsB = this.kmerMap.getSuppPositions( kmerB, pct ); // get positions of kmerB in pct
	// 	// see if this guide can be gripped to pct
	// 	if( !positionsA.isEmpty() && !positionsB.isEmpty() && !sequences.contains( guide.getSequence() ) ) { // if kmerA and kmerB were found in this target, try every kmerA/kmerB combinations
	// 	    for( Integer posA : positionsA ) { // for each position in kmer map A
	// 		for( Integer posB : positionsB ) { // for each position in kmer map B
	// 		    int bridgeLen = ( posA - ( 7 - lenKmerA ) ) - ( posB + lenKmerB ) - 1; // pos(facing g8) - pos(facing g12)
	// 		    if( bridgeLen >= 4 && bridgeLen <= this.distance )  { // B site within distance; bridgeLength: A-Box size of 3 + bulge
	// 			// apply guide/siRNA to kmerA and kmerB in pct
	// 			int tlast = posA + lenKmerA; // target's position facing g1 (last)
	// 			int t1 = tlast - ( this.guideSize + bridgeLen - 4 ); // target's position facing glast; bridge contains A-Box
	// 			int g2 = posA + lenKmerA - 1;
	// 			int g12 = posB + lenKmerB;
	// 			String newGuideId = guide.getId() + "[" + pct.getName() + "." + posA + "." + posB + "]." + guideReferenceNumber++;
	// 			Guide newGuide = new Guide( newGuideId, "design", pct, "", guide.getSequence(), g2, g12, kmerA, kmerB, bridgeLen, t1, tlast, true );
	// 			nbAdded++;
	// 			this.addGuide( newGuide );
	// 			sequences.add( guide.getSequence() );
	// 		    }
	// 		} // all posB
	// 	    } // all posA
	// 	} // end if kmerA and kmerB found in pct
	//     } // all pcts
	// System.out.println( " done in " + ( ( System.currentTimeMillis() - startTime ) / 1000 ) + " seconds" );
    }

    // the kmers are common to all input transcripts, for each transcript, design siRNAs for each pair A, B that satisfy the bridge constraint
    private void buildMapAB() {
	System.out.print( "building the initial set of siRNAs... " );
	long startTime = System.currentTimeMillis(); // to measure time spent in buildMapAB
	Set<ProteinCodingTranscript> excluded = new HashSet<ProteinCodingTranscript>(); // empty set
	for( String kmerA : this.kmerMap3p.getExclusive( this.targets, excluded ) ) { // for each kmer in the kmer map A
	    int lenKmerA = kmerA.length();
	    for( String kmerB : this.kmerMap5p.getExclusive( this.targets, excluded ) ) { // for each kmer in the kmer map B
		int lenKmerB = kmerB.length();
		for( ProteinCodingTranscript targeti : this.targets ) { // for each target
		    List<Integer> positionsA = this.kmerMap3p.getPositions( kmerA, targeti ); // get positions of kmerA in transcript
		    List<Integer> positionsB = this.kmerMap5p.getPositions( kmerB, targeti ); // get positions of kmerB in transcript
		    if( !positionsA.isEmpty() && !positionsB.isEmpty() ) { // if kmerA and kmerB were found in this target, try every combintions
			for( Integer posA : positionsA ) { // for each position in kmer map A
			    for( Integer posB : positionsB ) { // for each position in kmer map B
				int bridgeLen = ( posA - ( 7 - lenKmerA ) ) - ( posB + lenKmerB ); // pos(facing g8) - pos(facing g12)
				if( bridgeLen >= 4 && bridgeLen <= this.distance )  { // B site within distance; bridgeLength: A-Box size of 3 + bulge
				    // design an siRNA for kmerA and kmerB
				    int tlast = posA + lenKmerA; // target's position facing g1 (last)
				    int t1 = tlast - ( this.guideSize + bridgeLen - 4 ); // target's position facing glast; bridge contains A-Box
				    String sirna = targeti.getAntisense( tlast, t1, bridgeLen, this.kmerMap5p.getRegion(), this.guideSize ); // both kmerMap has the same region
				    // check if sirna satisfy sequence constraints
				    if( sirna != null && checkSequence( sirna ) ) {
					String grip = kmerA + "/" + kmerB;
					String guideId = grip + "[" + targeti.getName() + "." + posA + "." + posB + "]." + guideReferenceNumber++;
					int g2 = posA + lenKmerA - 1;
					int g12 = posB + lenKmerB;
					if( this.guides.containsKey( guideId ) ) Utils.stop( guideId + " already in map", 0 );
					Guide newGuide = new Guide( guideId, "design", targeti, "", sirna, g2, g12, kmerA, kmerB, bridgeLen, t1, tlast, true );
					this.addGuide( newGuide );
				    }
				}
			    } // all posB
			} // all posA
		    } // end if found in target i
		} // all targets
	    } // all kmerB
	} // all kmerA
	this.removeIncompleteGrips(); // remove the grips and guides of partial grips, ie, that are not in all transcripts
	System.out.println( " done in " + ( ( System.currentTimeMillis() - startTime ) / 1000 ) + " seconds" );
    }

    public void crossHybridize() {
	System.out.print("cross-hybridizing... ");
	long startTime = System.currentTimeMillis();

	// Thread-safe queue to store new duplexes
	ConcurrentLinkedQueue<Guide> newDuplexes = new ConcurrentLinkedQueue<>();
	AtomicInteger guideReferenceNumber = new AtomicInteger();

	// Process each grip in parallel
	this.getGrips().parallelStream().forEach(grip -> {

		// For each guide in the current grip
		//System.out.println( "number of guides under grip " + grip + " = " + this.guideSetFromGrip(grip) );
		for (Guide g : this.guideSetFromGrip(grip)) {
		    Set<String> locations = new HashSet<>(); // Track locations within each grip processing context
		    locations.add( g.getCT().getName() + "." + g.getLocation() );
		    for (Guide other : this.guideSetFromGrip(grip) ) {
			//System.out.println( "locations: " + locations + ", " + other.getSequence() );
			// Skip self-comparison and
			String location = other.getCT().getName() + "." + other.getLocation();
			if( ( !locations.contains( location ) ) ) {
			    // Track location to prevent duplicate binding at the same location
			    locations.add(location);
			    // Generate a unique ID and create a new guide based on `g` and `other`
			    Guide newGuide = new Guide(other, g.getId(), "miRDesign", g.getSequence());
			    newDuplexes.add(newGuide); // Thread-safe addition to newDuplexes
			}
		    }
		}
	    });

	// Add the new guides from newDuplexes to the main data structure
	newDuplexes.forEach(this::addGuide);

	// Filter unbinders and remove incomplete sequences if needed
	//this.foldAndFilter(); // Uncomment if required
	this.removeIncompleteSequences();
       	System.out.println("done in " + ((System.currentTimeMillis() - startTime) / 1000) + " seconds");
    }
    public void crossHybridize2() { // hybridize the guides designed on single sites to all sites
	System.out.print( "cross-hybridizing... " );
	long startTime = System.currentTimeMillis(); // to measure time spent in crossHybridize
	Set<Guide> newDuplexes = new LinkedHashSet<>(); // to collect the new duplexes
	for( String grip : this.getGrips() ) // for all grips
	    for( Guide g : this.guideSetFromGrip( grip ) ) { // for all duplexes under each grip and all transcripts
		Set<String> locations = new HashSet<>();
		for( Guide other : this.guideSetFromGrip( g.getGrip() ) ) // try to bind the sequence of g to all other duplexes (with same grip)
		    if( g != other ) { // don't duplicate same guide
			if( !locations.contains( other.getLocation() ) ) { // don't bind to the same location twice
			    locations.add( other.getLocation() );
			    String id = g.getId() + "[" + other.getLocation() + "]" + guideReferenceNumber++;
			    Guide newGuide = new Guide( other, id, "miRDesign", g.getSequence() );
			    newDuplexes.add( newGuide ); // add new guide to the set of new guides
			}
		    }
	    }
	for( Guide d : newDuplexes )
	    this.addGuide( d ); // add new guides to the data structure (maps)
	//this.foldAndFilter(); // filter unbinders
	this.removeIncompleteSequences(); // remove sequences that don't bind all required targets
	System.out.println( "done in " + ( ( System.currentTimeMillis() - startTime ) / 1000 ) + " seconds" );
    }

    // Check that each sequence targets all transcripts in targets
    public void checkTargetsAndFilter() {
	Map<String,Set<Guide>> invalidDesigns = new HashMap<>();
	for( String sequence : this.sequences.keySet() ) { // for each sequence
	    Set<CodingTranscript> targets = new HashSet<>();
	    Set<Guide> guides = new HashSet<>();
	    for( Guide g : this.guideSetFromSequence( sequence ) ) {
		targets.add( g.getCT() );
		guides.add( g );
	    }
	    if( targets.size() < this.targets.size() ) {
		invalidDesigns.put( sequence, guides );
	    }
	}
	System.out.println( "checkTargetsAndFilter..." );
	for( String s : invalidDesigns.keySet() ) {
	    for( Guide g : invalidDesigns.get( s ) )
		this.removeGuide( g );
	    this.sequences.remove( s );
	}
	System.out.println( "done" );
    }

    // Safe fold and filter
    //     fold just to make sure a guide can bind all requested targets
    public void safeFoldAndFilter() {
	// for all sequences
	for( String guide : this.getSequences() ) {
	    // for all required target
	    
	}
	    
    }

    // Fold the guides in this.sets
    public void foldAndFilter() {
	System.out.print( "folding and filtering... " );
	long startTime = System.currentTimeMillis(); // to measure time spent in crossHybridize
	// Set<String> failedGuides = new HashSet<>();
	// for( String key : this.guides.keySet() ) { // for all guides
	//     this.guides.get( key ).fold();
	//     if( !checkCondition( this.guides.get( key ) ) ) failedGuides.add( key );
	int count = 0;
	for( Guide g : this.guideSet() ) // counting number of unfolded guides
	    if( !g.isFolded() ) count++;
	System.out.print( count + " duplexes... " );
	Set<Guide> failedGuides = this.guideSet().parallelStream()
	    .peek( guide -> guide.fold() )
	    .filter( guide -> !checkCondition( guide ) )
	    .collect( Collectors.toSet() );

	failedGuides.forEach( this::removeGuide );
	System.out.print( failedGuides.size() + " removed... " );
	System.out.println( "done in " + ( ( System.currentTimeMillis() - startTime ) / 1000 ) + " seconds" );
    }

    private boolean guideInRequired( Guide guide ) {
	return this.targets.contains( guide.getCT() );
    }

    private boolean checkCondition( Guide guide ) {
	//if( !guide.isFolded() ) Util.stop( "guide not folded:\n" + guide, 0 );
	//System.out.println( "checkCondition" );
	//boolean rep = guide.isFolded() && ( guide.bind() || ( guide.initiate() && guide.compensate() ) );
	//System.out.println( "rep: " + rep );
	// recheck for brigde lenght, which may change after folding
	return guide.isFolded() && guide.getBridgeLength() <= this.distance && ( guide.seedBind() || ( guide.suppBind() ) );
    }

    // check no U at the end, absence of exclusions, and accept best GC contents only
    private boolean checkSequence( String guide ) {
	return
	    guide.charAt( guide.length() - 1 ) != 'U' &&
	    !containsExclusions( guide, this.exclusions ) &&
	    containsRightGCPercentage( guide, this.gcPercentMin, this.gcPercentMax );
    }

    private Set<String> generateAllExtensions( String sirna ) {
	Set<String> pool = new HashSet<>();
	List<String> seeds = new ArrayList<>(); // get kmers for seed extensions
	int k = this.kmerMap3p.getK(); // get the length of the seed prefix
     	if( k < 7 ) // else there is nothing to extend
	    StringSequence.generateKMers( 7 - k, this.kmerMap3p.getExclusions(), seeds ); // k = seed length (7) - prefix length (k)
	List<String> g12 = new ArrayList<>(); // get kmers for extensions
        k = 1;
	StringSequence.generateKMers( k, this.kmerMap5p.getExclusions(), g12 ); // k = 1
	List<String> supps = new ArrayList<>(); // get kmers for supp extensions
	k = this.kmerMap5p.getK(); // get the length of the supp region prefix
     	if( k < 5 ) // else there is nothing to extend
	    StringSequence.generateKMers( 5 - k, this.kmerMap5p.getExclusions(), supps ); // k = supp region length (5) - prefix length (k)

	for( String seed : seeds )
	    for( String g : g12 )
		for( String supp : supps )
		    pool.add( sirna.substring( 0, this.kmerMap3p.getK() + 1 ) + // g1+seed prefix
			      seed + // seed variable part
			      sirna.substring( 8, 11 ) + // box A
			      g + // g12
			      sirna.substring( 12, 12 + this.kmerMap5p.getK() ) + // supp prefix
			      supp + // supp variable part
			      sirna.substring( 17, sirna.length() ) );
	return pool;
    }

    // Extends the seeds when seed prefixes are under 7 nts
    public void extendSeeds( int seedPrefixLength ) {
	System.out.print( "mutating the seed... " );
	List<String> kMers = new ArrayList<>(); // get kmers for extensions
	int k = seedPrefixLength; // get the length of the seed prefix
     	if( k < 7 ) // else there is nothing to extend
	    StringSequence.generateKMers( 7 - k, this.exclusions, kMers ); // k = seed length (7) - prefix length (k)

	for( Guide guide : this.guideSet() ) { // for all guides
	    int index = guide.getId().lastIndexOf( "]." );
	    String guideIdRoot = guide.getId().substring( 0, index + 2 );
	    for( String kmer : kMers ) {
		// Extract the part before the closing bracket and append the incremented reference number
		String guideId = guideIdRoot + guideReferenceNumber++;
		if( guide.getSequence().indexOf( kmer, k + 1 ) != (k + 1) ) { // don't redo existing guide
		    Guide newDesign = guide.mutate( k + 1, kmer, guideId );
		    if( checkSequence( newDesign.getSequence() ) ) {
			this.addGuide( newDesign );
		    }
		}
	    }
	}
	System.out.println( " done" );
	//this.foldAndFilter(); // filter unbinders
    }

    // Extends g12
    public void extendG12() {
	System.out.print( "mutating g12... " );
	List<String> kMers = new ArrayList<>(); // get kmers for extensions
	int k = 1;
	StringSequence.generateKMers( 1, this.exclusions, kMers ); // k = 1

	for( Guide guide : this.guideSet() ) {
	    int index = guide.getId().lastIndexOf( "]." );
	    String guideIdRoot = guide.getId().substring( 0, index + 2 );
	    for( String kmer : kMers ) {
		// Extract the part before the closing bracket and append the incremented reference number
		String guideId = guideIdRoot + guideReferenceNumber++;
		if( guide.getSequence().charAt( 11 ) != kmer.charAt( 0 ) ) { // don't redo existing g12
		    Guide newDesign = guide.mutate( 11, kmer, guideId );
		    if( checkSequence( newDesign.getSequence() ) ) {
			this.addGuide( newDesign );
		    }
		}
	    }
	}
	System.out.println( " done" );
	//this.foldAndFilter(); // filter unbinders
    }

    public void extendSupp( int suppPrefixLength ) {
	System.out.print( "mutating the supplementary region... " );
	List<String> kMers = new ArrayList<>(); // get kmers for extensions
	int k = suppPrefixLength; // get the length of the supp region prefix
     	if( k < 5 ) // else there is nothing to extend
	    StringSequence.generateKMers( 5 - k, this.exclusions, kMers ); // k = supp region length (5) - prefix length (k)

	for( Guide guide : this.guideSet() ) { // for all guides
	    // extract part before the "]." and create guideIdRoot
	    int index = guide.getId().lastIndexOf( "]." );
	    String guideIdRoot = guide.getId().substring( 0, index + 2 );
	    for( String kmer : kMers ) {
		if( guide.getSequence().indexOf( kmer, k + 12 ) != k + 12 ) { // don't redo existing supp
		    String newSequence = guide.getSequence( 0, k+12 ) + kmer + guide.getSequence( k + 12 + kmer.length(), guide.getLength() );
		    if( checkSequence( newSequence ) ) {
			Guide newDesign = guide.mutate( newSequence, guideIdRoot + guideReferenceNumber++ );
			this.addGuide( newDesign );
		    }
		}
	    }
	}
	System.out.println( " done" );
	//this.foldAndFilter(); // filter unbinders
    }

    // // Extends A-Box
    // public void extendABox() {
    // 	List<String> kMers = new ArrayList<>(); // get kmers for extensions
    // 	int k = 3; // A-Box length
    // 	int g9 = 8; // index of g9
    // 	Utilities.generateKMers( 3, this.kmerMap.getExclusions(), kMers );
    // 	//System.out.println( kMers.size() + " A-Box extensions to be considered" );
    // 	for( Guide guide : this.guideSet() ) {
    // 	    int id = 0;
    // 	    for( String kmer : kMers ) {
    // 		String guideId = guide.getId().replaceAll( "\\]\\.(\\d+)$", "]." + guideReferenceNumber++;
    // 		if( guide.getSequence().indexOf( kmer, g9 ) != g9 ) { // don't redo existing guide
    // 		    Guide newDesign = guide.mutate( g9, kmer, guideId );
    // 		    if( checkSequence( newDesign.getSequence(), this.exclusions ) ) {
    // 			this.addGuide( newDesign );
    // 		    }
    // 		}
    // 	    }
    // 	}
    // }

    // // Extends D-Box
    // public void extendDBox() {
    // 	List<String> kMers = new ArrayList<>(); // get kmers for extensions
    // 	int g18 = 17;
    // 	int k = this.guideSize - g18; // D-Box length
    // 	Utilities.generateKMers( k, this.kmerMap.getExclusions(), kMers );
    // 	//System.out.println( kMers.size() + " D-Box extensions to be considered" );
    // 	for( Guide guide : this.guideSet() ) {
    // 	    int id = 0;
    // 	    for( String kmer : kMers ) {
    // 		String guideId = guide.getId().replaceAll( "\\]\\.(\\d+)$", "]." + guideReferenceNumber++;
    // 		if( guide.getSequence().indexOf( kmer, g18 ) != g18 ) { // don't redo existing guide
    // 		    Guide newDesign = guide.mutate( g18, kmer, guideId );
    // 		    if( checkSequence( newDesign.getSequence(), this.exclusions ) ) {
    // 			this.addGuide( newDesign );
    // 		    }
    // 		}
    // 	    }
    // 	}
    // }

    //    public void

    
    /**
       design siRNA for sequence:
       5'-AAUGUAAAUUUGUGUUUAUUU-3' target
       3'-UUACAUUUAAACACAAAUAAA-5' shRNA
       5'-             GUU     -3' MRE (GUU)
       3'-             ACG     -5' kMer (GCA reverse by countBPs)
       3'-             CAA     -5' currentExtension
       5'-AAAUAAACACAAAUUUACAUU-3'
     */

    // print this.matches, that is the guides by sequence
    public String toString() {
	String res = "";
	for( String design : this.getSequences() ) {
	    //System.out.println( "#guides under design " + design + " = " + this.sequences.get( design ).size() + ", expecting more or equal to " + this.targets.size() );
	    if( this.sequences.get( design ).size() >= this.targets.size() ) {
		res += design + "\n";
		for( Guide g : this.sequences.get( design ) )
		    res += g;
	    }
	}
	return res;
    }

    // KMerMap kmerMap;
    // Set<Integer> targets;
    // List<ProteinCodingTranscript> excludedSubset;
    // Map<String,Guide> guides = new HashMap<>(); // a map for <guide name,guide info> of the guides
    // Map<String,Set<Guide>> grips = new HashMap<>(); // a map to keep the guides under their kmer prefix grips: <seedPrefix/suppPrefix,Set<guide>>
    // Map<String,Set<Guide>> matches = new HashMap<>(); // a map to keep the guide sequences linked to their set of guides: <guide sequence,set of guides>
    // int distance = 0; // distance constraint between sites in maps A and B
    // int guideSize;
    // List<String> exclusions;
    // boolean siRNAInAllTargets = false;
    // static int guideReferenceNumber = 0; // global guide reference number to avoid guide key duplicates
    // private double gcPercentMin = 0.3;
    // private double gcPercentMax = 0.64;
    // private int guideAdded = 0;
    // private int guideRemoved = 0;
    
    public void dumpAll() {
	System.out.println( "dumpAll:" );
	// dump simple variables
	System.out.println( "distance: " + this.distance );
	System.out.println( "guideSize: " + this.guideSize );
	System.out.println( "siRNAInAllTargets: " + this.siRNAInAllTargets );
	System.out.println( "gcPercentMin: " + this.gcPercentMin );
	System.out.println( "gcPercentMax: " + this.gcPercentMax );
	System.out.println( "guideAdded: " + this.guideAdded );
	System.out.println( "guideRemoved: " + this.guideRemoved );
	System.out.println( "guideReferenceNumber: " + guideReferenceNumber );

	// dump lists and sets
	System.out.println( "targets: " + this.targets );
	System.out.println( "excludedSubset: " + this.excludedSubset );
	System.out.println( "exclusions: " + this.exclusions );

	// dump maps
	System.out.println( "KMerMap:\n" + this.kmerMap3p );
	System.out.println( "KMerMap:\n" + this.kmerMap5p );
	System.out.println( "guides:\n" + this.guides );
	System.out.println( "grips:\n" + this.grips );
	System.out.println( "matches:\n" + this.sequences );
    }
}
