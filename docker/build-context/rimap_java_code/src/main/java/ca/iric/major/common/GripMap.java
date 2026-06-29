/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * GripMap maintains the common kmer pairs for a list of required, optional and excluded transcripts
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 1.0 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.util.Map;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Set;
import java.util.HashSet;

public class GripMap {

    private static int SEEDLENGTH = 7;
    private static int BOXALENGTH = 3;
    private static Set<Grip> EMPTYSET = new HashSet<>();

    private boolean gripIsInAllTranscriptsOfMap( Map<String,Set<Grip>> map, String grip, List<ProteinCodingTranscript> pcts ) {
	Set<ProteinCodingTranscript> present = new HashSet<>();
	for( Grip g : map.get( grip ) )
	    present.add( g.getPCT() );
	if( present.containsAll( pcts ) ) return true;
	return false;
    }

    private void buildMap( Map<String,Set<Grip>> map, List<ProteinCodingTranscript> pcts, Set<String> constraintSet, boolean mustBeInAllTranscripts ) {
	Set<ProteinCodingTranscript> setOfTranscripts = new HashSet<>( pcts );
	KMerMap map3p = new KMerMap( pcts, this.region, this.exclusions, this.k3 );
    	KMerMap map5p = new KMerMap( pcts, this.region, this.exclusions, this.k5 );
	//System.out.println( "map 3':\n" + map3p );
	//System.out.println( "map 5':\n" + map5p );
	for( String kmer3p : map3p.getExclusive( setOfTranscripts ) ) {
	    for( String kmer5p : map5p.getExclusive( setOfTranscripts ) ) {
		if( constraintSet.isEmpty() || constraintSet.contains( kmer3p + "/" + kmer5p ) ) { // check only if in constraintSet, as those of the required transcripts for optional and excluded
		    for( ProteinCodingTranscript pct : pcts ) { // for each transcript
			//System.out.println( "for pct: " + pct );
			List<Integer> positions3p = map3p.getPositions( kmer3p, pct ); // get positions of 3' kmer in transcript
			List<Integer> positions5p = map5p.getPositions( kmer5p, pct ); // get positions of 5' kmer in transcript
			if( !positions3p.isEmpty() && !positions5p.isEmpty() ) {
			    for( Integer pos3p : positions3p ) { // for each position in 3' kmer map
				//System.out.println( "pos3p: " + pos3p );
				for( Integer pos5p : positions5p ) { // for each position in 5' kmer map
				    //System.out.println( "pos3p: " + pos3p );
				    int g8Tog12Distance = ( pos3p - ( SEEDLENGTH - this.k3 ) ) - ( pos5p + this.k5 ) - 1;
				    if( g8Tog12Distance >= BOXALENGTH && g8Tog12Distance <= this.maxDistance )  {
					//System.out.println( g8Tog12Distance + " is between " + BOXALENGTH + " and " + this.maxDistance );
					// add grip to map
					String key = kmer3p + "/" + kmer5p;
					Grip grip = new Grip( pct, pos3p, pos5p, g8Tog12Distance );
					Set<Grip> gripsForThisKey = map.computeIfAbsent( key, k -> new HashSet<>() );
					gripsForThisKey.add( grip );
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	//System.out.println( "map: " + map );
	if( mustBeInAllTranscripts ) {
	    // keep only the common grips
	    // make a set of invalid grips
	    Set<String> invalidGrips = new HashSet<>();
	    for( String grip : map.keySet() )
		if( !gripIsInAllTranscriptsOfMap( map, grip, pcts ) )
		    invalidGrips.add( grip );
	    // remove the invalid grips
	    for( String grip : invalidGrips )
		map.remove( grip );
	}
    }

    Map<String,Set<Grip>> requiredGripMap = new LinkedHashMap<>(); // a map to keep the grips, keys = "3'kmer/5'kmer"
    Map<String,Set<Grip>> optionalGripMap = new LinkedHashMap<>(); // a map to keep the grips, keys = "3'kmer/5'kmer"
    Map<String,Set<Grip>> excludedGripMap = new LinkedHashMap<>(); // a map to keep the grips, keys = "3'kmer/5'kmer"
    List<ProteinCodingTranscript> required;
    List<ProteinCodingTranscript> optional;
    List<ProteinCodingTranscript> excluded;
    int region;
    int maxDistance; // max number of nucleotides between 5' and 3' kmers
    List<String> exclusions; // kmers to be excluded
    int k3; // length of the 3' kmers
    int k5; // length of the 5' kmers

    public GripMap(
		   List<ProteinCodingTranscript> required,
		   List<ProteinCodingTranscript> optional,
		   List<ProteinCodingTranscript> excluded,
		   int region,
		   int maxDistance,
		   List<String> exclusions,
		   int k3,
		   int k5 ) {
	this.required = required;
	this.optional = optional;
	this.excluded = excluded;
	this.region = region;
	this.maxDistance = maxDistance;
	this.exclusions = exclusions;
	this.k3 = k3;
	this.k5 = k5;
	// build the grip map for the required transcripts
	this.buildMap( this.requiredGripMap, this.required, new HashSet<String>(), true );
	this.buildMap( this.optionalGripMap, this.optional, this.requiredGripMap.keySet(), false );
	this.buildMap( this.excludedGripMap, this.excluded, this.requiredGripMap.keySet(), false );
	//System.out.println( this.requiredGripMap );
	//System.exit( 0 );
    }

    // get 3' and 5' kmers from a grip map key
    public String get3pKey( String key ) { return key.substring( 0, this.k3 ); }
    public String get5pKey( String key ) { return key.substring( this.k3 + 1, key.length() ); }
    public int getK5() { return this.k5; }
    public int getK3() { return this.k3; }
    public int getMaxDistance() { return this.maxDistance; }
    public int getRegion() { return this.region; }
    public List<String> getExclusions() { return this.exclusions; }
    public List<ProteinCodingTranscript> getRequired() { return this.required; }
    public Map<String,Set<Grip>> getRequiredGripMap() { return this.requiredGripMap; }
    public Map<String,Set<Grip>> getOptionalGripMap() { return this.optionalGripMap; }
    public Map<String,Set<Grip>> getExcludedGripMap() { return this.excludedGripMap; }
    public Set<Grip> getRequiredGrips( String grip ) {
	if( this.requiredGripMap.get( grip ) == null ) return EMPTYSET;
	else return this.requiredGripMap.get( grip );
    }
    public Set<Grip> getOptionalGrips( String grip ) {
	if( this.optionalGripMap.get( grip ) == null ) return EMPTYSET;
	else return this.optionalGripMap.get( grip );
    }
    public Set<Grip> getExcludedGrips( String grip ) {
	if( this.excludedGripMap.get( grip ) == null ) return EMPTYSET;
	else return this.excludedGripMap.get( grip );
    }

    @Override
    public String toString() {
	String ret = "";
	ret += "Required " + this.requiredGripMap.keySet().size() + " grips:\n"; // + this.requiredGripMap.keySet().toString() + "\n";
	// for( String grip : this.requiredGripMap.keySet() )
	//     ret += grip + ": " + this.requiredGripMap.get( grip ) + "\n";
	ret += "Optional " + this.optionalGripMap.keySet().size() + " grips:\n"; // + this.optionalGripMap.keySet().toString() + "\n";
	// for( String grip : this.optionalGripMap.keySet() )
	//     ret += grip + ": " + this.optionalGripMap.get( grip ) + "\n";
	ret += "Excluded " + this.excludedGripMap.keySet().size(); // + " grips:\n" + this.excludedGripMap.keySet().toString() + "\n";
	// for( String grip : this.excludedGripMap.keySet() )
	//     ret += grip + ": " + this.excludedGripMap.get( grip ) + "\n";
	return ret;
    }
}
