package ca.iric.major.common;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Utility for generating miRNA seed candidates (RIMap-RISC rules) for a given 8-mer target.
 *
 * Rules:
 *
 */
public final class GenSeed {

    // guide index
    private final static int g2 = 0;
    private final static int g3 = 1;
    private final static int g4 = 2;
    private final static int g5 = 3;
    private final static int g6 = 4;
    private final static int g7 = 5;
    private final static int g8 = 6;
    // target index
    private final static int t2 = 7;
    private final static int t3 = 6;
    private final static int t4 = 5;
    private final static int t5 = 4;
    private final static int t6 = 3;
    private final static int t7 = 2;
    private final static int t8 = 1;
    private final static int t9 = 0;

    private GenSeed() {
        // utility class
    }

    /**
     * Generate all seed candidates given allowed nts at each positions.
     *
     * @param allowed array of arrays of nts
     * @param seeds list to which the seeds are added
     */
    private static void addSeeds( char[][] allowed, Set<String> seeds ) {
	if( allowed == null || allowed.length == 0 ) {
	    return;
	}

	// Sanity check: no empty positions
	for( int i = 0; i < allowed.length; i++ ) {
	    if( allowed[i] == null || allowed[i].length == 0 ) {
		throw new IllegalArgumentException( "allowed[" + i + "] is empty" );
	    }
	}
	StringBuilder sb = new StringBuilder( allowed.length );
	generate( allowed, 0, sb, seeds );
    }

    private static void generate( char[][] allowed, int pos, StringBuilder sb, Set<String> seeds) {
	if( pos == allowed.length ) {
	    // Reached full length => record one seed
	    seeds.add( sb.toString() );
	    return;
	}

	for( char c : allowed[pos] ) {
	    sb.append( c );
	    generate( allowed, pos + 1, sb, seeds );
	    sb.setLength( sb.length() - 1 ); // backtrack
	}
    }

    /**
     * Generate all seed candidates for an 8-mer according to the RIMap-RISC.
     *
     * @param kmer8 8-nt RNA sequence (A,C,G,U)
     * @return list of unique valid 7-nt seed strings
     */
    public static Set<String> generateSeedsFor8mer( String kmer8 ) {

	Set<String> seeds = new LinkedHashSet<>();
	char[][] allowed = new char[7][];
	
	// 4mer g2-g5 + sample g6-g8
	allowed[g2] = StringSequence.allowedNts( kmer8.charAt( t2 ), true ); // g2 = complement of t2
	allowed[g3] = StringSequence.allowedNts( kmer8.charAt( t3 ), true ); // g3 = complement of t3
	allowed[g4] = StringSequence.allowedNts( kmer8.charAt( t4 ), true ); // g4 = complement of t4
	allowed[g5] = StringSequence.allowedNts( kmer8.charAt( t5 ), true ); // g5 = complement of t5
	allowed[g6] = StringSequence.nucleotides; // sample g6
	allowed[g7] = StringSequence.nucleotides; // sample g7
	allowed[g8] = StringSequence.nucleotides; // sample g8
	addSeeds( allowed, seeds );

	// 4mer g3-g6 + sample g2,g7-g8
	allowed[g2] = StringSequence.nucleotides; // sample g2
	allowed[g6] = StringSequence.allowedNts( kmer8.charAt( t6 ), true ); // g6 = complement of t6
	addSeeds( allowed, seeds );

	// 4mer g4-g7 + sample g2-g3,g8
	allowed[g3] = StringSequence.nucleotides; // sample g3
	allowed[g7] = StringSequence.allowedNts( kmer8.charAt( t7 ), true ); // g7 = complement of t7
	addSeeds( allowed, seeds );

	// 7mer + sample g5
	allowed[g2] = StringSequence.allowedNts( kmer8.charAt( t2 ), true ); // g2 = complement of t4
	allowed[g3] = StringSequence.allowedNts( kmer8.charAt( t3 ), true ); // g3 = complement of t4
	allowed[g4] = StringSequence.allowedNts( kmer8.charAt( t4 ), true ); // g4 = complement of t4
	allowed[g5] = StringSequence.nucleotides; // sample g5
	allowed[g6] = StringSequence.allowedNts( kmer8.charAt( t6 ), true ); // g6 = complement of t6
	allowed[g7] = StringSequence.allowedNts( kmer8.charAt( t7 ), true ); // g7 = complement of t7
	allowed[g8] = StringSequence.allowedNts( kmer8.charAt( t8 ), true ); // g8 = complement of t8
	addSeeds( allowed, seeds );

	// 7mer + 1 deletion
	allowed[g2] = StringSequence.allowedNts( kmer8.charAt( t2 ), true ); // g2 = complement of t2
	allowed[g3] = StringSequence.allowedNts( kmer8.charAt( t3 ), true ); // g3 = complement of t3
	allowed[g4] = StringSequence.nucleotides; // sample g4
	allowed[g5] = StringSequence.allowedNts( kmer8.charAt( t4 ), true ); // g5 = complement of t4
	allowed[g6] = StringSequence.allowedNts( kmer8.charAt( t5 ), true ); // g6 = complement of t5
	allowed[g7] = StringSequence.allowedNts( kmer8.charAt( t6 ), true ); // g7 = complement of t6
	allowed[g8] = StringSequence.allowedNts( kmer8.charAt( t7 ), true ); // g8 = complement of t7
	addSeeds( allowed, seeds );

	// 7mer + 1 deletion
	allowed[g4] = StringSequence.allowedNts( kmer8.charAt( t4 ), true ); // g4 = complement of t4
	allowed[g5] = StringSequence.nucleotides; // sample g5
	addSeeds( allowed, seeds );

	// 7mer + 1 insertion
	allowed[g2] = StringSequence.allowedNts( kmer8.charAt( t2 ), true ); // g2 = complement of t2
	allowed[g3] = StringSequence.allowedNts( kmer8.charAt( t3 ), true ); // g3 = complement of t3
	allowed[g4] = StringSequence.allowedNts( kmer8.charAt( t4 ), true ); // g4 = complement of t4
	allowed[g5] = StringSequence.allowedNts( kmer8.charAt( t6 ), true ); // g5 = complement of t6
	allowed[g6] = StringSequence.allowedNts( kmer8.charAt( t7 ), true ); // g6 = complement of t7
	allowed[g7] = StringSequence.allowedNts( kmer8.charAt( t8 ), true ); // g7 = complement of t8
	allowed[g8] = StringSequence.allowedNts( kmer8.charAt( t9 ), true ); // g8 = complement of t9
	addSeeds( allowed, seeds );

	return seeds;
    }

}
