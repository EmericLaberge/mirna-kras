/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import ca.iric.major.common.*;

import java.util.*;
import java.util.stream.Collectors;
import java.util.Arrays;
//import java.util.List;
//import java.util.ArrayList;
//import java.util.Set;
//import java.util.HashSet;

import java.io.*;
import java.nio.file.*;
import java.nio.file.Files;
import java.nio.file.Path;

import java.text.DecimalFormat;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.core.type.TypeReference;

import org.apache.commons.csv.*;

public class RISCScan {

    public static String dropAfterDot( String input ) {
        int dotIndex = input.indexOf('.');
        if( dotIndex != -1 ) {
            return input.substring( 0, dotIndex );
        } else {
            return input; // No dot found, return the original string
        }
    }

    public static void writeStringToFile( String filePath, String content ) throws IOException {
        Files.write( Paths.get( filePath ), content.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING );
    }

    // Function to read sequence from bedTool FASTA (in DNA) and return the RNA version
    public static void convertSequenceToUpperCaseRNA( String inputFile ) throws IOException {
	String RNASequence = "";
        List<String> lines = Files.readAllLines( Paths.get( inputFile ) );
	for( String line : lines ) {
	    if( line.startsWith( ">" ) ) RNASequence += line + "\n";
	    else {
		String modifiedLine = line.toUpperCase().replace( 'T', 'U' );
		RNASequence += modifiedLine;
	    }
        }
	// replace the sequence in inputFile
	try {
            writeStringToFile( inputFile, RNASequence );
            //System.out.println( "File written successfully." );
        } catch( IOException e ) {
            e.printStackTrace();
        }
    }

    // Function to extract genomic location from GTF file
    public static Map<String, String> extractGenomicLocation( String gtfFile, String enstId ) throws IOException {
        Reader in = new FileReader( gtfFile );
        Iterable<CSVRecord> records = CSVFormat.TDF
	    .withCommentMarker('#')
	    .withHeader( "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute" )
	    .parse(in);

        for (CSVRecord record : records) {
            String feature = record.get( "feature" );
            String attribute = record.get( "attribute" );

            if (feature.equals( "transcript" ) && attribute.contains( "transcript_id \"" + enstId + "\"" ) ) {
                String seqname = record.get( "seqname" );
                String start = record.get( "start" );
                String end = record.get( "end" );
                String strand = record.get( "strand" );

                Map<String, String> location = new HashMap<>();
                location.put( "seqname", seqname );
                location.put( "start", start );
                location.put( "end", end );
                location.put( "strand", strand );
                return location;
            }
        }

        System.out.println( "ENST ID " + enstId + " not found in GTF file." );
        return null;
    }

    // Function to fetch sequence using bedtools
    public static void fetchSequence( String pctId, String miRId, String dataPath, String fastaFile, String seqname, int start, int end, String outputFile ) throws IOException, InterruptedException {
	//System.out.println( "fetchSequence( " + fastaFile + ", " + seqname + ", " + start + ", " + end + ", " + outputFile + " )" );
        String bedContent = seqname + "\t" + (start - 1) + "\t" + end + "\n";
        Path bedFile = Paths.get( dataPath + "region.bed" + "." + pctId + "." + miRId );

        // Write the BED content to a file
        Files.write( bedFile, bedContent.getBytes() );

        // Run bedtools getfasta
        ProcessBuilder pb = new ProcessBuilder( "bedtools", "getfasta", "-fi", fastaFile, "-bed", bedFile.toString(), "-fo", outputFile );
        pb.inheritIO();
        Process process = pb.start();
        int exitCode = process.waitFor();

        if (exitCode == 0) {
            //System.out.println("Sequence extraction completed. Output is in " + outputFile);
        } else {
	    Utils.stop( "Error running bedtools.", 5 );
            //System.out.println("Error running bedtools.");
        }

	// delete the region.bed file
	try {
	    Files.deleteIfExists( bedFile );
	} catch (IOException e) {
	    System.err.println( "Could not delete file: " + bedFile );
	    e.printStackTrace();
	}
    }

    // Function to find 31mer in the extracted sequence and adjust coordinates
    public static Map<String, Integer> search31merAndAdjustCoordinates( String sequence, int start, String query31mer ) throws IOException {
        int index = sequence.indexOf( query31mer );
        if( index == -1 ) {
            //Utils.debug("31mer not found in the extracted sequence: " + query31mer );
            return null;
        }
        int adjustedStart = start + index;
        int adjustedEnd = adjustedStart + query31mer.length();
        Map<String, Integer> coordinates = new HashMap<>();
        coordinates.put( "start", adjustedStart );
        coordinates.put( "end", adjustedEnd );
        return coordinates;
    }

    public static int findLongestStretch(double[] array, double threshold) {
        int longest = 0;
        int current = 0;
        for(double value : array) {
            if( value >= threshold ) {
                current++; // Increment the current stretch
                longest = Math.max( longest, current ); // Update the longest stretch if needed
            } else {
                current = 0; // Reset the current stretch
            }
        }
        return longest;
    }

    /**
     * length-weighted mean bedGraph value over the half-open range [qStart, qEnd).
     * returns -1f if the range has no coverage in the bedGraph.
     *
     * Assumes:
     *  - starts[], ends[] are sorted by starts ascending
     *  - intervals are 0-based, half-open [start, end)
     *  - vals[] has same length as starts[]/ends[]
     */
    static float fetchPhastConsScoresPrimitive(int[] starts, int[] ends, float[] vals,
					       int qStart, int qEnd) {
	if (qStart >= qEnd) return -1f;

	final int n = starts.length;
	if (n == 0) return -1f;

	// Find first candidate interval: start >= qStart
	int idx = Arrays.binarySearch(starts, qStart);
	int i = (idx >= 0) ? idx : (-idx - 1);
	// Also check the one just before, in case qStart falls inside it
	if (i > 0) i--;

	double wSum = 0.0;   // weighted sum of values
	long   wTot = 0L;    // total overlap length

	for (; i < n && starts[i] < qEnd; i++) {
	    // overlap between [starts[i], ends[i]) and [qStart, qEnd)
	    int s = Math.max(starts[i], qStart);
	    int e = Math.min(ends[i],   qEnd);
	    int overlap = e - s;
	    if (overlap > 0) {
		wSum += (double) overlap * (double) vals[i];
		wTot += overlap;
	    }
	}

	return (wTot > 0) ? (float) (wSum / wTot) : -1f;
    }
    
    // Function to fetch PhastCons scores from BedGraphLine array
    // public static double fetchPhastConsScoresOld( BedGraphLine[] bdLines, int start, int end ) throws IOException {
    //     double[] scores = new double[end - start];
    //     int index = getBedGraphStartIndex( bdLines, start, 0, bdLines.length - 1 );
    // 	do {
    // 	    int regionStart = bdLines[index].getStart();
    // 	    int regionEnd = bdLines[index].getEnd();
    // 	    double score = bdLines[index].getConservation();
    // 	    if( regionStart >= end )
    // 		break;  // Exit early if the current region is beyond the end position
    // 	    if( regionEnd <= start )
    // 		continue;  // Skip the current line if it is before the start position
    // 	    int effectiveStart = Math.max( regionStart, start );
    // 	    int effectiveEnd = Math.min( regionEnd, end );
    // 	    for( int i = effectiveStart; i < effectiveEnd; i++ )
    // 		scores[i - start] = score;
    // 	    index++;
    // 	} while( index < bdLines.length );
    // 	// Compute and print PhastCons scores
    // 	//double minScore = Arrays.stream(scores).min().orElse(0.0);
    // 	//double maxScore = Arrays.stream(scores).max().orElse(0.0);
    // 	//double averageScore = Arrays.stream(scores).average().orElse(0.0);
    // 	//System.out.println( "longest stretch: " + findLongestStretch( scores, 1.0 ) );
    // 	return Arrays.stream(scores).sum();
    // 	//int lenScores = scores.length;

    // 	//System.out.println("PhastCons scores for " + seqname + ":" + start + "-" + end + ":");
    // 	//System.out.println("Min: " + minScore + ", Average: " + averageScore + ", Max: " + maxScore + ", Sum: " + sumScores + ", Count: " + lenScores);
    //     //return sumScores;
    // }

    static float valueAt( int pos, int[] starts, int[] ends, float[] vals, float missing ) {
	// binary search in starts[]
	int idx = Arrays.binarySearch( starts, pos );
	if( idx >= 0 ) {
        // pos equals an interval start — check it is inside [start, end)
        return( pos < ends[idx] ) ? vals[idx] : missing;
	} else {
	    // insertion point - 1 is the interval that might contain pos
	    int i = -idx - 2;         // i = insertionPoint - 1
	    if( i >= 0 && pos < ends[i] ) {
		return vals[i];
	    }
	    return missing;
	}
    }

    private static void arrayizeBedGraphLines( List<String> bdEntries, BedGraphLine[] bdLines ) {
	int i = 0;
	for( String line: bdEntries ) {
	    String[] parts = line.split( "\t" );
	    bdLines[i++] = new BedGraphLine( parts[0], Integer.parseInt( parts[1] ), Integer.parseInt( parts[2] ), Double.parseDouble( parts[3] ) );
	}
    }

    private static int getBedGraphStartIndex( BedGraphLine[] bdLines, int start, int left, int right ) {
	while( left <= right ) {
            int mid = left + (right - left) / 2;  // Prevent potential overflow
            if( bdLines[mid].getStart() == start ) {
                return mid;  // Found the exact match
            } else if( bdLines[mid].getStart() < start ) {
                left = mid + 1;  // Search in the right half
            } else {
                right = mid - 1;  // Search in the left half
            }
        }
        // If the exact match is not found, return the insertion point
        return left == 0 ? left : left - 1;
    }

    private static boolean regionMatches( String region, boolean UTR5, boolean CDS, boolean UTR3 ) {
        switch( region ) {
	case "5UTR":
	    return UTR5;
	case "CDS":
	    return CDS;
	case "3UTR":
	    return UTR3;
	default: // catch the overapping cases
	    if( region.contains( "5UTR" ) && UTR5 ) {
		return true;
	    }
	    if( region.contains( "CDS" ) && CDS ) {
		return true;
	    }
	    if( region.contains( "3UTR" ) && UTR3 ) {
		return true;
	    }
	    return false; // Return false otherwise
        }
    }

    private static boolean accessMatches( double seedThreshold, double seedAccessibility, double suppThreshold, double suppAccessibility ) {
	seedAccessibility += 0.005;
	suppAccessibility += 0.005;
	//System.out.println( "seedAccessible: " + seedAccessible + ", seed accessibilit: " + seedAccessibility  );
	return seedAccessibility >= seedThreshold && suppAccessibility >= suppThreshold;
    }

    private static boolean suppMatches( String suppType, boolean supplementaryPaired ) {
        switch( suppType ) {
	case "suppless":
	case "supplementless":
	    return !supplementaryPaired;
	default:
	    return true; // Return true for other values
        }
    }

    private static boolean conservedMatches( int conserved, double phastconsScore ) {
	return conserved <= phastconsScore || phastconsScore == -1.0;
    }

    private static int printInteraction(
					List<JsonInteraction> toBePrinted,
					JsonInteraction jInt,
					String minSeed,
					double maxKd,
					boolean UTR5,
					boolean CDS,
					boolean UTR3,
					boolean supplementaryPaired,
					double seedAccessible,
					double supplementaryAccessible,
					int conserved ) {
	//System.out.println( "printInteraction( " + jInt.getPositionEnd() + ", " + minSeed + ", " + maxKd + ", " + UTR5 + ", " + CDS + ", " + UTR3 + ", " + supplementaryPaired + ", " + seedAccessible + ", " + supplementaryAccessible + ", " + conserved + " )" );
	String jIntSeedType = jInt.getSeedType();
	String jIntRegion = jInt.getRegion();
	String suppType = jInt.getSupplementaryType();
	boolean seedLess = jIntSeedType.equals( "seedless" );
	double kD = jInt.getKd();
	double seedAccessibility = jInt.getSeedAccessibility(); // seedAccessible => seedAccessibility >= 0.15
	double suppAccessibility = jInt.getSuppAccessibility(); // supplementaryAccessible => suppAccessibility >= 0.15
	double phastconsScore = jInt.getPhastconsScore(); // conserved => phastcons score >= 21.7; (70% of the target site nts are conserved)
	//System.out.println( "seedAccessible: " + seedAccessible + " seedAccessibility: " + seedAccessibility );
	//System.out.println( "supplementaryAccessible: " + supplementaryAccessible + " suppAccessibility: " + suppAccessibility );
	//System.out.println( "tentative: " + jInt.getPositionStart() );
	if( kD <= maxKd &&
	    minSeed.equals( "seedless" ) &&
	    regionMatches( jIntRegion, UTR5, CDS, UTR3 ) &&
	    suppMatches( suppType, supplementaryPaired ) &&
	    accessMatches( seedAccessible, seedAccessibility, supplementaryAccessible, suppAccessibility ) &&
	    conservedMatches( conserved, phastconsScore ) ) {
	    //System.out.println( "added: " + jInt.getPositionEnd() );
	    toBePrinted.add( jInt );
	    return 1;
	}
	else {
	    // extract core seed type
	    int pos = jInt.getSeedType().indexOf( "mer" );
	    if( kD <= maxKd &&
		!seedLess && minSeed.compareTo( jInt.getSeedType().substring( pos-1, pos+3 ) ) <=0 &&
		regionMatches( jIntRegion, UTR5, CDS, UTR3 ) &&
		suppMatches( suppType, supplementaryPaired ) &&
		accessMatches( seedAccessible, seedAccessibility, supplementaryAccessible, suppAccessibility ) &&
		conservedMatches( conserved, phastconsScore ) ) {
		toBePrinted.add( jInt );
		//System.out.println( "added: " + jInt.getPositionEnd() );
		return 1;
	    }
	}
	return 0;
    }

    private static List<Double> smoothWithMovingAverage( List<Double> data, int windowSize ) {
	if( windowSize < 1 ) throw new IllegalArgumentException( "Window size must be >= 1" );

	List<Double> smoothed = new ArrayList<>();
	int halfWindow = windowSize / 2;

	for( int i = 0; i < data.size(); i++ ) {
	    int start = Math.max( 0, i - halfWindow );
	    int end = Math.min( data.size() - 1, i + halfWindow );
	    
	    double sum = 0.0;
	    int count = 0;
	    for( int j = start; j <= end; j++ ) {
		sum += data.get( j );
		count++;
	    }
	    
	    smoothed.add( sum / count );
	}
	return smoothed;
    }

    private static List<Double> smoothWithGaussian( List<Double> data, double sigma ) {
	int kernelSize = (int) ( 6 * sigma ) + 1; // Define kernel size based on sigma
	if( kernelSize % 2 == 0 ) kernelSize++; // Ensure it's odd
	int halfKernel = kernelSize / 2;

	// Create Gaussian kernel
	double[] kernel = new double[kernelSize];
	double sum = 0.0;
	for( int i = 0; i < kernelSize; i++ ) {
	    double x = i - halfKernel;
	    kernel[i] = Math.exp( -0.5 * (x * x) / (sigma * sigma) );
	    sum += kernel[i];
	}
	for( int i = 0; i < kernelSize; i++ ) {
	    kernel[i] /= sum; // Normalize
	}

	// Apply convolution
	List<Double> smoothed = new ArrayList<>();
	for( int i = 0; i < data.size(); i++ ) {
	    double weightedSum = 0.0;
	    double weight = 0.0;

	    for( int j = -halfKernel; j <= halfKernel; j++ ) {
		int index = i + j;
		if( index >= 0 && index < data.size() ) {
		    weightedSum += data.get( index ) * kernel[j + halfKernel];
		    weight += kernel[j + halfKernel];
		}
	    }
	    smoothed.add( weightedSum / weight );
	}
	return smoothed;
    }

    private static List<Double> downsampleAndSmooth( List<Double> originalArray, int targetSize, double sigma ) {
	int originalSize = originalArray.size();
    
	// Down-sample the array to the target size
	List<Double> downsampledArray = new ArrayList<>();
	int binSize = originalSize / targetSize;
    
	for( int i = 0; i < targetSize; i++ ) {
	    int start = i * binSize;
	    int end = Math.min( start + binSize, originalSize );
	    double sum = 0.0;
	    for( int j = start; j < end; j++ ) {
		sum += originalArray.get( j );
	    }
	    downsampledArray.add( sum / ( end - start ) ); // Average values in each bin
	}
	// Apply Gaussian smoothing to the downsampled array
	return smoothWithGaussian( downsampledArray, sigma );
	//return smoothWithMovingAverage( downsampledArray, (int)sigma );
    }

    private static void computeConservationScore(
						 Guide g,
						 String query31mer,
						 String genomicSequence,
						 int gtfStart1Based, // original 'start' from GTF (1-based, inclusive)
						 int[] starts, int[] ends, float[] vals) {
	try {
	    // Find the 31-mer in the extracted genomic sequence, and map back to genomic coords
	    Map<String, Integer> adjustedCoordinates =
                search31merAndAdjustCoordinates(genomicSequence, gtfStart1Based, query31mer);

	    if (adjustedCoordinates != null) {
		// likely 1-based, closed
		int adjStart1 = adjustedCoordinates.get("start");
		int adjEnd1   = adjustedCoordinates.get("end");

		// convert to bedGraph coords: 0-based, half-open
		int qStart0 = adjStart1 - 1;
		int qEnd0   = adjEnd1;        // unchanged

		if (qStart0 < qEnd0) {
		    float score = fetchPhastConsScoresPrimitive( starts, ends, vals, qStart0, qEnd0 );
		    g.setPhastConsScore( Math.round( score * query31mer.length() ) );
		} else {
		    g.setPhastConsScore( -1 );
		}
	    } else {
		g.setPhastConsScore( -1 );
	    }
	} catch (Exception e) {
	    e.printStackTrace();
	    g.setPhastConsScore( -1 );
	}
    }

    // private static void computeConservationScoreOld( Guide g, String query31mer, String genomicSequence, int start, BedGraphLine[] bdLines ) {
    // 	// compute conservation score, need ENST without decimal id and 31mer target sequence
    // 	try {
    // 	    Map<String, Integer> adjustedCoordinates = search31merAndAdjustCoordinates( genomicSequence, start, query31mer );
    // 	    if( adjustedCoordinates != null ) {
    // 		int adjStart = adjustedCoordinates.get( "start" );
    // 		int adjEnd = adjustedCoordinates.get( "end" );
    // 		g.setPhastConsScore( fetchPhastConsScores( bdLines, adjStart, adjEnd ) );
    // 	    }
    // 	    else {
    // 		// Set<Character> setOfSymbolsInGS = new HashSet<>();
    // 		// for (char c : genomicSequence.toCharArray() )
    // 		//     setOfSymbolsInGS.add(c);
    // 		// Utils.debug( "symbols in genomic sequence: " + setOfSymbolsInGS );
    // 		// Utils.debug( "bad start: " + start + ", query31mer: " + query31mer + ", len( genomicSequence ): " + genomicSequence.length() );
    // 		g.setPhastConsScore( -1 ); // set to -1, no phastcons score extracted
    // 	    }
    // 	} catch( Exception e ) {
    // 	    e.printStackTrace();
    // 	}
    // }

    public static <K,V> void conditionalPut( Map<K, V> map, K key, V value, boolean condition, Runnable action ) {
	if( condition ) {
	    map.put( key, value ); // Add the key-value pair
	    action.run(); // Execute the independent action
	}
    }

    private static Guide makeGuide( int seedSite, int suppSite, int seedOffset, int suppOffset, ProteinCodingTranscript target, MatureMiRNA miRNA, int targetLen, boolean hasSeed, boolean hasSupp ) {
	int g2 = seedSite + 3 + seedOffset;
	//int g5 = seedSite;
	int g8 = seedSite - 3 + seedOffset;
	//int g12 = suppSite + 4;
	int g13 = suppSite + 3 + suppOffset;
	//int g16 = suppSite;
	int bridgeLength = g8 - g13 - 1;
	if( target.inTranscript( g2 - 29, g2 + 1 ) ) // if the target at g2 allows to define a valid 31mer...
	    if( bridgeLength >= 1 && bridgeLength <= 15 ) { // if the bridge is at most 15 nts...
		//Utils.debug( "buildin Guide... with " + miRNA.getMatureSequence() );
		Guide guide = new Guide( miRNA.getFullName(),
					 miRNA.getMIMATCode(),
					 target,
					 target.getGeneId(),
					 miRNA.getMatureSequence(),
					 g2,
					 g13,
					 StringSequence.reverseComplement( miRNA.getSeed() ),
					 StringSequence.reverseComplement( miRNA.getSupp() ),
					 bridgeLength,
					 g2 - targetLen + 2, // t31, aligned with 31mer + bridgeLength
					 g2 + 1, // t1
					 hasSupp ); // supplementary match provided?
		//Utils.debug( "folding guide..." );
		guide.fold();
		if( guide.bind() ) {
		    //Utils.debug( "afterfold, binding guide:\n" + guide );
		    return guide;
		}
		//else Utils.debug( "rejected by unbinding: " + guide.getSeedType() );
		//Utils.debug( "afterfold, unbinding guide:\n" + guide );
	    }
	//else Utils.debug( "rejected because of bridge length: " + bridgeLength );
	return null; // guide has too long a bridge or does not bind
    }
    
    public static void main( String[] args ) {

	// ***** JSON ObjectMapper
	ObjectMapper mapper = new ObjectMapper();

	// ***** READ ARGUMENTS

	if( args.length < 2  || args.length > 13 ) {
	    System.out.println( "usage: RISCScan <Gencode version> <ProteinCodingTranscript name> <MicroRNA name> [optional <Minimum seed> <Maximum Kd> <5'UTR> <CDS> <3'UTR> <Supplementary region paired?> <Seed accessibility> <Supplementary region accessibility> <Seed site conservation> <Force write>]" );
	    System.out.println( "   Gencode version: Gencode Protein Transcripts database version (only v46 is currently supported)" );
	    System.out.println( "   Transcript: Gencode Protein Transcript name" );
	    System.out.println( "   MicroRNA name or sequence: miRBase name or sequence" );
	    System.out.println( "   Minimum seed: minimum seed type to be shown in the output" );
	    System.out.println( "   max Kd: maximum Kd threshold (integer >= 0)" );
	    System.out.println( "   5'UTR: boolean to show (or not) the interactions in the 5'UTR" );
	    System.out.println( "   CDS: boolean to show (or not) the interactions in the CDS" );
	    System.out.println( "   3'UTR: boolean to show (or not) the interactions in the 3'UTR" );
	    System.out.println( "   Supplementary region paired: boolean to show the interactions with few or no supplementary pairing" );
	    System.out.println( "   Seed accessibility: real value (0 to 1) to show the interactions with accessibility equals or above the value" );
	    System.out.println( "   Supplementary region accessibility: real value (0 to 1) to show the interactions with accessibility equals or above the value" );
	    System.out.println( "   Seed site conservation: integer value (0 to 31) to show the interactions with conservation among 100 vertebrates equals or above the value" );
	    System.out.println( "   Force write: if 'flush', the results will be forced-written in the database (NOT IMPLEMENTED)" );
	    Utils.stop( "bye!" , 0 );
	}

	// ***** READ system variable for the data path
	String dataPath = System.getenv( "DATA_PATH" ); // path to the data directory (includes the final file separator)
	if( dataPath == null )
	    Utils.stop( "Unix DATA_PATH variable not set.", 17 );

	// get gencode version from args[0]
	String gencodeVersion = args[0];

	// get human genome file names
	final String HUMANGENOMEFASTAFILE = System.getenv( "HG_FASTA" );
	final String HUMANGENOMELOCATIONS = System.getenv( "HG_LOCATIONS" );
	final String HUMANGENOMEOUTPUTBUFFER = System.getenv( "HG_BUFFER" );

	// create a RIMapIO instance
	RIMapRISCIO riMapRISCIO = new RIMapRISCIO( dataPath, gencodeVersion );

	// set genomic utility files
	String fastaFile = HUMANGENOMEFASTAFILE; // human genome version (FASTA)

	// stuf related to computing conservation score
	String gtfFile = HUMANGENOMELOCATIONS;
	String outputFileRoot = HUMANGENOMEOUTPUTBUFFER;

	// set miRBase directory path
	String miRBasePath  = dataPath + RIMapRISCIO.MIRBASEDIRNAME;

	String pctName = args[1]; // Protein coding transcript name
	String miRName = args[2]; // MicroRNA name (MirBase) or guide sequence

	// none can be empty
	if( pctName.equals( "" ) || miRName.equals( "" ) )
	    Utils.stop( "Please, provide microRNA name or sequence, and transcript identifiers.", 4 );

	String minSeed = "seedless";
	if( args.length > 3 ) minSeed = args[3];

	double maxKd = 1500;
	if( args.length > 4 ) maxKd = Double.parseDouble( args[4] );

	boolean UTR5 = true;
	if( args.length > 5 ) UTR5 = Boolean.parseBoolean( args[5] );

	boolean CDS = true;
	if( args.length > 6 ) CDS = Boolean.parseBoolean( args[6] );

	boolean UTR3 = true;
	if( args.length > 7 ) UTR3 = Boolean.parseBoolean( args[7] );

	boolean supplementaryPaired = false;
	if( args.length > 8 ) supplementaryPaired = Boolean.parseBoolean( args[8] );

	double seedAccessible = 0.15;
	if( args.length > 9 ) seedAccessible = Double.parseDouble( args[9] );

	double supplementaryAccessible = 0.0;
	if( args.length > 10 ) supplementaryAccessible = Double.parseDouble( args[10] );

	int conserved = 0;
	if( args.length > 11 ) conserved = Integer.parseInt( args[11] );

	String forceWrite = "";
	if( args.length > 12 ) forceWrite = args[12];

	// ***** READ pct and miR

	// get the transcript and sequence
	ProteinCodingTranscript target = riMapRISCIO.getTranscript( pctName );
	String targetSequence = target.getStringSequence();

	// get the mirna and its sequence
	MatureMiRNA miRNA = riMapRISCIO.getGuide( miRName );
	miRName = miRNA.getName();
	String miRSequence = miRNA.getStringSequence();

	// interactions information
	List<Integer> targetRegions = new ArrayList<>();

	// transcript variables
	TiledSecondaryStructure target2D = null;

	// check if results exist
	if( riMapRISCIO.exist( pctName, miRSequence, targetSequence ) ) {
	    // process results
	    List<JsonInteraction> interactions = riMapRISCIO.getInteractions( pctName, miRSequence ); // read the JsonInteraction in a list
	    //System.out.println( "read " + interactions.size() + " interactions." );

	    // filter-out the interactions to be shown

	    int numberOfResults = 0; // to give feedback if no results are generated
	    List<JsonInteraction> toBePrinted = new ArrayList<>(); // interactions that pass selection criteria
	    if( target.getUTR5End() != -1 ) targetRegions.add( target.getUTR5End() );
	    targetRegions.add( target.getCDSEnd() );
	    if( target.getUTR3End() != -1 ) targetRegions.add( target.getUTR3End() );
	    // PRINT TARGET ID and GUIDE RNA
	    System.out.println( target.getName() + ":" + miRName );
	    // PRINT TARGET REGIONS
	    System.out.println( targetRegions );
	    for( JsonInteraction jInt: interactions )
		numberOfResults += printInteraction( toBePrinted, jInt, minSeed, maxKd, UTR5, CDS, UTR3, supplementaryPaired, seedAccessible, supplementaryAccessible, conserved );
	    if( numberOfResults == 0 ) {
		System.out.println( "No interaction to show!" );
	    }
	    else { // there are results
		// PRINT INTERACTIONS, POSITIONS and K_Ds
		toBePrinted.sort( Comparator.comparingDouble( JsonInteraction::getKd).thenComparingDouble( JsonInteraction::getDuplexEnergy ) );
		toBePrinted.forEach( i -> System.out.println( i ) );
		List<Integer> bindingPositions = new ArrayList<>();
		toBePrinted.forEach( i -> bindingPositions.add( i.getPositionStart() ) );
		System.out.println( "controler positions: " + bindingPositions );

		// get accessibility data; was previously computed
		String accessibilityFile = riMapRISCIO.getTranscriptFile( pctName ) + ".scan" + File.separator + "accessibilities"; // transcript accessibilities
		try {
		    target2D = new TiledSecondaryStructure( new FileReader( accessibilityFile ) );
		} catch( FileNotFoundException e ) {
		    Utils.stop( "Accessibility file not found: " + accessibilityFile, 8 );
		}
		// fill transcriptAccessibility with the values
		List<Double> transcriptAccessibility = Arrays.stream( target2D.getAccessibility() )
		    .boxed()
		    .collect( Collectors.toCollection(ArrayList::new));
		
		// smooth before sending
		DecimalFormat df = new DecimalFormat("#.####");
		List<Double> smoothedList = downsampleAndSmooth( transcriptAccessibility, 200, 2.0 );
		String formatted = smoothedList.stream()
		    .map( df::format )
		    .reduce( (a, b) -> a + "\t" + b )
		    .orElse("");
		System.out.println( "controler accessibilities: [" + formatted + "]" );
		//System.out.println( "controler accessibilities: " + downsampleAndSmooth( transcriptAccessibility, 200, 2.0 ) );
	    }
	    System.exit( 0 ); // no scan to be performed
	}
	else { // results must be generated, the file is already locked by riMapRISCIO.exist() and accessibility has been created if needed
	    // insert the JsonInteraction in a list
	    
	    // COMPUTE AND WRITE TO DB
	
	    // results file did not exist, but all variables are set to compute and write it
	    // ***** PHASE 1: find all seed and supp complementary sites *****
	    // seed and supp initiation sites of 4 nucleotides

	    List<JsonInteraction> interactions = new ArrayList<>();

	    // Prepare interaction information
	    targetRegions = new ArrayList<>();
	    if( target.getUTR5End() != -1 ) targetRegions.add( target.getUTR5End() );
	    targetRegions.add( target.getCDSEnd() );
	    if( target.getUTR3End() != -1 ) targetRegions.add( target.getUTR3End() );
	    // PRINT TARGET ID and GUIDE RNA
	    //System.out.println( target.getName() + ":" + miRNA.getFullName() + "(2)" );
	    // PRINT TARGET REGIONS
	    //System.out.println( targetRegions );

	    // Get genomic data
	    Map<String, String> location = null;
	    int start = 0; int end = 0;
	    String outputFile = outputFileRoot + "." + pctName + "." + miRName;
	    String strand = "";
	    String enstId = dropAfterDot( target.getGeneId() ); // ex) "ENST00000212015"
	    String seqname = "";
	    try {
		location = extractGenomicLocation( gtfFile, enstId );
		if( location != null ) {
		    seqname = "chr" + location.get( "seqname" );
		    start = Integer.parseInt(location.get( "start" ) );
		    end = Integer.parseInt(location.get( "end" ) );
		    strand = location.get( "strand" );

		    // Fetch sequence using bedtools
		    fetchSequence( pctName, miRName, dataPath, fastaFile, seqname, start, end, outputFile ); // save the sequence in outputFile
		}
	    } catch( Exception e ) {
		e.printStackTrace();
	    }

	    // read the lines of the bedfile into an array to allow for log n search
	    String bedGraphFile = dataPath + "bedGraph/" + seqname;

	    //List<String> bdEntries = new ArrayList<>();
	    //try {
	    //bdEntries = Files.readAllLines( Paths.get( bedGraphFile ) );
	    //} catch( IOException e ) {
	    //		e.printStackTrace();
	    //}
	    //BedGraphLine[] bdLines = new BedGraphLine[bdEntries.size()];
	    //arrayizeBedGraphLines( bdEntries, bdLines );

	    // avoiding all lines input

	    Path bgPath = Paths.get( bedGraphFile );

	    // --- PASS 1: count the valid bedGraph lines ---
	    int n = 0;
	    try( BufferedReader br = Files.newBufferedReader( bgPath ) ) {
		String line;
		while ((line = br.readLine()) != null) {
		    if (line.isEmpty() || line.charAt(0) == '#') continue;

		    // Expect: chrom <tab> start <tab> end <tab> value
		    int t1 = line.indexOf('\t');
		    if (t1 < 0) continue;
		    int t2 = line.indexOf('\t', t1 + 1);
		    if (t2 < 0) continue;
		    int t3 = line.indexOf('\t', t2 + 1);
		    if (t3 < 0) continue;

		    // Quick sanity parse (avoid counting garbage)
		    try {
			// We don’t keep values, just validate they parse
			Integer.parseInt(line, t1 + 1, t2, 10);
			Integer.parseInt(line, t2 + 1, t3, 10);
			Float.parseFloat(line.substring(t3 + 1));
			n++;
		    } catch (NumberFormatException ignore) {
			// skip malformed numeric line
		    }
		}
	    } catch (IOException e) {
		e.printStackTrace();
		n = 0;  // fallback; you can also return early
	    }

	    // allocate tight primitive arrays (memory-light)
	    int[]  starts = new int[n];
	    int[]  ends   = new int[n];
	    float[] vals  = new float[n];

	    // --- PASS 2: parse into primitives ---
	    int ii = 0;
	    try( BufferedReader br = Files.newBufferedReader( bgPath ) ) {
		String line;
		while( ( line = br.readLine() ) != null ) {
		    if( line.isEmpty() || line.charAt( 0 ) == '#' ) continue;

		    int t1 = line.indexOf( '\t' );
		    int t2 = ( t1 >= 0 ) ? line.indexOf( '\t', t1 + 1 ) : -1;
		    int t3 = ( t2 >= 0 ) ? line.indexOf( '\t', t2 + 1 ) : -1;
		    if( t1 < 0 || t2 < 0 || t3 < 0 ) continue;

		    try{
			int s = Integer.parseInt( line, t1 + 1, t2, 10 );
			int e = Integer.parseInt( line, t2 + 1, t3, 10 );
			float v = Float.parseFloat( line.substring( t3 + 1 ) );

			// guard against rare mismatch between pass1 and pass2
			if( ii < n ) {
			    starts[ii] = s;
			    ends[ii]   = e;
			    vals[ii]   = v;
			    ii++;
			} else {
			    // shouldn't happen, but protects from overflow
			    break;
			}
		    } catch( NumberFormatException ignore ) {
			// skip malformed numeric line
		    }
		}
	    } catch( IOException e ) {
		e.printStackTrace();
		// If you want, zero out arrays to signal failure
	    }

	    // If some lines were skipped in pass 2, shrink arrays (optional)
	    if( ii < n ) {
		starts = java.util.Arrays.copyOf( starts, ii );
		ends   = java.util.Arrays.copyOf( ends,   ii );
		vals   = java.util.Arrays.copyOf( vals,   ii );
	    }

	    // Convert the outputFile into RNA
	    try {
		convertSequenceToUpperCaseRNA( outputFile );
	    } catch( Exception e ) {
		e.printStackTrace();
	    }

	    // get the genomic sequence
	    String genomicSequence = "";
	    try( BufferedReader br = Files.newBufferedReader( Paths.get( outputFile ) ) ) {
		// skip FASTA header
		br.readLine();
		StringBuilder sb = new StringBuilder( Math.max( 1024, end - start + 32 ) );
		String line;
		while( ( line = br.readLine() ) != null ) {
		    sb.append( line.trim() );
		}
		genomicSequence = sb.toString();
	    } catch (IOException e) {
		e.printStackTrace();
	    }

	    int bedStart0 = start - 1;  // 0-based
	    int bedEnd0   = end;        // half-open

	    int L = genomicSequence.length();
	    float[] accessibility = new float[L];

	    // get the genomic sequence
	    // String genomicSequence = "";
	    // try {
	    // 	List<String> lines = Files.readAllLines( Paths.get( outputFile ) );
	    // 	genomicSequence = lines.stream().skip( 1 ).collect( Collectors.joining() ).replaceAll( "\\n", "" );
	    // } catch( Exception e ) {
	    // 	e.printStackTrace();
	    // }

	    // delete the output file
	    try {
		Files.deleteIfExists( Paths.get( outputFile ) );
	    } catch( IOException e ) {
		System.err.println( "Could not delete file: " + outputFile );
		e.printStackTrace();
	    }

	    // get accessibility data; was previously computed
	    String accessibilityFile = riMapRISCIO.getTranscriptFile( pctName ) + ".scan" + File.separator + "accessibilities"; // transcript accessibilities
	    try {
		target2D = new TiledSecondaryStructure( new FileReader( accessibilityFile ) );
	    } catch( FileNotFoundException e ) {
		Utils.stop( "Accessibility file not found: " + accessibilityFile, 11 );
	    }
	    // fill transcriptAccessibility with the values
	    List<Double> transcriptAccessibility = Arrays.stream( target2D.getAccessibility() )
		.boxed()
		.collect( Collectors.toCollection(ArrayList::new));

	    // compute grips: 4mer in the seed (4 WC pairs); 4mer in the supp (4 WC+Wobble pairs); 7mer with 1 default at g5, g6 or g7

	    // seeds
	    
	    List<List<Integer>> seedSites = new ArrayList<>();
	    //System.out.println( "computing seed complement for: " + miRNA.getSeed() );
	    List<Integer> seedInitiationSites = StringSequence.complementSites( targetSequence, miRNA.getSeed(), false ); // get g2s

	    // generate offset seed sites, keep new sites only
	    Set<Integer> siteBag = new HashSet<>( seedInitiationSites ); // transfer to a set for fast inclusion checks
	    //System.out.println( "computing offset seed complement for: " + miRNA.getOffsetSeed() );
	    List<Integer> seedOffsetSites = StringSequence.complementSites( targetSequence, miRNA.getOffsetSeed(), false ); // get g3s
	    Set<Integer> garbageBag = new HashSet<>(); // put the same sites in a garbage bag
	    for( Integer offsetSite : seedOffsetSites )
		if( siteBag.contains( offsetSite + 1 ) )
		    garbageBag.add( offsetSite );;
	    for( Integer garbage : garbageBag ) // remove the same sites from seedOffsetSites
		seedOffsetSites.remove( garbage );

	    siteBag.addAll( seedOffsetSites ); // add the offset sites to the seed set
	    // generate distant seed sites, keep new sites only
	    //System.out.println( "computing distant seed complement for: " + miRNA.getDistantSeed() );
	    List<Integer> seedDistantSites = StringSequence.complementSites( targetSequence, miRNA.getDistantSeed(), false ); // get g4s
	    garbageBag.clear();
	    for( Integer distantSite : seedDistantSites )
		if( siteBag.contains( distantSite + 2 ) || siteBag.contains( distantSite + 1 ) )
		    garbageBag.add( distantSite );
	    for( Integer garbage : garbageBag )
		seedDistantSites.remove( garbage );

	    // put the seed sites in the seedSites list of lists of seed sites
	    seedSites.add( seedInitiationSites );
	    seedSites.add( seedOffsetSites );
	    seedSites.add( seedDistantSites );

	    // supps
	    
	    List<List<Integer>> suppSites = new ArrayList<>();
	    //System.out.println( "computing supp sites for " + miRNA.getSupp() );
	    List<Integer> suppInitiationSites = StringSequence.complementSites( targetSequence, miRNA.getSupp(), true ); // get g13s

	    // generate offset supp sites, keep new sites only
	    //System.out.println( "computing supp sites for " + miRNA.getOffsetSupp() );
	    List<Integer> suppOffsetSites = StringSequence.complementSites( targetSequence, miRNA.getOffsetSupp(), true ); // get g14s

	    suppSites.add( suppInitiationSites );
	    //System.out.println( "supp sites: " + suppInitiationSites ); // DEBUG
	    suppSites.add( suppOffsetSites );
	    //System.out.println( "supp offset sites: " + suppOffsetSites ); // DEBUG

	    int targetLen = 31; // for all types

	    Map<Integer,Guide> optimalSites = new HashMap<>(); // map of optimal seed sites
	    Map<Integer,Double> optimalEnergy = new HashMap<>(); // map of optimal energy for optimal sites

	    // special 7mers with one default

	    // 7mers with 1 wobble or mismatch at G4 or G5 (matching 6mers in target)

	    //Utils.debug( "miRNA full seed: " + miRNA.getFullSeed() );
	    List<Integer> specialSites = StringSequence.complementSites( targetSequence, miRNA.getFullSeed() ); // get g2s
	    //Utils.debug( "special sites: " + specialSites ); // DEBUG

	    // iterate on all special and supp sites; adjust supp sites by offset
	    for( Integer specialSite : specialSites ) {
		for( int suppOffset = 0; suppOffset < 2; suppOffset++ ) for( Integer suppSite : suppSites.get( suppOffset ) ) {
			Guide g = makeGuide( specialSite, suppSite, 3, suppOffset, target, miRNA, targetLen, true, true ); // 2-grip sites; has seed and supp
			if( g != null ) {
			    if( g.getDeltaG() < optimalEnergy.getOrDefault( specialSite, 0.0 ) ) {
				optimalEnergy.put( specialSite, g.getDeltaG() );
				optimalSites.put( specialSite, g );
				//Utils.debug( "accept special site " + specialSite );
				//Utils.debug( g.toString() );
			    }
			}
			//else Utils.debug( "reject special site " + specialSite );
		    }
	    }

	    // iterate on all special sites (no supp); adjust supp from seed (seed-only sites)
	    for( Integer specialSite : specialSites ) {
		int suppSite = specialSite - 11;
		Guide g = makeGuide( specialSite, suppSite, 3, 0, target, miRNA, targetLen, true, false ); // single grip seed side; has seed no supp
		if( g != null ) {
		    if( g.getDeltaG() < optimalEnergy.getOrDefault( specialSite, 0.0 ) ) {
			optimalEnergy.put( specialSite, g.getDeltaG() );
			optimalSites.put( specialSite, g );
			//Utils.debug( "accept special site " + specialSite );
			//Utils.debug( g.toString() );
		    }
		    //else Utils.debug( "reject special site " + specialSite );
		}
	    }

	    // 7mers with 1 bulge at G4 or G5 (deletion in target)

	    //Utils.debug( "miRNA full seed: " + miRNA.getFullSeed() );
	    List<Integer> deleteSites = StringSequence.complementWithDeletionAt2or3( targetSequence, miRNA.getFullSeed() ); // get g2s
	    //Utils.debug( "delete sites: " + deleteSites ); // DEBUG

	    // iterate on all delete and supp sites; adjust supp sites by offset
	    for( Integer deleteSite : deleteSites ) {
		for( int suppOffset = 0; suppOffset < 2; suppOffset++ ) for( Integer suppSite : suppSites.get( suppOffset ) ) {
			Guide g = makeGuide( deleteSite, suppSite, 2, suppOffset, target, miRNA, targetLen, true, true ); // 2-grip sites; has seed and supp
			if( g != null ) {
			    if( g.getDeltaG() < optimalEnergy.getOrDefault( deleteSite, 0.0 ) ) {
				optimalEnergy.put( deleteSite, g.getDeltaG() );
				optimalSites.put( deleteSite, g );
				//Utils.debug( "accept delete site with supp" + deleteSite );
				//Utils.debug( g.toString() );
			    }
			}
			//else Utils.debug( "reject delete site with supp " + deleteSite );
		    }
	    }

	    // iterate on all delete sites (no supp); adjust supp from seed (seed-only sites)
	    for( Integer deleteSite : deleteSites ) {
		int suppSite = deleteSite - 11;
		Guide g = makeGuide( deleteSite, suppSite, 2, 0, target, miRNA, targetLen, true, false ); // single grip seed side; has seed no supp
		if( g != null ) {
		    if( g.getDeltaG() < optimalEnergy.getOrDefault( deleteSite, 0.0 ) ) {
			optimalEnergy.put( deleteSite, g.getDeltaG() );
			optimalSites.put( deleteSite, g );
			//Utils.debug( "accept delete site without supp " + deleteSite );
			//Utils.debug( g.toString() );
		    }
		    //else Utils.debug( "reject delete site without supp " + deleteSite );
		}
	    }
	    
	    // 7mers with 1 bulge in target b(4.5) (insertion in target)

	    //Utils.debug( "miRNA full seed: " + miRNA.getFullSeed() );
	    List<Integer> insertSites = StringSequence.complementWithInsertionAt3( targetSequence, miRNA.getFullSeed() ); // get g2s
	    //Utils.debug( "insert sites: " + insertSites ); // DEBUG

	    // iterate on all insert and supp sites; adjust supp sites by offset
	    for( Integer insertSite : insertSites ) {
		for( int suppOffset = 0; suppOffset < 2; suppOffset++ ) for( Integer suppSite : suppSites.get( suppOffset ) ) {
			Guide g = makeGuide( insertSite, suppSite, 4, suppOffset, target, miRNA, targetLen, true, true ); // 2-grip sites; has seed and supp
			if( g != null ) {
			    if( g.getDeltaG() < optimalEnergy.getOrDefault( insertSite, 0.0 ) ) {
				optimalEnergy.put( insertSite, g.getDeltaG() );
				optimalSites.put( insertSite, g );
				//Utils.debug( "accept delete site with supp" + deleteSite );
				//Utils.debug( g.toString() );
			    }
			}
			//else Utils.debug( "reject insert site with supp " + insertSite );
		    }
	    }

	    // iterate on all insert sites (no supp); adjust supp from seed (seed-only sites)
	    for( Integer insertSite : insertSites ) {
		int suppSite = insertSite - 11;
		Guide g = makeGuide( insertSite, suppSite, 4, 0, target, miRNA, targetLen, true, false ); // single grip seed side; has seed no supp
		if( g != null ) {
		    if( g.getDeltaG() < optimalEnergy.getOrDefault( insertSite, 0.0 ) ) {
			optimalEnergy.put( insertSite, g.getDeltaG() );
			optimalSites.put( insertSite, g );
			//Utils.debug( "accept insert site without supp " + deleteSite );
			//Utils.debug( g.toString() );
		    }
		    //else Utils.debug( "reject insert site without supp " + deleteSite );
		}
	    }

	    // grips

	    // iterate on all seed and supp sites; adjust seed and supp sites by offset
	    for( int seedOffset = 0; seedOffset < 3; seedOffset++ ) for( Integer seedSite : seedSites.get( seedOffset ) ) {
		    for( int suppOffset = 0; suppOffset < 2; suppOffset++ ) for( Integer suppSite : suppSites.get( suppOffset ) ) {
			    //Guide seedg = makeGuide( seedSite, seedSite - 11, seedOffset, target, miRNA, targetLen ); // single grip seed side; suppSite at g13 by default
			    Guide g = makeGuide( seedSite, suppSite, seedOffset, suppOffset, target, miRNA, targetLen, true, true ); // 2-grip sites; has seed and supp
			    //Guide suppg = makeGuide( suppSite + 10 + suppOffset, suppSite, seedOffset, target, miRNA, targetLen ); // single grip supp side (seedless); suppSites at g15, g16, and g17
			    if( g != null ) {
				//System.out.println( "seed+supp, seedSite added: " + seedSite );
				if( g.getDeltaG() < optimalEnergy.getOrDefault( seedSite, 0.0 ) ) {
				    optimalEnergy.put( seedSite, g.getDeltaG() );
				    optimalSites.put( seedSite, g );
				    //System.out.println( "saving deltaG = " + g.getDeltaG() + ", suppSite = " + suppSite );
				}
			    }
			}
		}

	    // iterate on all seed; adjust supp from seed (seed-only sites)
	    for( int seedOffset = 0; seedOffset < 3; seedOffset++ ) for( Integer seedSite : seedSites.get( seedOffset ) ) {
		    int suppSite = seedSite - 11 + seedOffset;
		    Guide g = makeGuide( seedSite, suppSite, seedOffset, 0, target, miRNA, targetLen, true, false ); // single grip seed side; has seed no supp
		    if( g != null ) {
			//System.out.println( "seed-only, seedSite added: " + seedSite );
			if( g.getDeltaG() < optimalEnergy.getOrDefault( seedSite, 0.0 ) ) {
			    optimalEnergy.put( seedSite, g.getDeltaG() );
			    optimalSites.put( seedSite, g );
			    //System.out.println( "saving deltaG = " + g.getDeltaG() + ", suppSite = " + suppSite );
			}
		    }
		}

	    // iterate on all supp; adjust seed from supp(seedless sites)
	    for( int suppOffset = 0; suppOffset < 2; suppOffset++ ) for( Integer suppSite : suppSites.get( suppOffset ) ) {
		    int seedSite = suppSite + 11 + suppOffset;
		    // boolean alreadyDone = false;
		    // for( int i = seedSite - 2; i <= seedSite + 2; i++ )
		    // 	if( optimalSites.containsKey( i ) ) {
		    // 	    alreadyDone = true;
		    // 	    break;
		    // 	}
		    // if( !alreadyDone ) {
		    //Utils.debug( "makeGuide( " + seedSite + ", " + suppSite + ", " + 0 + ", " + suppOffset + ", " + target + ", " + targetLen + " )" );
			Guide g = makeGuide( seedSite, suppSite, 0, suppOffset, target, miRNA, targetLen, false, true ); // single grip seed side; has no seed, has supp
			if( g != null ) {
			    //Utils.debug( "guide:\n" + g );
			    // check if new seed contains 4mer
			    if( !g.has4merInSeed() )
				if( g.getDeltaG() < optimalEnergy.getOrDefault( seedSite, 0.0 ) ) {
				    optimalEnergy.put( seedSite, g.getDeltaG() );
				    optimalSites.put( seedSite, g );
				    //System.out.println( "saving deltaG = " + g.getDeltaG() + ", suppSite = " + suppSite );
				}
			}
			//}
		}

	    // fill the interactions
	    for( Guide g : optimalSites.values() ) {
		//System.out.println( "g.g2: " + g.getg2() );
		computeConservationScore( g, targetSequence.substring( g.getg2() - targetLen + 2, g.getg2()+2 ), genomicSequence, start, starts, ends, vals );
		g.setSeedAccessibility( target2D.getAccessibility( g.getg5(), g.getg2() ) ); // g2-g5
		g.setSuppAccessibility( target2D.getAccessibility( g.getg16(), g.getg13() ) ); // g13-g16
		try {
		    String jsonString = mapper.writeValueAsString( g );
		    //scanResultsWriter.write( jsonString + "\n" ); // write human readable to stdout
		    ObjectMapper objectMapper = new ObjectMapper();
		    JsonInteraction jInt = objectMapper.readValue( jsonString, JsonInteraction.class ); // Convert JSON string to JsonInteraction
		    interactions.add( jInt );
		    //System.out.println( "successful g.g2: " + g.getg2() );
		} catch( Exception e ) {
		    e.printStackTrace();
		}
	    }
	
	    // IO to the database
	    riMapRISCIO.putInteractions( pctName, miRSequence, interactions ); // create the output file, putInteractions remove the lock file
	    //System.out.println( "put " + interactions.size() + " interactions." );
	    
	    // filter-out the interactions to be shown

	    int numberOfResults = 0; // to give feedback if no results are generated
	    List<JsonInteraction> toBePrinted = new ArrayList<>(); // interactions that pass selection criteria
	    if( target.getUTR5End() != -1 ) targetRegions.add( target.getUTR5End() );
	    targetRegions.add( target.getCDSEnd() );
	    if( target.getUTR3End() != -1 ) targetRegions.add( target.getUTR3End() );
	    // PRINT TARGET ID and GUIDE RNA
	    System.out.println( target.getName() + ":" + miRName );
	    // PRINT TARGET REGIONS
	    System.out.println( targetRegions );
	    for( JsonInteraction jInt: interactions )
		numberOfResults += printInteraction( toBePrinted, jInt, minSeed, maxKd, UTR5, CDS, UTR3, supplementaryPaired, seedAccessible, supplementaryAccessible, conserved );
	    if( numberOfResults == 0 ) {
		System.out.println( "No interaction to show!" );
	    }
	    else { // there are results
		// PRINT INTERACTIONS, POSITIONS and K_Ds
		toBePrinted.sort( Comparator.comparingDouble( JsonInteraction::getKd).thenComparingDouble( JsonInteraction::getDuplexEnergy ) );
		toBePrinted.forEach( i -> System.out.println( i ) );
		List<Integer> bindingPositions = new ArrayList<>();
		toBePrinted.forEach( i -> bindingPositions.add( i.getPositionStart() ) );
		System.out.println( "controler positions: " + bindingPositions );
		
		DecimalFormat df = new DecimalFormat("#.####");
		List<Double> smoothedList = downsampleAndSmooth( transcriptAccessibility, 200, 2.0 );
		String formatted = smoothedList.stream()
		    .map( df::format )
		    .reduce( (a, b) -> a + "\t" + b )
		    .orElse("");		
		System.out.println("controler accessibilities: [" + formatted + "]" );
		//System.out.println( "controler accessibilities: " + downsampleAndSmooth( transcriptAccessibility, 200, 2.0 ) );
	    }
	}

	// come back here to implement the flush option
	// 	    if( forceWrite.equals( "flush" ) ) {
	// 	System.out.println( "Option force write used: results will be written to the database" );
	// 	Files.deleteIfExists( filePath );
	// 	System.out.println( "Opening scan results file: " + scanResultsFile );
	// 	scanResultsWriter = new FileWriter( scanResultsFile ); // create result file
	// }
    }
}
