/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.List;
import java.util.LinkedList;
import java.util.stream.IntStream;
import java.util.concurrent.locks.ReentrantLock;
import java.util.Arrays;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.text.DecimalFormat;

public class TiledSecondaryStructure {

    private int[] tileSizes;
    private double[] tileE;
    private double threshold;

    protected String strand;
    //    protected SecondaryStructure[] tiles;
    protected int[] canonicalBasePairCount;
    protected int[] stateCount;
    protected double[] canonicalBasePairProbability;

    public char getNucleotide( int i ) { return this.strand.charAt( i ); }

    public double[] getBasePairProbabilities() {
	return this.canonicalBasePairProbability.clone();
    }

    public double getReactivity( int nt ) {
	return 1 - this.canonicalBasePairProbability[nt];
    }

    public double getAccessibility( int begin, int end ) {
	double res = 1.0;
	for( int i = begin; i <= end; i++ )
	    if( this.canonicalBasePairProbability[i] < res )
		res = this.canonicalBasePairProbability[i];
	return 1 - res;
    }

    public double[] getAccessibility() {
	double[] res = new double[this.canonicalBasePairProbability.length];
	for( int i = 0; i < this.canonicalBasePairProbability.length; i++ ) res[i] = 1 - this.canonicalBasePairProbability[i];
	return res;
    }

    // Tiled folding on a substring
    public TiledSecondaryStructure( String strand, int left, int right ) {
	this.strand = strand.substring( left, right + 1 );
	//System.out.println( this.strand );
	//System.out.println( "TiledSecondaryStructure for a strand of length: " + this.strand.length() + " right: " + right + " left: " + left );
	int tiny = this.strand.length() / 5;
	int smal = this.strand.length() / 4;
	int bigg = this.strand.length() / 3;
	this.tileSizes = new int[] { tiny, smal, bigg };
	this.tileE = new double[] { 3, 3, 3 };
        this.canonicalBasePairCount = new int[this.strand.length()];
	this.stateCount = new int[this.strand.length()];
	this.canonicalBasePairProbability = new double[this.strand.length()];


	IntStream.range( 0, tileSizes.length ).parallel().forEach( i -> {
		if( tileSizes[i] <= this.strand.length() ) {
		    IntStream.range( 0, this.strand.length() - tileSizes[i] + 1 ).forEach( start -> {
			    SecondaryStructure tile = new SecondaryStructure( this.strand.substring( start, start + tileSizes[i]), tileE[i] );
			    for( int k = 0; k < tileSizes[i]; k++ ) {
				synchronized( this ) {
				    this.canonicalBasePairCount[k + start] += tile.getCanonicalCount( k );
				    this.stateCount[k + start] += tile.getNumberOfStates();
				}
			    }
			});
		}
	    });
	// adjust the probabilities
	for( int i = 0; i < this.strand.length(); i++ )
	    this.canonicalBasePairProbability[i] = (double)this.canonicalBasePairCount[i] / this.stateCount[i];
	}

    // Tiled folding on a string
    public TiledSecondaryStructure( String strand, int[] tileSizes, double[] tileT, double threshold ) {
	//System.out.println( "TiledSecondaryStructure for a strand of length: " + strand.length() );
	this.threshold = threshold;
	this.strand = strand;
        this.canonicalBasePairCount = new int[this.strand.length()];
	this.stateCount = new int[this.strand.length()];
	this.canonicalBasePairProbability = new double[this.strand.length()];


	IntStream.range( 0, tileSizes.length ).parallel().forEach( i -> {
		if( tileSizes[i] <= this.strand.length() ) {
		    IntStream.range( 0, this.strand.length() - tileSizes[i] + 1 ).forEach( start -> {
			    SecondaryStructure tile = new SecondaryStructure( this.strand.substring( start, start + tileSizes[i]), tileT[i] );
			    for( int k = 0; k < tileSizes[i]; k++ ) {
				synchronized( this ) {
				    this.canonicalBasePairCount[k + start] += tile.getCanonicalCount( k );
				    this.stateCount[k + start] += tile.getNumberOfStates();
				}
			    }
			});
		}
	    });
	// adjust the probabilities
	for( int i = 0; i < this.strand.length(); i++ )
	    this.canonicalBasePairProbability[i] = (double)this.canonicalBasePairCount[i] / this.stateCount[i];
    }

    public TiledSecondaryStructure(String strand) {
	this.strand = strand;

	// example defaults – you can replace with whatever set you like
	// this.tileSizes = new int[] { 50, 65, 125 };
	// this.tileE     = new double[] { 2.0, 1.0, 1.0 };
	this.tileSizes = new int[] { 71 };       // 20 + 31 + 20
	this.tileE     = new double[] { 1.0 };

	if (tileSizes.length != tileE.length)
	    throw new IllegalArgumentException("tileSizes and tileE must be the same length");

	final int n = strand.length();

	this.canonicalBasePairCount       = new int[n];
	this.stateCount                   = new int[n];
	this.canonicalBasePairProbability = new double[n];

	if (n == 0) {
	    Arrays.fill(this.canonicalBasePairProbability, Double.NaN);
	    return;
	}

	// Striped locks for parallel merging
	final int NUM_LOCKS = 256;
	final ReentrantLock[] locks = new ReentrantLock[NUM_LOCKS];
	for (int i = 0; i < NUM_LOCKS; i++) locks[i] = new ReentrantLock();

	// Loop over each tile size / energy pair
	for (int t = 0; t < tileSizes.length; t++) {
	    final int tileSize     = tileSizes[t];
	    final double tileEnergy = tileE[t];

	    if (tileSize <= 0) continue;
	    if (n < tileSize)  continue;

	    IntStream.range(0, n - tileSize + 1).parallel().forEach(start -> {
		    final String subseq = strand.substring(start, start + tileSize);
		    final SecondaryStructure tile = new SecondaryStructure(subseq, tileEnergy);
		    

		    final int[] localCBC = new int[tileSize];
		    for (int k = 0; k < tileSize; k++) {
			localCBC[k] = tile.getCanonicalCount(k);
		    }
		    final int states = tile.getNumberOfStates();

		    for (int k = 0; k < tileSize; k++) {
			final int idx = start + k;
			final ReentrantLock L = locks[idx & (NUM_LOCKS - 1)];
			L.lock();
			try {
			    canonicalBasePairCount[idx] += localCBC[k];
			    stateCount[idx]            += states;
			} finally {
			    L.unlock();
			}
		    }
		});
	}

	// finalize probabilities
	for (int i = 0; i < n; i++) {
	    canonicalBasePairProbability[i] =
		stateCount[i] > 0 ? (double) canonicalBasePairCount[i] / stateCount[i]
		: Double.NaN;
	}
    }
    
    // Tiled folding on a string
    public void TiledSecondaryStructure( String strand ) {
	//this.tileSizes = new int[] { 50, 65, 125 };
	//this.tileE = new double[] { 2, 1, 1 };
	this.tileSizes = new int[] { 71 }; // 20 + 31 + 20 (20 nts around 31mers)
	this.tileE = new double[] { 1 };
	//System.out.println( "TiledSecondaryStructure for a strand of length: " + strand.length() );
	this.strand = strand;
        this.canonicalBasePairCount = new int[this.strand.length()];
	this.stateCount = new int[this.strand.length()];
	this.canonicalBasePairProbability = new double[this.strand.length()];


	IntStream.range( 0, tileSizes.length ).parallel().forEach( i -> {
		if( tileSizes[i] <= this.strand.length() ) {
		    IntStream.range( 0, this.strand.length() - tileSizes[i] + 1 ).forEach( start -> {
			    SecondaryStructure tile = new SecondaryStructure( this.strand.substring( start, start + tileSizes[i]), tileE[i] );
			    for( int k = 0; k < tileSizes[i]; k++ ) {
				synchronized( this ) {
				    this.canonicalBasePairCount[k + start] += tile.getCanonicalCount( k );
				    this.stateCount[k + start] += tile.getNumberOfStates();
				}
			    }
			});
		}
	    });
	// adjust the probabilities
	for( int i = 0; i < this.strand.length(); i++ )
	    this.canonicalBasePairProbability[i] = (double)this.canonicalBasePairCount[i] / this.stateCount[i];
    }

    // Tiled results from file
    public TiledSecondaryStructure( FileReader fileR ) {
        try( BufferedReader br = new BufferedReader( fileR ) ) {
            // Read the first line (assume all values are on a single line)
            String line = br.readLine();
            if( line != null ) {
                // Split the line by tab delimiter
                String[] values = line.split("\t");

                // Initialize the array with the size of the values
                this.canonicalBasePairProbability = new double[values.length];

                // Parse each value into a double and assign it to the array
                for( int i = 0; i < values.length; i++ ) {
                    this.canonicalBasePairProbability[i] = Double.parseDouble(values[i]);
                }
            }
        } catch (IOException e) {
            System.err.println( "Error reading the file: " + e.getMessage() );
        } catch( NumberFormatException e ) {
            System.err.println( "Error parsing a value as a double: " + e.getMessage() );
        }
    }

    @Override
    public String toString() {
	DecimalFormat df = new DecimalFormat( "#.####" );
	StringBuilder out = new StringBuilder();
	for( int i = 0; i < canonicalBasePairProbability.length; i++ ) {
	    out.append( df.format( canonicalBasePairProbability[i] ) );
	    if( i < canonicalBasePairProbability.length - 1 )
		out.append( "\t" );
	}
	return out.toString();
    }

    public String toCSV() {
	String out = "";
	for( int i = 0; i < this.canonicalBasePairProbability.length - 1; i++ )
	    out += this.canonicalBasePairProbability[i] + " ";
	out += this.canonicalBasePairProbability[this.canonicalBasePairProbability.length-1];
	return out;
    }
    
    public String toPretty() {
	String out = "";
	for( int i = 0; i < this.canonicalBasePairProbability.length; i++ )
	    out += this.canonicalBasePairProbability[i] > this.threshold ? "|" : "."; 
	return out;
    }
}
