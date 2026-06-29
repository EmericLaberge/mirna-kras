/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Interface for RNA sequences
 *     instances are immutable, but
 *     for being changed at once using setSequence
 *   
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2023 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

public interface Sequence {

    /**
     * Returns the RNA sequence as a String
     *
     * @return the RNA sequence
     */
    public String getSequence();

    /**
     * Returns the subsequence as a String
     *
     * @return the RNA subsequence
     * @param the interval of the subsequence [start,end[
     * @throws IllegalArgumentException if the interval is outside the sequence
     */
    public String getSequence( int start, int end );

    /**
     * Returns the nucleotide as a character
     *
     * @return the nucleotide
     * @param the position of the nucleotide, position
     * @throws IllegalArgumentException if the position is outside the sequence
     */
    public char getNucleotide( int position );

    /**
     * Returns the length of the RNA sequence
     *
     * @return the legnth of the RNA sequence
     */
    public int length();

    /**
     * Sets the RNA sequence to a new value
     *     all characters converted to upper case, A, C, G, and U only
     *
     * @param sequence the new RNA sequence
     * @throws IllegalArgumentException if the sequence contains invalid characters
    */
    public void setSequence( String sequence ) throws IllegalArgumentException;
    
    @Override
    String toString();
}
