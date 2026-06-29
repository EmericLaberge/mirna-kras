/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import ca.iric.major.common.Sequence;

/**
 * Interface for transcripts
 *   
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2023 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

public interface Transcript {

    /**
     * Returns the transcript type (based on vega.archive.ensembl)
     *
     * @return the transcript type
     */
    public Biotype getBiotype();

    /**
     * Returns the sequence of the transcript
     *
     * @return the sequence
     */
    public Sequence getSequence();

    /**
     * Returns the sequence of the transcript as a String
     *
     * @return the sequence as a String
     */
    public String getStringSequence();

    /**
     *
     * Set the Biotype of the transcript
     *
     * @param the Biotype
     */
    public void setBiotype( Biotype type );

    /**
     *
     * Set the sequence of the transcript
     *
     * @param the RNA Sequence
     */
    public void setSequence( Sequence sequence );

    @Override
    String toString();
}
