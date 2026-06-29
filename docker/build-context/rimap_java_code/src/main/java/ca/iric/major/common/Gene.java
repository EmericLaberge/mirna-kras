/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Interface for genes
 *   
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2023 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.lang.Iterable;
import java.util.Iterator;

public interface Gene extends Iterable<Transcript> {

    /**
     *
     * @return the iterator of corresponding transcripts
     */    
    public Iterator<Transcript> iterator();

    @Override
    String toString();
}
