/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Interface for mutable RNA sequences
 *     instances can be assigned multiple times using setSequence
 *     and can be modified
 *   
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2023 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

public interface MutableSequence extends Sequence {

    @Override
    String toString();
}
