/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Class for managing subsets of a set
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 2024 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

// Use LinkdedHashSet to preserve the order of the elements as they were entered.

import java.util.Set;
import java.util.LinkedHashSet;

public final class AdvancedSet {

    private AdvancedSet() { throw new UnsupportedOperationException( "Utility class AdvanceSet cannot be instantiated" ); }
    
    // Power set (all subsets) of the set
    public static <E> Set<Set<E>> getPowerSet( Set<E> list ) {
        Set<Set<E>> powerSet = new LinkedHashSet<>();
        powerSetHelper( list, new LinkedHashSet<E>(), powerSet );
        return powerSet;
    }

    // Helper method to generate the power set recursively
    private static <E> void powerSetHelper( Set<E> originalSet, Set<E> currentSet, Set<Set<E>> powerSet ) {
        if( originalSet.isEmpty() ) {
            powerSet.add( new LinkedHashSet<E>( currentSet ) );
            return;
        }
        E element = originalSet.iterator().next();
        Set<E> remainingSet = new LinkedHashSet<>( originalSet );
        remainingSet.remove( element );

        // Exclude the current element
        powerSetHelper( remainingSet, currentSet, powerSet );

        // Include the current element
        currentSet.add( element );
        powerSetHelper( remainingSet, currentSet, powerSet );
        currentSet.remove( element );
    }

    public static <E> Set<E> union( Set<E> A, Set<E> B ) {
	Set<E> unionSet = new LinkedHashSet<>( A );
	unionSet.addAll( B );
	return unionSet;
    }

    public static <E> Set<E> union( Set<E> A, E element ) {
	Set<E> unionSet = new LinkedHashSet<>( A );
	unionSet.add( element );
	return unionSet;
    }
}
