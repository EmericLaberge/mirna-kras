/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Objects;
import java.util.Set;
import java.util.stream.Stream;
import java.util.stream.Collectors;

/**
 * TranscriptSetKey produces a key from three Set<ProteinCodingTranscript>
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 1.0 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

public class TranscriptSetKey {
    private final Set<ProteinCodingTranscript> targets;
    private final Set<ProteinCodingTranscript> options;
    private final Set<ProteinCodingTranscript> excludes;

    public TranscriptSetKey( Set<ProteinCodingTranscript> targets,
                             Set<ProteinCodingTranscript> options,
                             Set<ProteinCodingTranscript> excludes) {
        this.targets  = Set.copyOf( targets );   // defensively copy
        this.options  = Set.copyOf( options );
        this.excludes = Set.copyOf( excludes );
    }

    @Override
    public boolean equals( Object o ) {
        if( this == o ) return true;
        if ( !(o instanceof TranscriptSetKey) ) return false;
        TranscriptSetKey that = (TranscriptSetKey) o;
        return Objects.equals( this.targets, that.targets ) &&
               Objects.equals( this.options, that.options) &&
               Objects.equals( this.excludes, that.excludes);
    }

    @Override
    public int hashCode() {
        return Objects.hash( this.targets, this.options, this.excludes );
    }

    @Override
    public String toString() {
	return "{" +
	    Stream.concat(
			  Stream.concat( this.targets.stream(), this.options.stream() ), this.excludes.stream() )
	    .map( ProteinCodingTranscript::toString ) // Or use .map(t -> t.getName()) if you prefer specific fields
	    .collect( Collectors.joining( ", " ) )
	    + "}";
    }
}
