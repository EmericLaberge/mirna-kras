/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Class for Reads (composed of an id and a Sequence)
 *
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2024 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.lang.IllegalArgumentException;

public class Read {

    // composed of Id and Sequence

    private String id;              // Read id
    private Sequence sequence;      // Result of this.getSequence()

    // getters
    public String   getId()           { return this.id; }
    public Sequence getSequence()     { return this.sequence; }

	
    // setters
    public void setId( String id ) { this.id = id; }
    public void setSequence( Sequence sequence ) { this.sequence = sequence; }

    @Override
    public String toString() {
	// Returns the Ensembl/GeneCode corresponding header
	// >10
	return "@" +
	    this.getId() + "\n" +
	    StringSequence.toFasta( this.sequence );
    }
    
    // Constructors
    public Read() {
	this.id = "";
	this.sequence = null;
    }
}
