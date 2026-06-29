/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Objects;

public class Site implements Comparable<Site> {
    protected final int position;
    protected final int length;
    
    public Site( int position, int length ) {
	this.position = position;
	this.length = length;
    }

    public int getStart() { return this.position; }
    public int getLength() { return this.length; }
    public int getEnd() { return this.position + this.length - 1; }

    @Override
    public int compareTo( Site other ) {
	return Integer.compare( this.position, other.position );
    }

    // check if the site is accessible (4 dots in a row) in a given conformational state
    public boolean accessible( String state ) {
	int minExp = Math.min( this.length, 4 );
	String subseq = state.substring( this.position, this.position + this.length );

	int consecutive = 0;
	for( int i = 0; i < subseq.length(); i++ ) {
	    if( subseq.charAt(i) == '.' ) {
		consecutive++;
		if( consecutive >= minExp ) {
		    return true; // found 4 consecutive dots
		}
	    } else {
		consecutive = 0; // reset streak
	    }
	}
	return false;
    }

    @Override
    public boolean equals( Object o ) { // assume o is a Site; discipline the users
        Site l = (Site)o;
        return this.position == l.position && this.length == l.length;
    }

    @Override
    public int hashCode() { return Objects.hash( position, length ); }


    @Override
    public String toString() {
	return "Site( " + (this.position + 1) + " " + (this.position + this.length) + " )";
    }
}
