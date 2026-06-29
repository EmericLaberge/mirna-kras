/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Objects;

public class Loop implements Comparable<Loop> {
    protected final int position;
    protected final int length;
    
    public Loop( int position, int length ) {
	this.position = position;
	this.length = length;
    }

    @Override
    public int compareTo( Loop other ) {
	return Integer.compare( this.position, other.position );
    }

    @Override
    public boolean equals( Object o ) { // assume o is a Loop; discipline the users
        Loop l = (Loop)o;
        return this.position == l.position && this.length == l.length;
    }

    @Override
    public int hashCode() { return Objects.hash( position, length ); }


    @Override
    public String toString() {
	return "i j = " + (this.position + 1) + " " + (this.position + this.length) + " : ";
    }
}
