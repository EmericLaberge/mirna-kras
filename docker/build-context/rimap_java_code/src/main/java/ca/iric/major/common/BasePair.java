/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Objects;

public final class BasePair implements Comparable<BasePair> {
    protected final int i;
    protected final int j;

    public BasePair( int i, int j ) { // assume i < j for consistency; discipline the users
	this.i = i;
	this.j = j;
    }

    public int getI() { return this.i; }
    public int getJ() { return this.j; }

    @Override
    public int compareTo( BasePair other ) {
	return Integer.compare( this.i, other.i );
    }

    @Override
    public boolean equals( Object o ) { // assume o is a BasePair; discipline the users
        BasePair bp = (BasePair)o;
        return this.i == bp.i && this.j == bp.j;
    }

    @Override
    public int hashCode() { return Objects.hash( i, j ); }

    @Override
    public String toString() {
        return "i j = " + (this.i+1) + " " + (this.j+1) + " : ";
    }
}
