/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.ArrayList;

public class Positions {
    // attributes
    ArrayList<Integer> positions = new ArrayList<>();
    int count = 0;

    public Positions( Integer p ) { this.positions.add( p ); this.count++; }

    // getters
    public int getCount() { return this.count; }
    public ArrayList<Integer> getPositions() { return this.positions; }

    // setters
    public void addPosition( Integer p ) { this.positions.add( p ); this.count++; }
}
