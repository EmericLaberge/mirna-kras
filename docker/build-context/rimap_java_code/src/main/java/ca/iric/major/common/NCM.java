/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Set;

public class NCM {
    String innerSequence;
    Set<BasePair> basePairs;
    String structure;
    
    public NCM( String innerSequence, Set<BasePair> basePairs, String structure ) {
        this.innerSequence = innerSequence;
        this.basePairs = basePairs;
        this.structure = structure;
    }
    
    public String toString() {
        return "(" + innerSequence + ", " + basePairs + ", \"" + structure + "\")";
    }
}
