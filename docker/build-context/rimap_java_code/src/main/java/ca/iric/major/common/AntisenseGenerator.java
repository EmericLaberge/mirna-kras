/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.common;

/**
 * Antisense genertor produces a String generator for antisense with 4-nt and 3-nt grips for, respectively, seed and supp.
 *    this generator is lazy and cheap (does not generate all possible sequences), keeps one sequence at a time.
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 1.0 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Predicate;

/**
 * AntisenseGenerator
 *
 * Generates antisense (guide) RNA sequences of length 21 nt with the following layout:
 *
 *  g1         = C
 *  g2–g5      = 4-nt seed grip (already reverse-complemented and passed in as argument)
 *  g6–g8      = fully sampled (A,C,G,U)
 *  g9–g12     = fixed central region "UUAA"
 *  g13–g15    = 3-nt supplementary grip (already reverse-complemented and passed as argument)
 *  g16–g17    = fully sampled (A,C,G,U)
 *  g18–g21    = fixed tail "AAUU"
 *
 * Free positions: g6,g7,g8,g16,g17 (5 positions, 4 nucleotides each) → 4^5 = 1024 sequences max.
 *
 * The generator is lazy: it only constructs the next sequence on demand, and at most one
 * candidate is stored in memory at a time. A Predicate<String> filter is applied on each
 * candidate; only sequences passing the filter are returned.
 */

public class AntisenseGenerator implements Iterator<String> {

    private static final int SEQUENCE_LENGTH = 21;
    private static final int MAX_COMBINATIONS = 1024; // 4^5

    private static final char[] NT = { 'A', 'C', 'G', 'U' };

    private final String seedGrip;   // length 4 → g2–g5
    private final String suppGrip;   // length 3 → g13–g15
    private final Predicate<String> filter;

    // Enumerates all combinations for g6,g7,g8,g16,g17
    private int index = 0;

    // Lazy state
    private String nextSequence = null;
    private boolean prepared = false;

    public AntisenseGenerator(String seedGrip, String suppGrip, Predicate<String> filter) {
        if (seedGrip == null || seedGrip.length() != 4) {
            throw new IllegalArgumentException("seedGrip must be a 4-nt string");
        }
        if (suppGrip == null || suppGrip.length() != 3) {
            throw new IllegalArgumentException("suppGrip must be a 3-nt string");
        }
        this.seedGrip = seedGrip.toUpperCase();
        this.suppGrip = suppGrip.toUpperCase();
        this.filter = (filter != null) ? filter : s -> true;
    }

    @Override
    public boolean hasNext() {
        if (!prepared) {
            prepareNext();
        }
        return nextSequence != null;
    }

    @Override
    public String next() {
        if (!hasNext()) {
            throw new NoSuchElementException("AntisenseGenerator exhausted");
        }
        String result = nextSequence;
        // Reset lazy state so that the next call to hasNext() will advance enumeration
        nextSequence = null;
        prepared = false;
        return result;
    }

    private void prepareNext() {
        prepared = true; // we are about to try to find the next; if we don't, nextSequence stays null

        while (index < MAX_COMBINATIONS) {
            int current = index++;
            // Decode "current" in base 4 into 5 digits corresponding to the 5 free positions
            int d0 = current % 4; current /= 4;   // g6
            int d1 = current % 4; current /= 4;   // g7
            int d2 = current % 4; current /= 4;   // g8
            int d3 = current % 4; current /= 4;   // g16
            int d4 = current % 4;                 // g17

            char[] seq = new char[SEQUENCE_LENGTH];

            // g1 = C (index 0)
            seq[0] = 'C';

            // g2–g5 = seedGrip (indices 1–4)
            for (int i = 0; i < 4; i++) {
                seq[1 + i] = seedGrip.charAt(i);
            }

            // g6–g8 sampled (indices 5–7)
            seq[5] = NT[d0];
            seq[6] = NT[d1];
            seq[7] = NT[d2];

            // g9–g12 = "UUAA" (indices 8–11)
            seq[8]  = 'U';
            seq[9]  = 'U';
            seq[10] = 'A';
            seq[11] = 'A';

            // g13–g15 = suppGrip (indices 12–14)
            for (int i = 0; i < 3; i++) {
                seq[12 + i] = suppGrip.charAt(i);
            }

            // g16–g17 sampled (indices 15–16)
            seq[15] = NT[d3];
            seq[16] = NT[d4];

            // g18–g21 = "AAUU" (indices 17–20)
            seq[17] = 'A';
            seq[18] = 'A';
            seq[19] = 'U';
            seq[20] = 'U';

            String candidate = new String(seq);

            // Apply external filter
            if (filter.test(candidate)) {
                nextSequence = candidate;
                return;
            }
            // else: continue to next index
        }

        // No more candidates
        nextSequence = null;
    }
}
