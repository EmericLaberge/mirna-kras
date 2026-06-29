package ca.iric.major.rinexus.rimaprisc;

import java.util.*;
import ca.iric.major.common.Duplex;

public class GenSeed {

    private static final char[] BASES = { 'A', 'C', 'G', 'U' };

    /** Mapping from seed index (0..6) to target index (0..8), or -1 for seed bulge. */
    private static final class Alignment {
        final int[] seedToTarget;  // length 7; -1 means that seed base is a bulge (unpaired)
        final boolean targetBulge; // true iff there is one bulged target nt in this alignment

        Alignment(int[] seedToTarget, boolean targetBulge) {
            this.seedToTarget = seedToTarget;
            this.targetBulge = targetBulge;
        }
    }

    private GenSeed() {
        // utility class
    }

    /**
     * Main entry point: generate all seeds for a given 9-nt target.
     */
    public static Set<String> generateSeedsForTarget9mer(String target9) {
        if (target9 == null || target9.length() != 9) {
            throw new IllegalArgumentException("target9 must be a 9-nt string");
        }

        target9 = normalizeNucleotides(target9);

        List<Alignment> alignments = buildAlignmentsFor9mer();
        Set<String> result = new LinkedHashSet<>();

        char[] seed = new char[7];
        for (Alignment alignment : alignments) {
            generateAllSeedsRecursive(seed, 0, target9, alignment, result);
        }

        return result;
    }

    /** Normalize to upper-case RNA alphabet (T -> U). */
    private static String normalizeNucleotides(String s) {
        StringBuilder sb = new StringBuilder(s.length());
        for (char c : s.toUpperCase().toCharArray()) {
            if (c == 'T') c = 'U';
            sb.append(c);
        }
        return sb.toString();
    }

    // -------------------------------------------------------------------------
    // Alignment generation
    // -------------------------------------------------------------------------

    /**
     * Build all alignments of a 7-nt seed onto a 9-nt target with:
     *  - no bulge (7 vs 7),
     *  - one seed bulge (7 seed nt to 6 target nt),
     *  - one target bulge (7 seed nt to 7 of 8 target nt).
     *
     * All indices refer to the target 9-mer positions [0..8].
     */

    private static List<Alignment> buildAlignmentsFor9mer() {
        List<Alignment> list = new ArrayList<>();

        // 1) No bulge: 7 seed nts onto 7-nt windows of target: [0..6], [1..7], [2..8]
        for (int tStart = 0; tStart <= 2; tStart++) {
            int[] map = new int[7];
            for (int i = 0; i < 7; i++) {
                map[i] = tStart + i;
            }
            list.add(new Alignment(map, false));
        }

	// 2) Seed bulge: one unpaired seed base, 6 target nts (windows [0..5], [1..6], [2..7], [3..8])
        for (int tStart = 0; tStart <= 3; tStart++) { // 6-nt windows
            for (int bulgeSeedIndex = 0; bulgeSeedIndex < 7; bulgeSeedIndex++) {
                int[] map = new int[7];
                for (int i = 0; i < 7; i++) {
                    if (i == bulgeSeedIndex) {
                        map[i] = -1; // seed bulge: unpaired
                    } else if (i < bulgeSeedIndex) {
                        map[i] = tStart + i;
                    } else {
                        // i > bulgeSeedIndex
                        map[i] = tStart + (i - 1);
                    }
                }
                list.add(new Alignment(map, false));
            }
        }

	// 3) Target bulge: one unpaired target base in an 8-nt window
        //    So we slide an 8-nt window: [0..7], [1..8] on the 9-mer
        //    and choose which of the 8 positions is bulged.
        for (int tStart = 0; tStart <= 1; tStart++) { // 8-nt windows
            for (int bulgeOffset = 0; bulgeOffset < 8; bulgeOffset++) {
                int[] map = new int[7];
                int tIdx = tStart;
                for (int i = 0; i < 7; i++) {
                    // Skip bulge position
                    if (tIdx == tStart + bulgeOffset) {
                        tIdx++;
                    }
                    map[i] = tIdx;
                    tIdx++;
                }
                list.add(new Alignment(map, true));
            }
        }

        return list;
    }

    // -------------------------------------------------------------------------
    // Seed enumeration and filtering
    // -------------------------------------------------------------------------
    
    private static void generateAllSeedsRecursive(
            char[] seed,
            int pos,
            String target9,
            Alignment alignment,
            Set<String> out
    ) {
        if (pos == 7) {
            if (isValidSeedForAlignment(seed, target9, alignment)) {
                out.add(new String(seed));
            }
            return;
        }

        for (char base : BASES) {
            seed[pos] = base;
            generateAllSeedsRecursive(seed, pos + 1, target9, alignment, out);
        }
    }

    /**
     * Check if this particular seed satisfies the binding rules for THIS alignment
     * against the given 9-nt target.
     */
    private static boolean isValidSeedForAlignment(char[] seed, String target9, Alignment alignment) {
        // Arrays indexed 0..6 for g2..g8
        boolean[] wc = new boolean[7]; // WC pair at this seed position?
        boolean[] paired = new boolean[7]; // is this seed base paired (not a bulge)?

	int wcCount = 0;
        int mismatchCount = 0;
        int seedBulgeCount = 0;

        for (int i = 0; i < 7; i++) {
            int tIdx = alignment.seedToTarget[i];
            if (tIdx < 0) {
                // Seed bulge
                seedBulgeCount++;
                paired[i] = false;
                wc[i] = false;
            } else {
                paired[i] = true;
                char sBase = seed[i];
                char tBase = target9.charAt(tIdx);
                boolean isWC = isWatsonCrick(sBase, tBase);
                wc[i] = isWC;
                if (isWC) {
                    wcCount++;
                } else {
                    mismatchCount++;
                }
            }
        }

	// At most one seed bulge in any alignment by construction,
        // but code will tolerate bad alignments as well.
        int targetBulgeCount = alignment.targetBulge ? 1 : 0;
        int deviationCount = mismatchCount + seedBulgeCount + targetBulgeCount;

        // ---- Rule A: >= 4 contiguous WC pairs between g2..g7 (indices 0..5) ----
        if (hasContiguous4MerInG2toG7(wc, paired)) {
            return true;
        }

        // ---- Rule B: >= 6 WC pairs and exactly 1 deviation (mismatch or bulge) ----
        if (wcCount >= 6 && deviationCount == 1) {
            return true;
        }

        return false;
    }

    /** Watson-Crick pair: A-U, U-A, C-G, G-C. */
    private static boolean isWatsonCrick(char seedBase, char targetBase) {
        seedBase = Character.toUpperCase(seedBase);
        targetBase = Character.toUpperCase(targetBase);
        if (seedBase == 'T') seedBase = 'U';
        if (targetBase == 'T') targetBase = 'U';

        return (seedBase == 'A' && targetBase == 'U') ||
               (seedBase == 'U' && targetBase == 'A') ||
               (seedBase == 'C' && targetBase == 'G') ||
               (seedBase == 'G' && targetBase == 'C');
    }

    /**
     * Check if there exists a contiguous run of 4 WC pairs between g2..g7
     * seed indices 0..5 (so windows [0..3], [1..4], [2..5]).
     * Bulges or mismatches break contiguity.
     */
    private static boolean hasContiguous4MerInG2toG7(boolean[] wc, boolean[] paired) {
        // Possible start positions for a 4-mer within indices 0..5: 0,1,2
        for (int start = 0; start <= 2; start++) {
            boolean ok = true;
            for (int k = 0; k < 4; k++) {
                int idx = start + k;
                if (!paired[idx] || !wc[idx]) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                return true;
            }
        }
        return false;
    }

    // -------------------------------------------------------------------------
    // Small demo main (optional)
    // -------------------------------------------------------------------------
    public static void main(String[] args) {
        if (args.length != 1) {
            System.err.println("Usage: GenSeed <9-mer target>");
            System.exit(1);
        }
        String target = args[0];
        Set<String> seeds = generateSeedsForTarget9mer(target);
        System.out.println("Target: " + target);
        System.out.println("Number of compatible seeds: " + seeds.size());
        for (String s : seeds) {
            System.out.println(s);
	    //Duplex d = new Duplex( target, s, "", 2, false, false );
	    //System.out.println( d );
        }
    }
}
