package ca.iric.major.tools;

public class SeedToEightMerCache {
    // Simple 2-bit encoding: A=0, C=1, G=2, U=3
    // 7mer = 14 bits fits in an int
    public static int encode7mer(String seed) {
        int encoded = 0;
        for (char c : seed.toCharArray()) {
            encoded <<= 2;
            switch (c) {
                case 'A': encoded |= 0; break;
                case 'C': encoded |= 1; break;
                case 'G': encoded |= 2; break;
                case 'U': encoded |= 3; break;
            }
        }
        return encoded;
    }

    public static String decode7mer(int encoded) {
        StringBuilder sb = new StringBuilder(7);
        for (int i = 6; i >= 0; i--) {
            int bits = (encoded >> (i * 2)) & 3;
            switch (bits) {
                case 0: sb.append('A'); break;
                case 1: sb.append('C'); break;
                case 2: sb.append('G'); break;
                case 3: sb.append('U'); break;
            }
        }
        return sb.toString();
    }
}
