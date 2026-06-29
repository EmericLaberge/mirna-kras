package ca.iric.major.tools;

public class SeedDuplexCache {
    // Simple 2-bit encoding: A=0, C=1, G=2, U=3
    // 8mer = 16 bits fits in an int
    public static String decode8mer(int encoded) {
        StringBuilder sb = new StringBuilder(8);
        for (int i = 7; i >= 0; i--) {
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
