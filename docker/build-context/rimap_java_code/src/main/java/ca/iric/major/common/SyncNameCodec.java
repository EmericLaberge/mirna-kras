package ca.iric.major.common;

public final class SyncNameCodec {

    private SyncNameCodec() {}

    // Map A/C/G/U(T) -> 2-bit code
    private static int baseToCode(char b) {
        switch (Character.toUpperCase(b)) {
            case 'A': return 0b00;
            case 'C': return 0b01;
            case 'G': return 0b10;
            case 'U':
            case 'T': return 0b11;  // treat T as U
            default:
                throw new IllegalArgumentException("Invalid base: " + b);
        }
    }

    // Map 2-bit code -> base (we'll always emit U, not T)
    private static char codeToBase(int code) {
        switch (code & 0b11) {
            case 0b00: return 'A';
            case 0b01: return 'C';
            case 0b10: return 'G';
            case 0b11: return 'U';
            default:   throw new IllegalArgumentException("Invalid code: " + code);
        }
    }

    /**
     * Encode a mature SYNC guide as:
     *   sync-SEED-BODY-L
     * where:
     *   SEED = g2–g8 (7 nt, letters)
     *   BODY = hex for [g1, g9..g_end] (2 bits per nt)
     *   L    = #nts after g17 (1..9)
     */
    public static String encodeSyncName(String guide) {
        if (guide == null) {
            throw new IllegalArgumentException("Guide cannot be null");
        }
        guide = guide.trim().toUpperCase();
        int len = guide.length();
        if (len < 18 || len > 26) {
            throw new IllegalArgumentException("Guide length must be between 18 and 26, got " + len);
        }

        // SEED = g2–g8
        String seed = guide.substring(1, 8);  // indices 1..7

        // Tail length L = # nts after g17
        int tailLen = len - 17;
        if (tailLen < 1 || tailLen > 9) {
            throw new IllegalArgumentException("Tail length must be between 1 and 9, got " + tailLen);
        }

        // BODY encodes: g1 + g9..g_end
        // g1 at index 0, g9 at index 8
        long bodyVal = 0L;

        // First g1
        bodyVal = (bodyVal << 2) | baseToCode(guide.charAt(0));

        // Then g9..g_end
        for (int i = 8; i < len; i++) {
            bodyVal = (bodyVal << 2) | baseToCode(guide.charAt(i));
        }

        String bodyHex = Long.toHexString(bodyVal).toUpperCase();

        return "sync-" + seed + "-" + bodyHex + "-" + tailLen;
    }

    /**
     * Decode sync-SEED-BODY-L back into the mature sequence.
     * Not strictly needed for your tools (since you start from sequence),
     * but useful for sanity checks.
     */
    public static String decodeSyncName(String syncName) {
        if (syncName == null || !syncName.startsWith("sync-")) {
            throw new IllegalArgumentException("Invalid SYNC name: " + syncName);
        }

        String[] parts = syncName.split("-");
        if (parts.length != 4) {
            throw new IllegalArgumentException("SYNC name must have 4 parts: sync-SEED-BODY-L, got: " + syncName);
        }

        String seed  = parts[1]; // g2–g8
        String body  = parts[2]; // hex for g1 + g9..g_end
        String tailS = parts[3]; // L (tail length)

        if (seed.length() != 7) {
            throw new IllegalArgumentException("SEED must be 7 nt, got: " + seed);
        }

        int tailLen;
        try {
            tailLen = Integer.parseInt(tailS);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Tail length L must be integer 1..9, got: " + tailS, e);
        }
        if (tailLen < 1 || tailLen > 9) {
            throw new IllegalArgumentException("Tail length must be between 1 and 9, got: " + tailLen);
        }

        int len = 17 + tailLen;              // total guide length
        int bodyNtCount = 1 + (len - 8);     // g1 + (g9..g_end) = 10 + tailLen

        long bodyVal = body.isEmpty() ? 0L : Long.parseLong(body, 16);

        // Decode from LSB to MSB
        char[] bodyBases = new char[bodyNtCount];
        for (int i = bodyNtCount - 1; i >= 0; i--) {
            int code = (int) (bodyVal & 0b11);
            bodyBases[i] = codeToBase(code);
            bodyVal >>= 2;
        }

        // Assemble:
        // bodyBases[0]     -> g1
        // seed (7 nt)      -> g2..g8
        // bodyBases[1..9]  -> g9..g17
        // bodyBases[10..]  -> tail (g18..g_end)
        StringBuilder sb = new StringBuilder(len);
        sb.append(bodyBases[0]);          // g1
        sb.append(seed);                  // g2..g8
        sb.append(bodyBases, 1, 9);       // g9..g17
        sb.append(bodyBases, 10, tailLen);// g18..g_end

        return sb.toString();
    }
}
