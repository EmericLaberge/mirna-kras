/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.JsonSerializer;
import com.fasterxml.jackson.databind.SerializerProvider;

import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;

public class TwoDecimalDoubleSerializer extends JsonSerializer<Double> {

    private static final DecimalFormat DF;
    static {
        DF = new DecimalFormat( "#0.00" );
        DF.setRoundingMode( RoundingMode.HALF_UP );
    }

    @Override
    public void serialize( Double value, JsonGenerator gen, SerializerProvider serializers )
            throws IOException {
        gen.writeNumber(DF.format( value ) );
    }
}
