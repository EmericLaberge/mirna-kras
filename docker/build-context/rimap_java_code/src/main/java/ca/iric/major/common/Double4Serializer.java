/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import com.fasterxml.jackson.core.JsonGenerator;
import com.fasterxml.jackson.databind.JsonSerializer;
import com.fasterxml.jackson.databind.SerializerProvider;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;

public class Double4Serializer extends JsonSerializer<Double> {
    @Override
    public void serialize( Double value, JsonGenerator gen, SerializerProvider serializers ) throws IOException {
        if( value == null ) { gen.writeNull(); return; }
        // Exactly 4 decimals:
        // BigDecimal bd = BigDecimal.valueOf(value).setScale(4, RoundingMode.HALF_UP);
        // Up to 4 decimals (no unnecessary trailing zeros):
        BigDecimal bd = BigDecimal.valueOf( value ).setScale( 4, RoundingMode.HALF_UP ).stripTrailingZeros();
        gen.writeNumber( bd );
    }
}
