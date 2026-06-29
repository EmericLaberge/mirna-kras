/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Json class for guide RNAs
 *     
 * @version %I% %G%
 * @author Francois Major
 * @copyright 2025 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import com.fasterxml.jackson.annotation.JsonProperty;

// manage json information of a guide
public class JsonGuide {

    @JsonProperty("cna")
    private String commonName;  // e.g., hsa-miR-101-5p
    @JsonProperty("mid")
    private String id;          // e.g., MIMAT0004513
    @JsonProperty("spc")
    private String species;     // e.g., Homo sapiens
    @JsonProperty("seq")
    private String sequence;    // e.g., UACAGUACUGUGAUAACUGAA
    @JsonProperty("dbs")
    private double disturbance; // e.g., disturbance score 0 to 1.

    public JsonGuide() {
        // Required for JSON deserialization
    }

    public JsonGuide( String commonName, String id, String species, String sequence, double disturbance ) {
        this.commonName = commonName;
        this.id = id;
        this.species = species;
        this.sequence = sequence;
	this.disturbance = disturbance;
    }

    public String getCommonName() {
        return this.commonName;
    }

    public void setCommonName( String commonName ) {
        this.commonName = commonName;
    }

    public String getId() {
        return this.id;
    }

    public void setId( String id ) {
        this.id = id;
    }

    public String getSpecies() {
        return this.species;
    }

    public void setSpecies( String species ) {
        this.species = species;
    }

    public String getSequence() {
        return this.sequence;
    }

    public void setSequence( String sequence ) {
        this.sequence = sequence;
    }

    public double getDisturbance() {
        return this.disturbance;
    }

    public void setDisturbance( double disturbance ) {
        this.disturbance = disturbance;
    }

    @Override
    public String toString() {
        return "JsonGuide{" +
                "cna='" + this.commonName + '\'' +
                ", mid='" + this.id + '\'' +
                ", spc='" + this.species + '\'' +
                ", seq='" + this.sequence + '\'' +
	        ", dbs=" + this.disturbance +
                '}';
    }
}
