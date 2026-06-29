/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */

package ca.iric.major.rinexus.rimaprisc;

import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.List;

// for iric installation
import lombok.Getter;                                                                                                                          
import lombok.NoArgsConstructor;                                                                                                               
import lombok.Setter;

@Getter
@Setter
@NoArgsConstructor // ensures Jackson can deserialize empty JSON requests
public class RISCScanRequest {
    private String gencodeVersion = "v46";
    private String transcriptName;
    private String miRName;
    private String minSeed;
    private double maxKd;
    private boolean UTR5;
    private boolean CDS;
    private boolean UTR3;
    private boolean supplementaryPaired;
    private double seedAccessible;
    private double supplementaryAccessible;
    private int conserved;

    private String miRNATargetPair;
    private String bindingLandscape;
    private List<Integer> transcriptRegions;
    private List<Integer> bindingPositions;
    private List<Double> transcriptAccessibility;

    // Getters and setters
    @JsonProperty("miRName") // Ensure correct JSON property mapping
    public String getMiRName() {
        return this.miRName;
    }

    public void setMiRName(String miRName) {
        this.miRName = miRName;
    }

    @JsonProperty("gencode") // Ensure correct JSON property mapping
    public String getGencodeVersion() {
        return this.gencodeVersion;
    }

    public void setGencodeVersion(String gencodeVersion) {
        this.gencodeVersion = gencodeVersion;
    }

    @JsonProperty("minSeed") // Ensure correct JSON property mapping
    public String getMinSeed() {
	return this.minSeed;
    }

    public void setMinSeed( String minSeed ) {
	this.minSeed = minSeed;
    }

    @JsonProperty("maxKd") // Ensure correct JSON property mapping
    public double getMaxKd() {
	return this.maxKd;
    }

    public void setMaxKd( double maxKd ) {
	this.maxKd = maxKd;
    }
    
    @JsonProperty("UTR5") // Ensure correct JSON property mapping
    public boolean isUTR5() {
	return this.UTR5;
    }

    public void setUTR5( boolean UTR5 ) {
	this.UTR5 = UTR5;
    }

    @JsonProperty("CDS") // Ensure correct JSON property mapping
    public boolean isCDS() {
	return this.CDS;
    }

    public void setCDS( boolean CDS ) {
	this.CDS = CDS;
    }

    @JsonProperty("UTR3") // Ensure correct JSON property mapping
    public boolean isUTR3() {
	return this.UTR3;
    }

    public void setUTR3( boolean UTR3 ) {
	this.UTR3 = UTR3;
    }

    @JsonProperty("transcriptName") // Ensure correct JSON property mapping
    public String getTranscriptName() {
        return this.transcriptName;
    }

    public void setTranscriptName( String transcriptName ) {
        this.transcriptName = transcriptName;
    }

    @JsonProperty("supplementaryPaired") // Ensure correct JSON property mapping
    public boolean isSupplementaryPaired() {
        return this.supplementaryPaired;
    }

    public void setSupplementaryPaired( boolean supplementaryPaired ) {
	this.supplementaryPaired = supplementaryPaired;
    }

    @JsonProperty("seedAccessible") // Ensure correct JSON property mapping
    public double getSeedAccessible() {
        return this.seedAccessible;
    }

    public void setSeedAccessible( double seedAccessible ) {
	this.seedAccessible = seedAccessible;
    }

    @JsonProperty("supplementaryAccessible") // Ensure correct JSON property mapping
    public double getSupplementaryAccessible() {
        return this.supplementaryAccessible;
    }

    public void setSupplementaryAccessible( double supplementaryAccessible ) {
	this.supplementaryAccessible = supplementaryAccessible;
    }

    @JsonProperty("conserved") // Ensure correct JSON property mapping
    public int getConserved() {
        return this.conserved;
    }

    public void setConserved( int conserved ) {
	this.conserved = conserved;
    }

    public String getMiRNATargetPair() {
        return this.miRNATargetPair;
    }

    public void setMiRNATargetPair( String miRNATargetPair ) {
        this.miRNATargetPair = miRNATargetPair;
    }

    public String getBindingLandscape() {
        return this.bindingLandscape;
    }

    public void setBindingLandscape( String bindingLandscape ) {
        this.bindingLandscape = bindingLandscape;
    }
    
    public List<Integer> getTranscriptRegions() {
        return transcriptRegions;
    }

    public void setTranscriptRegions(List<Integer> transcriptRegions) {
        this.transcriptRegions = transcriptRegions;
    }

    public List<Integer> getBindingPositions() {
        return this.bindingPositions;
    }

    public void setBindingPositions(List<Integer> bindingPositions) {
        this.bindingPositions = bindingPositions;
    }

    public List<Double> getTranscriptAccessibility() {
	return this.transcriptAccessibility;
    }

    public void setTranscriptAccessibility(List<Double> transcriptAccessibility) {
        this.transcriptAccessibility = transcriptAccessibility;
    }
}
