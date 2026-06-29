/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonFormat;
import com.fasterxml.jackson.databind.annotation.JsonSerialize;

public class JsonInteraction {

    @JsonProperty("gna")
    private String guideRNA;

    @JsonProperty("gseq")
    private String guideRNASequence;

    @JsonProperty("brl")
    private int bridgeLength;

    @JsonProperty("pst")
    private int positionStart;

    @JsonProperty("pen")
    private int positionEnd;
    
    @JsonProperty("alt")
    private String alignedTarget;
    
    @JsonProperty("alp")
    private String alignedPairs;
    
    @JsonProperty("alg")
    private String alignedGuide;
    
    @JsonProperty("pcs")
    private double phastconsScore;

    @JsonProperty("sac")
    @JsonSerialize(using = Double4Serializer.class)
    private double seedAccessibility;

    @JsonProperty("suc")
    @JsonSerialize(using = Double4Serializer.class)
    private double suppAccessibility;
    
    @JsonProperty("rgn")
    private String region;
    
    @JsonProperty("isc")
    private int seedScore;
    
    @JsonProperty("ity")
    private String seedType;

    @JsonProperty("sty")
    private String siteType;
    
    @JsonProperty("sen")
    private double seedEnergy;
    
    @JsonProperty("sut")
    private String supplementaryType;
    
    @JsonProperty("sue")
    private double supplementaryEnergy;
    
    @JsonProperty("tid")
    private String transcriptId;
    
    @JsonProperty("kdv")
    private double Kd;

    @JsonProperty("gid")
    private String guideId;
    
    @JsonProperty("trn")
    private String transcript;
    
    @JsonProperty("tseq")
    private String transcriptSequence;
    
    @JsonProperty("den")
    private double duplexEnergy;

    @JsonProperty("itb")
    private boolean inTarBase;

    // Getters
    public String getGuideRNA() { return this.guideRNA; }
    public String getGuideRNASequence() { return this.guideRNASequence; }
    public int getBridgeLength() { return this.bridgeLength; }
    public int getPositionStart() { return this.positionStart; }
    public int getPositionEnd() { return this.positionEnd; }   
    public String getAlignedTarget() { return this.alignedTarget; }
    public String getAlignedPairs() { return this.alignedPairs; }
    public String getAlignedGuide() { return this.alignedGuide; }
    public double getPhastconsScore() { return this.phastconsScore; }
    public double getSeedAccessibility() { return seedAccessibility; }
    public double getSuppAccessibility() { return suppAccessibility; }
    public String getRegion() { return this.region; }
    public int getSeedScore() { return this.seedScore; }
    public String getSeedType() { return this.seedType; }
    public String getSiteType() { return this.siteType; }
    public double getSeedEnergy() { return this.seedEnergy; }
    public String getSupplementaryType() { return this.supplementaryType; }
    public double getSupplementaryEnergy() { return this.supplementaryEnergy; }
    public String getTranscriptId() { return this.transcriptId; }
    @JsonProperty("kdv")
    public double getKd() { return this.Kd; }
    public String getGuideId() { return this.guideId; }
    public String getTranscript() { return this.transcript; }
    public String getTranscriptSequence() { return this.transcriptSequence; }
    public double getDuplexEnergy() { return this.duplexEnergy; }
    public boolean getInTarBase() { return this.inTarBase; }

    public boolean inRegion( int region ) { // region int values defined in common/CodingTranscript
	if( region == 7 ) return true;
	if( this.region.equals( "5UTR" ) )
	    return region == 4 || region == 5 || region == 6;
	if( this.region.equals( "CDS" ) )
	    return region == 2 || region == 3 || region == 6;
	if( this.region.equals( "3UTR" ) )
	    return region == 1 || region == 3 || region == 5;
	return false;
    }
    
    // Setters
    public void setGuideRNA( String x ) { this.guideRNA = x; }
    public void setGuideRNASequence( String x ) { this.guideRNASequence = x; }
    public void setBridgeLength( int x ) { this.bridgeLength = x; }
    public void setPositionStart( int x ) { this.positionStart = x; }
    public void setPositionEnd( int x ) { this.positionEnd = x; }
    public void setAlignedTarget( String x ) { this.alignedTarget = x; }
    public void setAlignedPairs( String x ) { this.alignedPairs = x; }
    public void setAlignedGuide( String x ) { this.alignedGuide = x; }
    public void setPhastconsScore( double x ) { this.phastconsScore = x; }
    public void setSeedAccessibility( double x ) { this.seedAccessibility = x; }
    public void setSuppAccessibility( double x ) { this.suppAccessibility = x; }
    public void setRegion( String x ) { this.region = x; }
    public void setSeedScore( int x ) { this.seedScore = x; }
    public void setSeedType( String x ) { this.seedType = x; }
    public void setSiteType( String x ) { this.siteType = x; }
    public void setSeedEnergy( double x ) { this.seedEnergy = x; }
    public void setSupplementaryType( String x ) { this.supplementaryType = x; }
    public void setSupplementaryEnergy( double x ) { this.supplementaryEnergy = x; }
    public void setTranscriptId( String x ) { this.transcriptId = x; }
    public void setKd( double x ) { this.Kd = x; }
    public void setGuideId( String x ) { this.guideId = x; }
    public void setTranscript( String x ) { this.transcript = x; }
    public void setTranscriptSequence( String x ) { this.transcriptSequence = x; }
    public void setDuplexEnergy( double x ) { this.duplexEnergy = x; }
    public void setInTarBase( boolean inTarBase ) { this.inTarBase = inTarBase; }

    @Override
    public String toString() {
        return
	    "5'-" + this.getAlignedTarget() + "-3' " + this.getTranscript() + " [" + this.getPositionStart() + ".." + this.getPositionEnd() + "] " + "(" + this.getRegion() + ") target site: " + this.getSiteType() + (getInTarBase() ? " in: TarBase-v9.0" : "") + "\n" +
	    "   " + this.getAlignedPairs() + "    " + this.getSeedType() + "(" + String.format("%,.2f", this.getSeedAccessibility() ) + ") " + this.getSupplementaryType() + "(" + String.format("%,.2f", this.getSuppAccessibility() ) + ") bridge: " + this.getBridgeLength() + " phastcons: " + String.format("%,.2f", this.getPhastconsScore() ) + " dG: " + String.format("%,.2f", this.getDuplexEnergy() ) + " Kd: " + String.format( "%,.0f", this.getKd() ) + "\n" +
	    "3'-" + this.getAlignedGuide() + "-5' " + this.getGuideRNA();
    }
}
