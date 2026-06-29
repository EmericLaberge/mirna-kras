/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.IOException;
import java.io.File;
import java.io.ObjectInputStream;
import java.io.FileInputStream;
import java.nio.file.*;
import java.util.Map;
import java.util.HashMap;

import java.util.stream.Collectors;

public final class DisturbanceManager {

    private String disturbanceFile;
    private Map<String,Double> kmerDisturbance = new HashMap<>();

    public DisturbanceManager( String disturbanceFile ) {
	this.disturbanceFile = disturbanceFile;
	try( ObjectInputStream ois = new ObjectInputStream( new FileInputStream( disturbanceFile ) ) ) {
	    kmerDisturbance = (Map<String,Double>) ois.readObject();
	} catch( IOException | ClassNotFoundException e ) {
	    e.printStackTrace();
	}
    }

    public double getDisturbance( String guide ) {
	String sevenMer = guide.substring( 1, 8 ); // g1-g8
	String sixMer   = guide.substring( 1, 7 ); // g1-g7
	String fiveMer  = guide.substring( 1, 6 ); // g1-g6
	return ( kmerDisturbance.get( sevenMer ) + kmerDisturbance.get( sixMer ) / 10 + kmerDisturbance.get( fiveMer ) / 100 ) / 3.0;
    }
}
