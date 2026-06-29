/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.LinkedHashSet;

import java.util.Arrays;
import java.util.Collections;
import java.util.stream.Collectors;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.lang.InterruptedException;

/**
 * Design instances are used to generate and keep the RNA guides associated to a given set of targets.
 */
public class Design {

    // Attributes
    Set<ProteinCodingTranscript> targets;
    KMerMap seedKmers;
    KMerMap suppKmers;
    GuideMap guides;
    
    // Constructor
    public Design( Set<ProteinCodingTranscript> targets, KMerMap seedKmers, KMerMap suppKmers, GuideMap guides ) {
	this.targets = targets;
	this.seedKmers = seedKmers;
	this.suppKmers = suppKmers;
	this.guides = guides;
    }

    // public methods
    public Set<ProteinCodingTranscript> getTargets() { return this.targets; }

    public String toString() {
	return this.targets.toString();
    }
}
