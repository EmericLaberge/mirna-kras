/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

/**
 * Interface for reporter transcripts (composed of a Sequence)
 *
 * @version 1.0
 * @author Francois Major
 * @copyright 2024 - MajorLab, IRIC, Universite de Montreal
 * @license MIT
*/

import java.lang.IllegalArgumentException;
import java.util.List;

public class ReporterTranscript extends CodingTranscript implements Transcript {

    // Constructor
    public ReporterTranscript() {
	super();
    }
}
