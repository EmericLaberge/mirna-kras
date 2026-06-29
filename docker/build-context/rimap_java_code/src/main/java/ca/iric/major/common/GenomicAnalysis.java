/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.Collectors;
import org.apache.commons.csv.*;

public class GenomicAnalysis {

    // Function to extract genomic location from GTF file
    public static Map<String, String> extractGenomicLocation(String gtfFile, String enstId) throws IOException {
        Reader in = new FileReader(gtfFile);
        Iterable<CSVRecord> records = CSVFormat.TDF
                .withCommentMarker('#')
                .withHeader("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
                .parse(in);

        for (CSVRecord record : records) {
            String feature = record.get("feature");
            String attribute = record.get("attribute");

            if (feature.equals("transcript") && attribute.contains("transcript_id \"" + enstId + "\"")) {
                String seqname = record.get("seqname");
                String start = record.get("start");
                String end = record.get("end");
                String strand = record.get("strand");

                Map<String, String> location = new HashMap<>();
                location.put("seqname", seqname);
                location.put("start", start);
                location.put("end", end);
                location.put("strand", strand);
                return location;
            }
        }

        System.out.println("ENST ID " + enstId + " not found in GTF file.");
        return null;
    }

    // Function to fetch sequence using bedtools
    public static void fetchSequence(String fastaFile, String seqname, int start, int end, String outputFile) throws IOException, InterruptedException {
        String bedContent = seqname + "\t" + (start - 1) + "\t" + end + "\n";
        Path bedFile = Paths.get("region.bed");

        // Write the BED content to a file
        Files.write(bedFile, bedContent.getBytes());

        // Run bedtools getfasta
        ProcessBuilder pb = new ProcessBuilder("bedtools", "getfasta", "-fi", fastaFile, "-bed", bedFile.toString(), "-fo", outputFile);
        pb.inheritIO();
        Process process = pb.start();
        int exitCode = process.waitFor();

        if (exitCode == 0) {
            System.out.println("Sequence extraction completed. Output is in " + outputFile);
        } else {
            System.out.println("Error running bedtools.");
        }
    }

    // Function to find 31mer in the extracted sequence and adjust coordinates
    public static Map<String, Integer> search31merAndAdjustCoordinates(String fastaFile, int start, String query31mer) throws IOException {
        List<String> lines = Files.readAllLines(Paths.get(fastaFile));
        String sequence = lines.stream().skip(1).collect(Collectors.joining()).replaceAll("\\n", "");

        int index = sequence.indexOf(query31mer);
        if (index == -1) {
            System.out.println("31mer not found in the extracted sequence.");
            return null;
        }

        int adjustedStart = start + index;
        int adjustedEnd = adjustedStart + query31mer.length();
        Map<String, Integer> coordinates = new HashMap<>();
        coordinates.put("start", adjustedStart);
        coordinates.put("end", adjustedEnd);
        return coordinates;
    }

    // Function to fetch PhastCons scores from BedGraph file
    public static double[] fetchPhastConsScores(String bedGraphFile, String seqname, int start, int end) throws IOException {
        double[] scores = new double[end - start];
        int index = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(bedGraphFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields[0].equals(seqname)) {
                    int regionStart = Integer.parseInt(fields[1]);
                    int regionEnd = Integer.parseInt(fields[2]);
                    double score = Double.parseDouble(fields[3]);
                    if (regionStart >= end) {
                        break;  // Exit early if the current region is beyond the end position
                    }
                    if (regionEnd <= start) {
                        continue;  // Skip the current line if it is before the start position
                    }
                    int effectiveStart = Math.max(regionStart, start);
                    int effectiveEnd = Math.min(regionEnd, end);
                    for (int i = effectiveStart; i < effectiveEnd; i++) {
                        scores[i - start] = score;
                    }
                }
            }
        }
        return scores;
    }

    // Function to fetch PhastCons scores from BedGraph file
    public static double[] fetchPhastConsScores2(String bedGraphFile, String seqname, int start, int end) throws IOException {
	System.out.println( "fetching PhastCons scores ... " );
        double[] scores = new double[end - start];
        int index = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(bedGraphFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                if (fields[0].equals(seqname)) {
                    int regionStart = Integer.parseInt(fields[1]);
                    int regionEnd = Integer.parseInt(fields[2]);
                    double score = Double.parseDouble(fields[3]);
                    for (int i = regionStart; i < regionEnd; i++) {
                        if (i >= start && i < end) {
                            scores[index++] = score;
                        }
                    }
                }
            }
        }
        return scores;
    }

    public static void main(String[] args) {
        try {
            String gtfFile = "data/Homo_sapiens.GRCh38.110.gtf";
            String enstId = "ENST00000212015";  // ENST ID
            String fastaFile = "/usr/local/hg38.fa";
            String outputFile = "data/output.fa";
            String query31mer = "ATAAAACACCCAGCTAGGACCATTACTGCCA";

            // Extract genomic location
            Map<String, String> location = extractGenomicLocation(gtfFile, enstId);
            if (location != null) {
                String seqname = "chr" + location.get("seqname");
		System.out.println( "seqname: " + seqname );
                int start = Integer.parseInt(location.get("start"));
                int end = Integer.parseInt(location.get("end"));
                String strand = location.get("strand");

                // Fetch sequence using bedtools
                fetchSequence(fastaFile, seqname, start, end, outputFile);
                // Search for 31mer and adjust coordinates
                Map<String, Integer> adjustedCoordinates = search31merAndAdjustCoordinates(outputFile, start, query31mer);
		System.out.println( "adjusted coordinates: " + adjustedCoordinates.get( "start" ) + ", " + adjustedCoordinates.get( "end" ) );
                if (adjustedCoordinates != null) {
                    int adjStart = adjustedCoordinates.get("start");
                    int adjEnd = adjustedCoordinates.get("end");

                    // Fetch PhastCons scores
                    String bedGraphFile = "data/hg38.phastCons100way.bedGraph";
                    double[] scores = fetchPhastConsScores(bedGraphFile, seqname, adjStart, adjEnd);

                    // Compute and print PhastCons scores
                    double minScore = Arrays.stream(scores).min().orElse(0.0);
                    double maxScore = Arrays.stream(scores).max().orElse(0.0);
                    double averageScore = Arrays.stream(scores).average().orElse(0.0);
                    double sumScores = Arrays.stream(scores).sum();
                    int lenScores = scores.length;

                    System.out.println("PhastCons scores for " + seqname + ":" + adjStart + "-" + adjEnd + ":");
                    System.out.println("Min: " + minScore + ", Average: " + averageScore + ", Max: " + maxScore + ", Sum: " + sumScores + ", Count: " + lenScores);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
