/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CommonKmers {
    
    public static void main(String[] args) {
	String[] strings = {"ATCGATCGGGA", "GATCGGATAC", "ATCGGATACC"};
        int minK = 3; // Minimum size of the kmer
        int maxK = 5; // Maximum size of the kmer
        int minCount = 1; // Minimum number of strings in which a kmer must appear
        Map<String, List<List<Integer>>> kmerIndices = findKmerIndices(strings, minK, maxK, minCount);
	System.out.println( kmerIndices );
        for (Map.Entry<String, List<List<Integer>>> entry : kmerIndices.entrySet()) {
            System.out.println("Kmer: " + entry.getKey() + ", Positions: " + entry.getValue());
        }
    }

    public static Map<String, List<List<Integer>>> findKmerIndices(String[] strings, int minK, int maxK, int minCount) {
        Map<String, List<List<Integer>>> kmerMap = new HashMap<>();

        // Iterate over each kmer length from minK to maxK
        for (int k = minK; k <= maxK; k++) {
            // Process each string
            for (int s = 0; s < strings.length; s++) {
                String str = strings[s];
                // Initialize a map to track kmers for the current kmer length
                Map<String, List<List<Integer>>> tempMap = new HashMap<>();

                for (int i = 0; i <= str.length() - k; i++) {
                    String kmer = str.substring(i, i + k);
                    if (!tempMap.containsKey(kmer)) {
                        tempMap.put(kmer, new ArrayList<>());
                        // Ensure all kmers have lists initialized for all strings
                        for (int j = 0; j < strings.length; j++) {
                            tempMap.get(kmer).add(new ArrayList<>());
                        }
                    }
                    tempMap.get(kmer).get(s).add(i);
                }

                // Merge temporary map into main kmer map
                tempMap.forEach((key, value) -> {
                    if (!kmerMap.containsKey(key)) {
                        kmerMap.put(key, new ArrayList<>());
                        for (int j = 0; j < strings.length; j++) {
                            kmerMap.get(key).add(new ArrayList<>());
                        }
                    }
                    for (int j = 0; j < value.size(); j++) {
                        kmerMap.get(key).get(j).addAll(value.get(j));
                    }
                });
            }
        }

        // Filter to retain only kmers found in at least minCount strings with proper indexing
        kmerMap.entrySet().removeIf(entry -> {
            int count = (int) entry.getValue().stream().filter(list -> !list.isEmpty()).count();
            return count < minCount;
        });

        return kmerMap;
    }
}
