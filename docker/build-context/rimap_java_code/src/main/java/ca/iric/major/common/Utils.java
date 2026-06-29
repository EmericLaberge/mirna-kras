/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.*;

public class Utils {

    // stop: print a message and exit using error
    public static void stop( String message, int error ) {
	System.out.println( "Error: " + message );
	System.exit( error );
    }

    // debug: print a message and continue
    public static void debug( String message ) {
	System.out.println( message );
    }

    // reverse a String
    public static String reverse( String s ) {
	return new StringBuilder( s ).reverse().toString();
    }

    // count the number of occurrences of a char in a String
    public static int count( String s, char c ) {
	int count = 0;
	for( int i = 0; i < s.length(); i++ )
            if( s.charAt( i ) == c ) count++;
	return count;
    }

    public static boolean isValidDouble( String input ) {
    if( input == null || input.trim().isEmpty() ) return false;
    try {
        Double.parseDouble( input.trim() );
        return true;
    } catch( NumberFormatException e ) {
        return false;
    }
    }

    public static boolean isValidInt( String input ) {
    if( input == null || input.trim().isEmpty() ) return false;
    try {
        Integer.parseInt( input.trim() );
        return true;
    } catch( NumberFormatException e ) {
        return false;
    }
    }

    public static boolean containsOnlyAllowedChars( String input, String allowedChars ) {
	return input != null && input.matches( allowedChars );
    }

    public static boolean isValidAbstractShape( String shape ) {
        if( shape == null ) return false;

        // Remove all whitespace characters
        shape = shape.replaceAll( "\\s", "" );

        if( shape.isEmpty() ) return true;

        Stack<Character> stack = new Stack<>();

        for( char c : shape.toCharArray() ) {
            if( c == '(' ) {
                stack.push( c );
            } else if( c == ')' ) {
                if( stack.isEmpty() ) return false;
                stack.pop();
            } else {
                return false; // Invalid character
            }
        }
        return stack.isEmpty(); // True if all open parens are closed
    }

    public static String convertToNrFormat( String input, String nrValue ) {
	StringBuilder sb = new StringBuilder();
	for( int i = 0; i < input.length(); i++ ) {
	    if( input.charAt( i ) == 'x' ) {
		sb.append( " -nr " ).append( i + 1 ).append( " " ).append( nrValue );
	    }
	}
	return sb.toString().trim(); // remove trailing space
    }
}
