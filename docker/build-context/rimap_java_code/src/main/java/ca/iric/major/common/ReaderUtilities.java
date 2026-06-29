/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

public class ReaderUtilities {

    protected static void processException( Throwable e, String fileName ) {
	System.out.println( "Something's wrong while reading " + fileName + "!" );
	e.printStackTrace();
    }
}
