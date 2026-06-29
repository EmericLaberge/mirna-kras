/*
 * Copyright (c) 2025 François Major, Major Lab (Université de Montréal)
 * Licensed under the MIT License. See LICENSE file in the project root for details.
 */
package ca.iric.major.common;

import java.util.Scanner;
import java.util.ArrayList;

public class Test {

    public ArrayList<Integer> globalList = new ArrayList<>();
    public ArrayList<Integer> copyOfGlobalList = globalList;
    public int attribut = 4;

    public static void uneMethode( Test a ) {
	System.out.println( a.attribut );
    }

    public void uneMethode() {
	ArrayList<Integer> newList = globalList;
    }

    public static void main( String[] args ) {
	int[] myArray1 = new int[5];
	int[] myArray2 = {1, 2, 3, 4, 5};
	int myArray3[][] = new int[10][10];

	myArray1[5] = 5;
	
	System.out.println( "Bonjour".substring(3) );
	String s = "Java"; s.concat("Langue");
	System.out.println( s );
	s = "";
	System.out.println( s.isEmpty() );
	System.out.println( s.length() == 0 );
	System.out.println( s.equals("") );

	char[] t = { 'a', 'l', 'l', 'o' };
	System.out.println( String.valueOf( t ) );
        //System.out.println( String.toString( t ) );
	String u = new String( t );
	System.out.println( u );
	String v = "Java"; v.toUpperCase(); System.out.println(v);

	uneMethode( new Test() );
    }
}
