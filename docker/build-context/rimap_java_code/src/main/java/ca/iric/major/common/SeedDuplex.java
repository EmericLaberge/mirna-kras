package ca.iric.major.common;
import ca.iric.major.tools.SeedToEightMerCache;
import java.util.List;
import java.io.Serializable;

public class SeedDuplex implements Serializable {

    private static final long serialVersionUID = 1L;

    private int seed = -1;
    private double kdA1 = Guide.SEEDLESSKD;
    private double kdB1 = Guide.SEEDLESSKD;

    public int getSeed()   { return this.seed; }
    public double getKdB1() { return this.kdB1; }
    public double getKdA1() { return this.kdA1; }

    // constructor
    // target is a 8mer; seed is a 7mer
    public SeedDuplex( String target, String seed ) {
	this.seed = SeedToEightMerCache.encode7mer( seed ); // save the design
	// tight seed folding
	String mask = "p".repeat( target.length() ) + Duplex.LOOPMASK + ")" + "q".repeat( seed.length() - 1 );
	Duplex duplex = new Duplex( target, seed, mask, 3, true ); // collect all states
	
	// conformational space and min kd search for B1
	double minKdA1 = Guide.SEEDLESSKD, minKdB1 = Guide.SEEDLESSKD;
	for( Duplex d: duplex.getConformationalSpace() ) { // iterate all states
	    if( !( d.isShifted() || d.isBulgy() ) ) {
		Utils.debug( d.toString() );
		double kdB1 = Guide.determineKd( d, false ); // determineKd, false for not A1
		double kdA1 = Guide.determineKd( d, true );  // determineKd, true for A1
		Utils.debug( "kdB1: " + kdB1 + ", kdA1: " + kdA1 );
		if( kdB1 < minKdB1 ) {
		    minKdB1 = kdB1;
		    this.kdB1 = kdB1;
		}
		if( kdA1 < minKdA1 ) {
		    minKdA1 = kdA1;
		    this.kdA1 = kdA1;
		}
	    }
	}
    }

    @Override
    public String toString() {
	return
	    SeedToEightMerCache.decode7mer( this.seed ) + ", KdB1 = " + this.kdB1 + ", KdA1 = " + this.kdA1;
    }
}
