package ca.iric.major.common;
import ca.iric.major.tools.SeedDuplexCache;
import java.util.List;
import java.io.Serializable;

public class EightMerDuplex implements Serializable {

    private static final long serialVersionUID = 1L;

    private int   eightMer = -1;
    private double kdA1 = Guide.SEEDLESSKD;
    private double kdB1 = Guide.SEEDLESSKD;

    public int    get8mer() { return this.eightMer; }
    public double getKdB1() { return this.kdB1; }
    public double getKdA1() { return this.kdA1; }

    // constructor
    // target is a 8mer
    public EightMerDuplex( int target, double kdA1, double kdB1 ) {
	this.eightMer = target; // save the target site
	this.kdA1 = kdA1;
	this.kdB1 = kdB1;
    }

    @Override
    public String toString() {
	return
	    SeedDuplexCache.decode8mer( this.eightMer ) + ", KdB1 = " + this.kdB1 + ", KdA1 = " + this.kdA1;
    }
}
