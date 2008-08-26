package ca.utoronto.cs.kdtree.basic;
import java.io.*;

/**
 * 
 * @author geri
 *
 * An instance of this class represents a keypoint singleton. It extends Keypoint class. It has two 
 * attributes: descriptor array that holds the Keypoint.DESCRIPTOR_LENGTH long 
 * keypoint descriptor values (this is inhereted from Keypoint parent class and 
 * location array that contain the 4 locaition values.
 * 
 */

public class KeypointSingleton extends Keypoint implements Serializable
{	
	/**
	 * For serialization.
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Array that stores the keypoint location values:
	 * x and y coordinates, scale and orientation.
	 */
	private double[] location = new double[4];

	/**
	 * Default constructor.
	 */
	public KeypointSingleton()
	{
		
	}
	
	/**
	 * Constructor that calls super with 'descriptor' as parameter and sets the lcoation parameter
	 * @param location - location of the keypoint
	 * @param descriptor - descriptor verctor
	 */
	public KeypointSingleton(double[] location, short[] descriptor)
	{		
		super(descriptor);
		this.location = location;
	}
	
	/**
	 * Setter for the keypoint location array.
	 * @param loc
	 */	
	public void setKeyptLocation(double[] loc)
	{
		this.location[0] = loc[0];
		this.location[1] = loc[1];
		this.location[2] = loc[2];
		this.location[3] = loc[3];
	}
	
	/**
	 * Get method for the keypoint location array.
	 * @return 
	 */
	public double[] getKeyptLocation()
	{
		return location;
	}

	
}
