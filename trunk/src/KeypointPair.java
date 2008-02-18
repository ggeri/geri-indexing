import java.io.*;

/**
 * 
 * @author geri
 *
 * An instance of this class represents a pair of keypoints. It extends the Keypoint class.
 * It has 3 attributes: a descriptor array that holds the Keypoint.DESCRIPTOR_LENGHT long
 * keypoint descriptor values - Keypoint.DESCRIPTOR_LENGHT/2 for each point of the pair and
 * two location arrays that contain the 4 locaition values - one
 * array for each of the points of the pair and  
 * 
 */

public class KeypointPair extends Keypoint implements Serializable
{	
	/**
	 * For serialization.
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Array that stores the keypoint location values:
	 * x and y coordinates, scale and orientation 
	 * of the first point in the pair.
	 */
	private double[] locationOne = new double[4];
	
	/**
	 * Array that stores the keypoint location values: 
	 * x and y coordinates, scale and oriantation
	 * of the second point in the pair.
	 */	
	private double[] locationTwo = new double[4];
	
	/**
	 * Default constructor.
	 */
	public KeypointPair()
	{
		
	}
	
	/**
	 * Constructor that calls super with 'descriptor' as parameter and sets the lcoationOne and 
	 * locationTwo parameters
	 * @param locationOne - location of the first element of the pair
	 * @param locationTwo - location of the second element of the pair
	 * @param descriptor - descriptor vector of this pair
	 */
	public KeypointPair(double[] locationOne, double[] locationTwo, short[] descriptor)
	{		
		super(descriptor);
		this.locationOne = locationOne;
		this.locationTwo = locationTwo;		
	}
	
	/**
	 * Setter for the pair's location arrays.
	 * @param loc
	 */
	public void setKeyptLocations(double[] locOne, double[] locTwo)
	{
		this.locationOne[0] = locOne[0];
		this.locationOne[1] = locOne[1];
		this.locationOne[2] = locOne[2];
		this.locationOne[3] = locOne[3];
		
		this.locationTwo[0] = locTwo[0];
		this.locationTwo[1] = locTwo[1];
		this.locationTwo[2] = locTwo[2];
		this.locationTwo[3] = locTwo[3];
	}
	
	/**
	 * Get method for the location array for the first keypooint
	 * of the pair.
	 * @return 
	 */
	public double[] getKeyptLocationOne()
	{
		return this.locationOne;
	}
	
	/**
	 * Get method for the location array for the second keypooint
	 * of the pair.
	 * @return 
	 */
	public double[] getKeyptLocationTwo()
	{
		return this.locationTwo;
	}
}
