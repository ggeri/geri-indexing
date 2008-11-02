import java.io.*;

/**
 * 
 * @author geri
 *
 * An instance of this class represents a keypoint. It has one 
 * attribute: descriptor array that holds the DESCRIPTOR_LENGTH long keypoint descriptor
 * values.
 * 
 */

public class Keypoint implements Serializable
{	
	/**
	 * For serialization
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * Array that stores the 36 keypoint descriptor values
	 */
	private short[] descriptor = new short[Keypoint.DESCRIPTOR_LENGTH];
	
	/**
	 * The length of the descriptor
	 */
	public static final int DESCRIPTOR_LENGTH = 36; // put back to 36
	
	/**
	 * Default constructor.
	 */
	public Keypoint()
	{
		
	}
	
	/**
	 * A constructor that takes as a parameter descriptor array of short.
	 * @param descriptor
	 */
	public Keypoint(short[] descriptor)
	{
		this.descriptor = descriptor;
	}
	
	/**
	 * Setter for the keypoint descriptor array.
	 * @param descriptor
	 */
	public void setDescriptor(short[] descriptor)
	{
		for(int i = 0; i < Keypoint.DESCRIPTOR_LENGTH; i++)
		{
			this.descriptor[i] = descriptor[i];
		}
	}
	
	/**
	 * Get method for the keypoint descriptor array.
	 * @return
	 */
	public short[] getDescriptor()
	{
		return this.descriptor;
	}
	
}
