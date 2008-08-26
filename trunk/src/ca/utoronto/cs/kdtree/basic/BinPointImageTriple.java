package ca.utoronto.cs.kdtree.basic;
/**
 * An instance of this class represents the triple that consists of the bin
 * number, image number and point in that image. It has three attributes: binID, 
 * pointID and imageID - one for each member of this triple.
 * @author geri
 *
 */
public class BinPointImageTriple {
	
	/**
	 * The id for the bin where the point belongs to.
	 */
	private int binID;
	
	/**
	 * The id of the point.
	 */
	private int pointID;
	
	/**
	 * The id of the image to which the point belongs to.
	 */
	private int imageID;
	
	/**
	 * The constructor that takes binID, pointID and imageID as ints
	 * and constructs an instance of this class with the passed values.
	 * @param binID - the bin where the point belongs to.
	 * @param pointID - the number of the keypoint we are looking at.
	 * @param imageID - the number of the image that contains this point.
	 */
	public BinPointImageTriple(int binID, int pointID, int imageID)
	{
		this.binID = binID;
		this.pointID = pointID;
		this.imageID = imageID;
	}

	/**
	 * The getter method for the binID attribute.
	 * @return - it returns the binID of this point.
	 */
	public int getBinID()
	{
		return this.binID;
	}
	
	/**
	 * The getter method for the pointID of this point.
	 * @return - it returns the pointID of this point.
	 */
	public int getPointID()
	{
		return this.pointID;
	}
	
	/**
	 * The getter for the imageID of the image to which this point belongs to.
	 * @return - it returns the imageID of this point, i.e. the number of the 
	 * image that contains this point.
	 */
	public int getImageID()
	{
		return this.imageID;
	}
}
