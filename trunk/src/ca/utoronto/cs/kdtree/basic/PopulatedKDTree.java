package ca.utoronto.cs.kdtree.basic;
/**
 * The instance of this class represents the populated KD tree. It consists of three arrays: these are the points
 * ordered by the bins to which they belong. Points are represented as elements of the three arrays: first array 
 * gives binID where the point belongs, the second one the imageID where the point is coming from and the third one
 * gives the pointID.
 * Each point is essentially a triple binID-imageID-pointID, where each element comes from one of the three arrays.  
 * @author geri
 *
 */

public class PopulatedKDTree {

	private int[] binIDs;
	private int[] pointIDs;
	private int[] imageIDs;
	
	public PopulatedKDTree(int numImages)
	{
		int avrgPtsPerImage = 5000;	// my approximation
		int arraySize = avrgPtsPerImage * numImages;
		binIDs = new int[arraySize];
		pointIDs = new int[arraySize];
		imageIDs = new int[arraySize];
	}
	
	
	public PopulatedKDTree(int[] binIDs, int[] pointIDs, int[] imageIDs)
	{
		this.binIDs = binIDs;
		this.pointIDs = pointIDs;
		this.imageIDs = imageIDs;
	}
	
	public int[] getBinIDs()
	{
		return this.binIDs;
	}
	
	public int[] getPointIDs()
	{
		return this.pointIDs;
	}
	
	public int[] getImageIDs()
	{
		return this.imageIDs;
	}
	
	public void setBinIDs(int[] binIDs)
	{
		this.binIDs = binIDs;
	}
	
	public void setPointIDs(int[] pointIDs)
	{
		this.pointIDs = pointIDs;
	}
	
	public void setImageIDs(int[] imageIDs)
	{
		this.imageIDs = pointIDs;
	}
}
