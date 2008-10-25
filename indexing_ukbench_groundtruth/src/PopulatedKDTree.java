import java.util.Vector;

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
	
	public void setBinID(int position, int binIDValue)
	{
		this.binIDs[position] = binIDValue;
	}
	
	public void setPointID(int position, int pointIDValue)
	{
		this.pointIDs[position] = pointIDValue;
	}
	
	public void setImageID(int position, int imageIDValue)
	{
		this.imageIDs[position] = imageIDValue;
	}
	
	public void sortPointsByBinID(Vector<CountIndexPair> binCounts) 
	{
		// place one by one point in its place in populatedTree arrays				
		int totalPoints = this.pointIDs.length;
		
		int[] populatedBinIDs = new int[totalPoints];
		int[] populatedPointIDs = new int[totalPoints];
		int[] populatedImageIDs = new int[totalPoints];
		
		// each cell gives where exactly we got with filling this bin
		int[] counter = new int[binCounts.size()];	//index in counter corresponds to binID
		
		for(int i = 0; i < totalPoints; i++)
		{								
			int bin = this.binIDs[i];	
			
			CountIndexPair pair = binCounts.get(bin);
			int index = pair.getIndex();	// this is the index of first element of this bin ('bins' vector)
			
			// the index where the next element goes is index+counter[bin]
			int nextElSlot = index + counter[bin];
			
			counter[bin]++;
			
			populatedBinIDs[nextElSlot] = this.binIDs[i];											
			populatedPointIDs[nextElSlot] = this.pointIDs[i];
			populatedImageIDs[nextElSlot] = this.imageIDs[i];
		
		}
	}
}
