package ca.utoronto.cs.kdtree.basic;
/**
 * An instance of this class represents the pair that consists of the 'count' 
 * and 'index'. The count gives the number of the points contained in a bin and the index
 * gives the index of the first element of this bin.
 * @author geri
 *
 */
public class CountIndexPair {

	/**
	 * The number of points that belong to this bin.
	 */
	private int count;
	
	/**
	 * The index of the first point that belongs to this bin.
	 */
	private int index;
	
	/**
	 * The default constructor - it sets the count to 0 and the index of the first
	 * point that belongs to this bin to -1 (i.e. invalid index-has to be changed later).
	 *
	 */
	public CountIndexPair()
	{
		this.count = 0;
		this.index = -1;
	}
	
	/**
	 * The getter method for the count attribute.
	 * @return - it returns the value of the count attribute contained in this object.
	 */
	public int getCount()
	{
		return this.count;
	}
	
	/**
	 * The getter method for the index attribute.
	 * @return - it returns the value of the index attribute contained in this object.
	 */
	public int getIndex()
	{
		return this.index;
	}
	
	public void incrementCount()
	{
		this.count++;
	}
	
	public void setIndex(int index)
	{
		this.index = index;
	}
}
