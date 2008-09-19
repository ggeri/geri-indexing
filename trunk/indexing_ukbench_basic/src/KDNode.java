
/**
 * This class represents a kd-tree inner node. Each inner node contains
 * the dimension along which the set is split, the split point in that dimension
 * and the pointers to its left and right child. 
 * @author geri
 *
 */
public class KDNode extends Node
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * The dimension along which the set at this node is split.
	 */
	private int dimension;
	
	/**
	 * The solit point that is used to split this set along the given dimension.
	 */
	private double splitPt;
	
	/**
	 * The pointer to the left child of this node.
	 */
	private int leftChild;
	
	/**
	 * The pointer to the right child of this node.
	 */
	private int rightChild;
	
	/**
	 * This array holds the upper and lower bounds on each dimension of points in this internal node.
	 */
	private double[] dimensionBounds = new double[2];
	
	/**
	 * Lower bound on the dimension.
	 */
	public static final int LO = 0;
	
	/**
	 * Higher bound of the dimension.
	 */
	public static final int HI = 1;
	/**
	 * The constructor. It takes dimension, split point and pointers to the left
	 * and right child as parameters and sets these values in this new node.
	 * @param dim
	 * @param splitPt
	 * @param leftChild
	 * @param rightChild
	 */
	public KDNode(int dim, double splitPt, int leftChild, int rightChild, double[] dimBounds)
	{
		this.dimension = dim;
		this.splitPt = splitPt;
		this.leftChild = leftChild;
		this.rightChild = rightChild;
		this.dimensionBounds[KDNode.LO] = dimBounds[dim];
		this.dimensionBounds[KDNode.HI] = dimBounds[KeypointSingleton.DESCRIPTOR_LENGTH + dim];	
	}
	
	/**
	 * Setter fot the dimension attribute.
	 * @param dim
	 */
	public void setDimension(int dim)
	{
		this.dimension = dim;
	}
	
	/**
	 * Setter for the split point attribute.
	 * @param splitPt
	 */
	public void setSplitPt(double splitPt)
	{
		this.splitPt = splitPt;
	}
	
	/**
	 * Setter for the left child attribute.
	 * @param left
	 */
	public void setLeftChild(int left)
	{
		this.leftChild = left;
	}
	
	/**
	 * Setter for the right child attribute.
	 * @param right
	 */
	public void setRightChild(int right)
	{
		this.rightChild = right;
	}
	
	/**
	 * Setter for the dimensionsBounds attribute that contains lower and 
	 * upper bounds for each dimension in this node.
	 *
	 */
	public void setDimensionsBounds(double[] bounds)
	{
		this.dimensionBounds = bounds;
	}
	
	/**
	 * Getter for the dimension attribute.
	 * @return
	 */
	public int getDimension()
	{
		return this.dimension;
	}
	
	/**
	 * Getter for the split point attribute.
	 * @return
	 */
	public double getSplitPt()
	{
		return this.splitPt;
	}
	
	/**
	 * Getter for the left child attribute.
	 * @return
	 */
	public int getLeftChild()
	{
		return this.leftChild;
	}
	
	/**
	 * Getter for the right child attribute.
	 * @return
	 */
	public int getRightChild()
	{
		return this.rightChild;
	}
	
	/**
	 * Getter for the attribute that represents the array of
	 * the dimension bounds.
	 * @return
	 */
	public double[] getDimensionsBounds()
	{
		return this.dimensionBounds;
	}

}
