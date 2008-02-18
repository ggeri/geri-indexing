/**
 * This class represents a kd-tree leaf node. Each leaf node contains
 * the number that maps the leaf nodes to the bins that contain actual
 * points. 
 * @author geri
 *
 */
public class KDLeaf extends Node
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * The mapping to the bins of points
	 */
	private int binIDMapping;		
		
	/**
	 * The default constructor.
	 *
	 */
	public KDLeaf()
	{		
	}
	
	/**
	 * Setter fot the mapping attribute.
	 * @param mapping
	 */
	public void setBinIDMapping(int mapping)
	{
		this.binIDMapping = mapping;
	}
		
	/**
	 * Getter for the mapping attribute.
	 * @return
	 */
	public int getBinIDMapping()
	{
		return this.binIDMapping;
	}

}
