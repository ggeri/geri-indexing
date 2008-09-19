/**
 * An instance of this class represents an element that we put into the 
 * priority queue that is used for best-bin-first search. PQElement instance has 
 * a nodeID and priority attributes. NodeID is reference to the element we are storing
 * in the queue and priority is the priority of that element.
 * @author geri
 *
 */
public class PQElement implements Comparable
{
	/**
	 * Node id is the reference to the node we are storing in this queue element.
	 */
	private Node nodeID;
	
	/**
	 * The priority of the stored node.
	 */
	private double priority;	

	/**
	 * The constructor that takes the node and its priority and
	 * creates the pq element containing them.
	 * @param nodeID
	 * @param priority
	 */
	public PQElement(Node nodeID, double priority)
	{
		this.nodeID = nodeID;
		this.priority = priority;
	}
	
	/**
	 * Getter for the node id.
	 * @return
	 */
	public Node getNodeID()
	{
		return this.nodeID;
	}
	
	/**
	 * Getter for the priority of the node contained in this pq element.
	 * @return
	 */
	public double getPriority()
	{
		return this.priority;
	}
	
	/**
	 * This method is implemented to make the PQElement objects comparable and
	 * therefore useful in the priority queue data structure.
	 * @param el
	 * @return
	 */
	public int compareTo(Object el)
	{
		if(el instanceof PQElement)
		{
			if(this.priority < ((PQElement) el).priority)
			{
				return -1;
			}else if(this.priority > ((PQElement) el).priority)
			{
				return +1;
			}else
			{
				return 0;
			}	
		}else
		{
			throw new ClassCastException();	
		}	
	}
}
