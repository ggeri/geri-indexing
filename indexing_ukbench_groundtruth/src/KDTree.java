import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;


/**
 * This class creates and train a kd-tree on the passed set
 * of points.
 * @author geri
 *
 */
public class KDTree
{
	/**
	 * The threshold for the size of a kd-tree node.
	 * If the size is less or equal this treshold, the node will
	 * not be split further. 
	 */	
	public static final int THRESHOLD = 10;	
	
	
	/**
	 * This method is used for creating/training the kd-tree.
	 * It takes a set of points as a paramters and it trains the tree on these points.
	 * It returns the reference to the root of the built tree.
	 * @param args
	 */

// Semantic change: leftInd and rightInd refer to the index of a point. To access the 3rd element of the first 
// point in a range, use pointsSet[leftInd*Keypoint.DESCRIPTOR_LENGTH + 2]
// leftInd is the index of the first point in the range
// rightInd is the index of the last point in the range
// so if the range is {0,1,2}, leftInd = 0 and rightInd = 2
// on initialization, leftInd should be 0 and rightInd should be  (size(pointsSet)/Keypoint.DESCRIPTOR_LENGTH)-1	 
	public int trainKDTree(short[] pointsSet, int leftInd, int rightInd, Vector<Node> trainedTree, double[] dimBounds)	
			throws Exception
	{						
		// rightInd - leftInd + 1  = number of points in the current pointsSet array
		if((rightInd-leftInd+1) > KDTree.THRESHOLD)  
		{	
			// find the dimension with the greatest range
			int dim = this.pickDimension(pointsSet, leftInd, rightInd);
			if (dim < 0) // this means that all points in the range are identical: very unlikely
			{
				// just make ourself into a leaf node
				int index = trainedTree.size();
				trainedTree.add(trainedTree.size(), new KDLeaf());	
				return index;
			}
			else if (dim >= Keypoint.DESCRIPTOR_LENGTH)
			{
				throw new Exception("Dimension picked is outside expected bounds (bounds of the descriptor length");
			}
		
			// find the mean value along this dimension
			double splitPt = this.pickSplitPoint(pointsSet, leftInd, rightInd, dim);			

			// find the sliding midpoint along this dimension
			//double splitPt = this.pickSplitPointSliding(pointsSet, leftInd, rightInd, dim, dimBounds);

			// find the highest value <= splitPt
			int splitPtInd = this.getSplitPtIndex(pointsSet, leftInd, rightInd, dim, splitPt);
			if (splitPtInd < leftInd || splitPtInd > rightInd)
			{
				throw new Exception("splitPtInd is outside expected bounds");
			}
			
			// sort the values such that [left,newSplitPtInd] <= splitPt and [newSplitPtInd+1, right] > splitPt
			int newSplitPtInd = this.inPlacePartition(pointsSet, leftInd, rightInd, dim, splitPtInd);
			if (newSplitPtInd < leftInd || newSplitPtInd >= rightInd)
			{
				throw new Exception("newSpitPtInd is outside expected bounds");
			}
			
			// add ourself to the tree as an internal node
			int index = trainedTree.size();
			trainedTree.add(trainedTree.size(), new KDNode(dim, splitPt, -1, -1, dimBounds));
			
			// change the bounds and pass to children
			double[] dimBoundsLeft = new double[2*Keypoint.DESCRIPTOR_LENGTH];			
			System.arraycopy(dimBounds, 0, dimBoundsLeft, 0, dimBounds.length);
			double[] dimBoundsRight = new double[2*Keypoint.DESCRIPTOR_LENGTH]; 
			System.arraycopy(dimBounds, 0, dimBoundsRight, 0, dimBounds.length);						
			dimBoundsLeft[dim + Keypoint.DESCRIPTOR_LENGTH] = splitPt;					
			dimBoundsRight[dim] = splitPt;								

			// further subdivide our left and right children
			int left = trainKDTree(pointsSet,leftInd, newSplitPtInd, trainedTree, dimBoundsLeft);
			int right = trainKDTree(pointsSet, newSplitPtInd+1, rightInd, trainedTree, dimBoundsRight);
			
			// set up the references for our left and right children
			KDNode node = (KDNode) trainedTree.elementAt(index);
			node.setLeftChild(left);
			node.setRightChild(right);
			return index;			
			
		}else
		{	
			// just make ourself into a leaf node
			int index = trainedTree.size();
			trainedTree.add(trainedTree.size(), new KDLeaf());	
			return index;
		}			
	}
		
	/**
	 * This method takes an array, the left and right index of the subarray we want partitioned
	 * and the index of the split point (that is the pivotInd). It partitions the subarray into
	 * left and right part, where the left part has all the elements <= splitPt and
	 * the right one has all the elements > splitPt. This partitioning is done in
	 * place, there are no new allocations of memory. The method returns the index of the split 
	 * after partitioning. The index of the left and right end of the subarray stay the same.
	 * so, [left, doneUpToHere] is all <= the split value, and [doneUpToHere+1, right] is all > the split value
	 * @param dataSet
	 * @param leftInd - index of the first descriptor vector in this subarray
	 * @param rightInd - index of the last descriptor vector in this subarray
	 * @param dim - dimension on which we will split
	 * @param splitPtInd- index of the descriptor vector on which to split
	 * @return
	 */
	public int inPlacePartition(short[] dataSet, int leftInd, int rightInd, int dim, int splitPtInd)
	{	
		//pivotInd is the index of the splitPt
		int pivotInd = splitPtInd;
		int pivotValue = dataSet[pivotInd*Keypoint.DESCRIPTOR_LENGTH + dim];
		
		// we put the 36-vector that contains pivot value off to the side - all the way to the right
		this.swap(dataSet, pivotInd, rightInd);
		
		// we start at left index and move towards right
		int doneUpToHere = leftInd;
		for(int i = leftInd; i < rightInd; i++)
		{
			if(dataSet[i*Keypoint.DESCRIPTOR_LENGTH + dim] <= pivotValue)
			{								
				this.swap(dataSet, doneUpToHere, i);
				doneUpToHere++;
			}
		}
		
		// put pivot (splitPt) back in its place		
		this.swap(dataSet, rightInd, doneUpToHere);
		
		// done
		return doneUpToHere;
	}
	
	private void swap(short[] dataSet, int firstInd, int secondInd)
	{
		short tempValue;
		for(int i = 0; i < Keypoint.DESCRIPTOR_LENGTH; i++)
		{
			tempValue = dataSet[secondInd*Keypoint.DESCRIPTOR_LENGTH + i];
			dataSet[secondInd*Keypoint.DESCRIPTOR_LENGTH + i] = dataSet[firstInd*Keypoint.DESCRIPTOR_LENGTH + i];
			dataSet[firstInd*Keypoint.DESCRIPTOR_LENGTH + i] = tempValue;
		}
	}
	
	
    /**
     * Chose the dimension along which to split the given set of points
     * We chose the dimension with the largest difference between the largest 
	 * and smallest coordinate for this dimension.
	*/
	private int pickDimension(short[] pointsSet, int leftInd, int rightInd)
	{
		int maxDiff = 0;
		int maxDim = -1;
		
		// loop over dimensions and find the max and the min and their diff - then pick largest difference
		for(int currDim = 0; currDim < Keypoint.DESCRIPTOR_LENGTH; currDim++)
		{
			int dimMin = pointsSet[leftInd*Keypoint.DESCRIPTOR_LENGTH + currDim];	//initialize it to the first element of this dimension
			int dimMax = pointsSet[leftInd*Keypoint.DESCRIPTOR_LENGTH + currDim];	//initialize it to the first element of this dimension			
			
			for(int currInd = leftInd+1; currInd <= rightInd; currInd++)
			{				
				short currValue = pointsSet[currInd*Keypoint.DESCRIPTOR_LENGTH + currDim];
				
				// loop over all points (feature vectors of all points) and find min
				if(currValue < dimMin)
				{
					dimMin = currValue;
				}
				
				// loop over all points and find max						
				if(currValue > dimMax)
				{
					dimMax = currValue;
				}								 								
			}
			
			// get the difference
			int diff = dimMax - dimMin;
			
			//System.out.println("dimMax, dimMin and diff: " + dimMax + " " + dimMin + " " + diff);
			
			// check against biggest difference so far
			if (diff > maxDiff)
			{
				maxDiff = diff;
				maxDim = currDim;
			}
		}
		
		if(maxDiff == 0)
		{
			System.out.println("maxDiff is : " + maxDiff + " !!!!!!!!!!!!! ");
		}
				
		return maxDim;
	}
	
	
	/**
	 * This method is used to pick the split point along the chosen dimension.
	 * It takes as arguments the set of points and the dimension to split along
	 * and it returns the split point. We choose for the split point the mean of
	 * the points along this dimension.
	 * @param pointsSet
	 * @param dim
	 * @return
	 */
	private double pickSplitPoint(short[] pointsSet, int leftInd, int rightInd, int dim)
	{
		double sum = 0;				
		for(int i = leftInd; i <= rightInd; i++)
		{			
			sum = sum + pointsSet[i*Keypoint.DESCRIPTOR_LENGTH + dim];
		}
		double mean = sum /(rightInd - leftInd + 1);
		return mean;
	}
	
	/**
	 * This method is used to pick the split point along the chosen dimension.
	 * It takes as arguments the set of points and the dimension to split along
	 * and it returns the split point. The split is in the center, unless that
	 * would leave all the points to one side. In that case, the split is shifted
	 * so that there is at least one point on each side.
	 * @param pointsSet
	 * @param dim
	 * @return
	 */
	private double pickSplitPointSliding(short[] pointsSet, int leftInd, int rightInd, int dim, double[] dimBounds)
	{
		
		// find the min and max values in this dimension	
		double currVal, minVal, maxVal;
		minVal = pointsSet[leftInd*Keypoint.DESCRIPTOR_LENGTH + dim];
		maxVal = minVal;				
		for(int i = leftInd+1; i <= rightInd; i++)
		{
			currVal = pointsSet[i*Keypoint.DESCRIPTOR_LENGTH + dim];
			if (currVal < minVal)
			{
				minVal = currVal+0.5;
			}
			else if (currVal > maxVal)
			{
				maxVal = currVal-0.5;
			}
		}

		// find the midpoint
		double midpt = (dimBounds[dim] + dimBounds[dim + Keypoint.DESCRIPTOR_LENGTH])/2;

		// shift if needed
		if (midpt <= minVal)
		{
			midpt = minVal + 0.5;
		}
		else if (midpt >= maxVal)
		{
			midpt = maxVal - 0.5;
		}
		
		// all done
		return midpt;
	}

	/**
	 * Given the pointsSet, the dimension and the splitPoint, this method returns the index of 
	 * a point in the pointsSet along the given dimension that is the closest value less than or equal to splitPt
	 */
	public int getSplitPtIndex(short[] dataSet, int leftInd, int rightInd, int dim, double splitPt)
	{		
		int bestInd = -1;
		double bestDiff = Double.MAX_VALUE;
		double currVal;
		double currDiff;
		for(int i = leftInd; i <= rightInd; i++)
		{
			currVal = dataSet[i*Keypoint.DESCRIPTOR_LENGTH + dim];
			currDiff = splitPt - currVal;
			if (currDiff >= 0 && currDiff < bestDiff)
			{
				bestDiff = currDiff;
				bestInd = i;
			}
		}
		
		try {
		if((bestInd < leftInd) || (bestInd > rightInd))
			throw new Exception("Could not find point <= splitPt in range");
		} catch (Exception e) {				
			e.printStackTrace();
		}
			
		return bestInd;
	}
	
	/**
	 * Assigns the mapping form the leaf nodes in the tree to the 
	 * bins that contain the points. The passed vector is the kd-tree
	 * we built and we are assigning mapping to its leaf nodes. This method
	 * returns the number of leaf nodes.
	 * @param tree
	 */
	public int mapLeafsToBins(Vector tree)
	{
		int m = 0;		
		for(int i = 0; i < tree.size(); i++)
		{
			if(tree.get(i) instanceof KDLeaf)
			{
				KDLeaf leaf = (KDLeaf)tree.get(i);
				leaf.setBinIDMapping(m);
				m++;			
			}
		}
		return m;
	}
	
	
	/**
	 * This method takes a query point and returns the a leaf node 
	 * that contain points close to it i.e. in the same bin with it.
	 * @param tree
	 * @param dataSetPopulating
	 * @param keyptStart
	 * @return
	 */
	public KDLeaf nearestBin(Vector tree, short[] dataSetPopulating, int keyptStart)
	{						
		int index = 0;
		while(tree.get(index) instanceof KDNode)
		{
			KDNode innerNode = (KDNode)tree.get(index);
			int splitDim = innerNode.getDimension();
			double splitPt = innerNode.getSplitPt();
			
			if(dataSetPopulating[keyptStart + splitDim] < splitPt)		// left of splitting place
			{
				index = innerNode.getLeftChild();
			}else
			{
				index = innerNode.getRightChild();
			}
		}
		KDLeaf leafNeighbour = (KDLeaf)tree.get(index);		
		
		return leafNeighbour;		
	}	
	
	/**
	 * This method takes a query point and returns the a leaf node 
	 * that contain points close to it i.e. in the same bin with it.
	 * @param tree
	 * @param dataSetPopulating
	 * @param keyptStart
	 * @return
	 */
	public KDLeaf nearestBin(Vector trainedTree, short[] pointData)
	{						
		int index = 0;
		while(trainedTree.get(index) instanceof KDNode)
		{
			KDNode innerNode = (KDNode)trainedTree.get(index);
			int splitDim = innerNode.getDimension();
			double splitPt = innerNode.getSplitPt();
			
			if(pointData[splitDim] < splitPt)		// left of splitting place
			{
				index = innerNode.getLeftChild();
			}else
			{
				index = innerNode.getRightChild();
			}
		}
		KDLeaf leafNeighbour = (KDLeaf)trainedTree.get(index);		
		
		return leafNeighbour;		
	}
	
	
	
	/**
	 * This method takes as a paramenter the tree, the query point and 2 tresholds.
	 * It returns the bins up to maximum 'maxBins' bins that contain the points at most 'threshold' away 
	 * from the query point.
	 * @param tree
	 * @param queryPoint
	 * @param threshold
	 * @param maxBins
	 * @return 
	 */
	public Vector<Integer> bestBinFirst(Vector tree, double[] distances, Keypoint query, double threshold, int maxBins)
	{
		short[] queryPt = query.getDescriptor();
		Vector<Integer> binIDs = new Vector<Integer>(maxBins);	
		
		PriorityQueue<PQElement> pq = new PriorityQueue<PQElement>();
		Node root = (Node)tree.get(0); 
		pq.add(new PQElement(root, 0));
		int i = 0;
		while(binIDs.size() < maxBins && !pq.isEmpty())		//size() starts at 0 so we need < maxBins
		{
			PQElement el = (PQElement)pq.poll();
			Node nodeID = el.getNodeID();
			double distance = el.getPriority();
			int binId = this.bestBinFirstNode(tree, queryPt, pq, nodeID, distance, threshold);
			binIDs.add(new Integer(binId));
			distances[i] = distance;
			i++;
		}
		// trim vector binIDs to its actual size
		binIDs.trimToSize();
		return binIDs;
	}
	
	
	
	/**
	 * This method is called from the bestBinFirst method. If the current node is a leaf node, it returns
	 * the bin id of this leaf and if the current node is not a leaf node, it is called recursivelly on it.
	 * It takes the following parameters: 
	 * @param tree - the kd-tree that we are searching on
	 * @param queryPt - the query point
	 * @param pq - priority queue which is used for priority search
	 * @param node - the current node
	 * @param oldDist - distance to the current node
	 * @param threshold - threshold on the distance - only the points that are less then 'threshold' away are accounted for
	 * @return - the bin id of the leaf node found
	 */
	private int bestBinFirstNode(Vector tree, short[] queryPt, PriorityQueue<PQElement> pq, 
			Node node, double oldDist, double threshold)	//dist is already squared
	{
		if(node instanceof KDLeaf)
		{
			KDLeaf leaf = (KDLeaf)node;					
			return leaf.getBinIDMapping();
		}else
		{
			KDNode kdnode = (KDNode)node;
			int l = kdnode.getLeftChild();
			int r = kdnode.getRightChild();
			Node leftChild = (Node)tree.get(l);
			Node rightChild = (Node)tree.get(r);
			
			int splitDim = kdnode.getDimension();
			double splitPt = kdnode.getSplitPt();
			
			double newDist;
			double dist;
			double takeOutDist;
			
			if(queryPt[splitDim] < splitPt)		//query point is left of split pt - left child is closer
			{				
				if(queryPt[splitDim] < kdnode.getDimensionsBounds()[KDNode.LO])	//it is outside of the parents bounds (left of lower bound)
				{
					dist = (splitPt - queryPt[splitDim]) * (splitPt - queryPt[splitDim]);				
					takeOutDist = (queryPt[splitDim] - kdnode.getDimensionsBounds()[KDNode.LO]) *
									(queryPt[splitDim] - kdnode.getDimensionsBounds()[KDNode.LO]);
					//substract already accounted distance
					newDist = dist + oldDist - takeOutDist;	
				}else
				{
					dist = (splitPt - queryPt[splitDim]) * (splitPt - queryPt[splitDim]);
					newDist = dist + oldDist;	
				}
				if(newDist <= threshold)
				{
					pq.add(new PQElement(rightChild, newDist));
				}				
				return bestBinFirstNode(tree, queryPt, pq, leftChild, oldDist, threshold);
				
			}else	//query point is to the right of split pt - right child is closer
			{			
				if(queryPt[splitDim] > kdnode.getDimensionsBounds()[KDNode.HI])	//it is outside of the paresnt bounds (right of higher bound)
				{
					dist = (queryPt[splitDim] - splitPt) * (queryPt[splitDim] - splitPt);
					takeOutDist = (kdnode.getDimensionsBounds()[KDNode.HI] - queryPt[splitDim]) * 
					                (kdnode.getDimensionsBounds()[KDNode.HI] - queryPt[splitDim]);
					//substract already accounted distance
					newDist = dist + oldDist - takeOutDist;		
				}else
				{
					dist = (queryPt[splitDim] - splitPt) * (queryPt[splitDim] - splitPt);
					newDist = dist + oldDist;				
				}
				if(newDist <= threshold)
				{
					pq.add(new PQElement(leftChild, newDist));
				}
				return bestBinFirstNode(tree, queryPt, pq, rightChild, oldDist, threshold);
			}					
		}
	}
	
	
	
	/**
	 * A new new version of the populateKDTree method which avoids numerous allocations
	 * of memory. It also avoids reading images in a loop and getting the Keypoint[] keyptArr in that loop.
	 * 
	 * The change: we represent the populated tree with three arrays of ints. The size of all three arrays is
	 * the same - it is equal: average_num_points_per_image * num_images = 5000 * numImages
	 * First array holds all the binIDs, second one holds all the pointIDs and the third one 
	 * holds all the imageIDs. Combination of three cells (one from each array) with the same index
	 * will represent one element in the populated tree.	
	 * 
	 * Another change: Instead of reading images one by one in a loop and getting the Keypoint[] array, it instead 
	 * has an array of all the keypoints that it will be using passed to it. This array dataSetPopulating is just a big
	 * array of short. The imageIndeces array is also passed to the method - it gives the indeces of 1st element
	 * of each image in this big dataSetPopulating array.	 
	 */
	
	
	public PopulatedKDTree populateKDTree(Vector trainedTree, short[] dataSetPopulating, int[] imageIdFromMap, int[] imageIndeces,
			Vector<CountIndexPair> counts)
	{
		
		// populatedTree object holds three arrays - these are the points ordered by the bins to which they belong -
		// points are represented as elements of the three arrays - first array gives the binID where the point belongs,
		// the second one the imageID where the point is coming from and the third one the pointID
		// (each point is as triples binID-imageID-pointID , where each element comes from one of the three arrays)
		// the number of points in each bin and the index of first point of each bin are given in the 'counts' vector 'bin' 
		
		PopulatedKDTree populatedTree;
		
		// Temp arrays of all the points - same as populated tree but with the different order
		int arraySize = dataSetPopulating.length / Keypoint.DESCRIPTOR_LENGTH;
		int[] binIDs = new int[arraySize];			
		int[] pointIDs = new int[arraySize];	
		int[] imageIDs = new int[arraySize];
		
		// populate 'counts' vector in advance
		int numBins = this.mapLeafsToBins(trainedTree);
		// System.out.println("Number of bins (leaf nodes) in this tree " + numBins + "\n");
		for(int k = 0; k < numBins; k++)
		{
			CountIndexPair pair = new CountIndexPair();
			counts.add(pair);
		}
		// trim the counts vector to its actual size
		counts.trimToSize();	
		
		int nextElIndx = 0;
		// loop over dataSetPopulating array and process each keypoint (keypt size is 36 - step is 36)
		// we loop over dataSetPopulating array image by image - we can do this by following indeces in imageIndeces
		for(int im = 0; im < imageIndeces.length-1; im++)	
		{												
			// imageID
			int imageID = imageIdFromMap[im]; 
			
            // loop over all keypoints in dataSetPopulating that belong to this curent image
			int pt = 0;
			for(int i = imageIndeces[im]; i < imageIndeces[im+1]; i = i + Keypoint.DESCRIPTOR_LENGTH)
			{
				// pointID - which pt this is (inside of this image)
				int pointID = pt;
				pt++;
				
				// binID (pass to nearestBin the subarray of dataSetPopulating)
				int keyptStart = imageIndeces[im] + (pointID * Keypoint.DESCRIPTOR_LENGTH);					
				KDLeaf leaf = this.nearestBin(trainedTree, dataSetPopulating, keyptStart);
				int binID = leaf.getBinIDMapping();	
				
                // now we put the thee values in the arrays					
				binIDs[nextElIndx] = binID;
				
				pointIDs[nextElIndx] = pointID;				
				
				imageIDs[nextElIndx] = imageID;				
				
				nextElIndx++;					
									
				// get the count-index pair in this bin and increment its count
				CountIndexPair pair = counts.get(binID);
				pair.incrementCount();
			}
		}
		
        // done iterating through all images and all keypoints
		// iterate through the 'counts' vector and set the indeces of first elements for each bin
		int firstElIndx = 0;
		for(int i = 0; i < counts.size(); i++)
		{
			CountIndexPair pair = counts.get(i);
			pair.setIndex(firstElIndx);
			firstElIndx = firstElIndx + pair.getCount();			
		}
		
		// place one by one point in its place in populatedTree arrays				
		int[] populatedBinIDs = new int[nextElIndx];
		int[] populatedPointIDs = new int[nextElIndx];
		int[] populatedImageIDs = new int[nextElIndx];
		
		// each cell gives where exactly we got with filling this bin
		int[] counter = new int[counts.size()];	//index in counter corresponds to binID
		
		//test if binIDs.length = length i gave it or the size of the filled array
		//and accordingly change the loop below
		
		for(int i = 0; i < nextElIndx; i++)
		{								
			int bin = binIDs[i];	
			
			CountIndexPair pair = counts.get(bin);
			int index = pair.getIndex();	// this is the index of first element of this bin ('bins' vector)
			
			// the index where the next element goes is index+counter[bin]
			int nextElSlot = index + counter[bin];
			
			counter[bin]++;
			populatedBinIDs[nextElSlot] = binIDs[i];											
			populatedPointIDs[nextElSlot] = pointIDs[i];
			populatedImageIDs[nextElSlot] = imageIDs[i];
			
		}
		
		populatedTree = new PopulatedKDTree(populatedBinIDs, populatedPointIDs, populatedImageIDs);
		
	return populatedTree;
}

	
	
	/**
	 * A new new version of the populateKDTree method which avoids gathering all the points for populating at once.
	 * This method is called per image from the IndexingPairsTest1 to populate the tree one image at a time.
	 * 
	 * We still represent the populated tree with three arrays of ints. The size of all three arrays is
	 * the same - it is set to the number of points in the database - we get this number by using matlab to count them.
	 * First array holds all the binIDs, second one holds all the pointIDs and the third one 
	 * holds all the imageIDs. Combination of three cells (one from each array) with the same index
	 * will represent one element in the populated tree.	
	 * 
	 * The imageIndeces array is also passed to the method - it gives the indeces of 1st element
	 * of each image in this big dataSetPopulating array.	
	 * trainedTree - trained tree that is to be populated
	 * populatedTree - the tree to populate with points from this image
	 * imageDataSet - points from this image - they are to be put into the tree
	 * imageID - id of the image that we are currently populating with
	 * currPos - where did we stop with populating in this tree
	 * binCounts - this array consists of CountIndexPair elements - it contains as many of those elements as there 
	 * is bins in the tree. Each CountIndexPair (one per bin) gives count - number of points in this bin, and index - 
	 * index of the first point of this bin - it is used later for voting 
	 */
	
	public int populateKDTreePerImage(Vector trainedTree, PopulatedKDTree populatedTree, short[] imageDataSet, int imageID, int currPos,
			Vector<CountIndexPair> binCounts)
	{
		
		// populatedTree object holds three arrays - these are the points ordered by the bins to which they belong -
		// points are represented as elements of the three arrays - first array gives the binID where the point belongs,
		// the second one the imageID where the point is coming from and the third one the pointID
		// (each point is as triples binID-imageID-pointID , where each element comes from one of the three arrays)
			
		// the number of points in each bin and the index of first point of each bin are given in the 'counts' vector 'bin' 
		
		// total points in this image
		int totalPts = imageDataSet.length / Keypoint.DESCRIPTOR_LENGTH;
		
		for(int pointID = 0; pointID < totalPts; pointID++)
		{
			KDLeaf leaf = this.nearestBin(trainedTree, imageDataSet);
			int binID = leaf.getBinIDMapping();
			
			// put the three values in the populatedTree arrays
			populatedTree.setImageID(currPos, imageID);
			populatedTree.setPointID(currPos, pointID);
			populatedTree.setBinID(currPos, binID);
			currPos++;

			// get the count-index pair in this bin and increment its count
			CountIndexPair pair = binCounts.get(binID);
			pair.incrementCount();
			
		}
		
		// Done iterating through all images and all keypoints - 
		// iterate through the 'counts' vector and set the indeces of first elements for each bin
		int firstElIndx = 0;
		for(int i = 0; i < binCounts.size(); i++)
		{
			CountIndexPair pair = binCounts.get(i);
			pair.setIndex(firstElIndx);
			firstElIndx = firstElIndx + pair.getCount();			
		}
		
		return currPos;
}
	
	/**
	 * This method takes a keypoint and queries the trained tree with it.
	 * It returns the binID of the bin to which it belongs.
	 * @param trainedTree - 
	 * @param queryPoint
	 * @return
	 */
	public int query(Vector trainedTree, Keypoint queryPoint)
	{
		
		// 0 passed to nearestBin here does not mean anything - nearestBin had to have that parameter
		// for usage in  populateKDTree
		KDLeaf leaf = this.nearestBin(trainedTree, queryPoint.getDescriptor(), 0);
		return leaf.getBinIDMapping();
	}
	
	
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the query
	 * point and the number of images.
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public void voteForPointOneVotePerImage(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Keypoint queryPoint, int queryPtID, Vector<CountIndexPair> counts, double[] votes, int[] lastQueryPtVoted)
	{				
		// threshold for max distance between this query pt and points that we consider close enough
		double threshold = 100000000000000000000000000.0;
		// threshold on number of bins we want to look into to find reasonably close points to this query point
		int maxBins = 25;
		
		// get all the bins tha contain points 'reasonably' close to the query point		
		Vector<Integer> binIDs = this.bestBinFirst(trainedTree, new double[maxBins], queryPoint, threshold, maxBins);
		
		// iterate over all returned bins (the bins with 'close' points)				
		for(int i = 0; i < binIDs.size(); i++)
		{
			// for each bin, iterate over its points and vote for 
			// the images where those points belong to
			
			// get a bin
			int binID = binIDs.get(i);
			// get its count and index of its first element
			int count = counts.get(binID).getCount();
			int firstElIndex = counts.get(binID).getIndex();
			
			// go into the populated tree, go to the point at index firstElIndex and 
			// iterate over next count points										
			for(int j = firstElIndex; j < firstElIndex + count; j++)
			{
				// get the imageID of the retrieved point				
				int imageID = (populatedTree.getImageIDs())[j];
				
				// vote for the image imageID				
				if(lastQueryPtVoted[imageID] != queryPtID)
				{
					votes[imageID]++;
					lastQueryPtVoted[imageID] = queryPtID;					
				}						
			}			
		}	
	}
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the query
	 * point and the number of images.
	 * It does not care if one point votes for multiple points in an image
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public void voteForPoint(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Keypoint queryPoint, int queryPtID, Vector<CountIndexPair> counts, double[] votes)
	{				
		// threshold for max distance between this query pt and points that we consider close enough
		double threshold = 0.7396;
		// threshold on number of bins we want to look into to find reasonably close points to this query point
		int maxBins = 25;
		
		// get all the bins tha contain points 'reasonably' close to the query point		
		Vector<Integer> binIDs = this.bestBinFirst(trainedTree, new double[maxBins], queryPoint, threshold, maxBins);
		
		// iterate over all returned bins (the bins with 'close' points)				
		for(int i = 0; i < binIDs.size(); i++)
		{
			// for each bin, iterate over its points and vote for 
			// the images where those points belong to
			
			// get a bin
			int binID = binIDs.get(i);
			// get its count and index of its first element
			int count = counts.get(binID).getCount();
			int firstElIndex = counts.get(binID).getIndex();
			
			// go into the populated tree, go to the point at index firstElIndex and 
			// iterate over next count points										
			for(int j = firstElIndex; j < firstElIndex + count; j++)
			{
				// get the imageID of the retrieved point				
				int imageID = (populatedTree.getImageIDs())[j];
				
				// vote for the image imageID				
				//if(lastQueryPtVoted[imageID] != queryPtID)
				//{
					votes[imageID]++;
				//	lastQueryPtVoted[imageID] = queryPtID;					
				//}						
			}			
		}	
	}
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the query
	 * point and the number of images.
	 * It does not care if one point votes for multiple points in an image
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 * It also opens the file with name outputFileName and adds to it the info about retrieved bins, imgs, pt indeces
	 */	
	public void voteForPointWithGroundTruth(Vector<Node> trainedTree, PopulatedKDTree populatedTree, int queryImageID,
			Keypoint queryPoint, int queryPtID, Vector<CountIndexPair> counts, double[] votes, BufferedWriter writer)
	{				
		// threshold for max distance between this query pt and points that we consider close enough
		// this threshold is squared since that is how we use it in KD-tree
		// threshold=0.86 so (0.86*max-short)^2
		// The threshold for KD-tree distance also has to be scaled to 
		// short range as the features descriptor values.
		double threshold = (0.4 * Short.MAX_VALUE) * (0.4 * Short.MAX_VALUE);
		// threshold on number of bins we want to look into to find reasonably close points to this query point
		int maxBins = IndexingGroundTruthTest2.NUM_LEAFS_IN_TRAINED_TREE; // we don't want to limit the number of bins
		
		// get all the bins tha contain points 'reasonably' close to the query point	
		double[] distances = new double[maxBins];
		Vector<Integer> binIDs = this.bestBinFirst(trainedTree, distances, queryPoint, threshold, maxBins);
		
		// Calculate how many points in total (accross all bins) did we retrieve
		int totalPts = 0;
		for(int i = 0; i < binIDs.size(); i++) 
		{
			int binID = binIDs.get(i);
			int count = counts.get(binID).getCount();
			totalPts = totalPts + count;			
		}
		
		// For ground truth
		int[] bins = new int[totalPts];
		int[] points = new int[totalPts];
		int[] images = new int[totalPts];
		double[] dists = new double[totalPts];
		int gtPos = 0;
		
		// iterate over all returned bins (the bins with 'close' points)				
		for(int i = 0; i < binIDs.size(); i++)
		{
			// for each bin, iterate over its points and vote for 
			// the images where those points belong to
			// also for each point in each bin save bin number, point index within its image and its image - these
			// things will be saved in three arrays: bins, points, images - all of same length
			// when all bins visited, write this info to a file along witht he queryPoint number
			
			// get a bin
			int binID = binIDs.get(i);
			// get its count and index of its first element
			int count = counts.get(binID).getCount();
			int firstElIndex = counts.get(binID).getIndex();
			
			// go into the populated tree, go to the point at index firstElIndex and 
			// iterate over next count points										
			for(int j = firstElIndex; j < firstElIndex + count; j++)
			{
				// get the imageID of the retrieved point				
				int imageID = (populatedTree.getImageIDs())[j];
				
				// vote for the image imageID				
				//if(lastQueryPtVoted[imageID] != queryPtID)
				//{
					votes[imageID]++;
				//	lastQueryPtVoted[imageID] = queryPtID;					
				//}		
					
				// save the info for ground truth file
				bins[gtPos] = populatedTree.getBinIDs()[j];
				points[gtPos] = populatedTree.getPointIDs()[j];
				images[gtPos] = populatedTree.getImageIDs()[j];
				// We record distances for each point even though the
				// the distances for points in the same bin are the same
				// Also the distance returned by best-bin-first is squared and it is claculated from the
				// points that were normalized and scaled to short space - so I will scale them back to 
				// double space by dividing them by Short.MAX^2 AND I will leave them squared.
				dists[gtPos] = distances[i] / (Short.MAX_VALUE * Short.MAX_VALUE); 
				gtPos++;
			}			
		}

		
		if(gtPos != totalPts) 
		{
			try {
				throw new Exception("Thrown in KDTree class. Something is wrong with ground truth info.");
			} catch (Exception e) {
				System.out.println(e.getMessage());
			}
		}

		
//		System.out.println("Writing to the file the following info:");
//		System.out.println(queryImageID);
//		System.out.println(queryPtID);
//		System.out.println(binIDs.size());
//		System.out.println(points.length);		
		
		try {
			writer.write("" + queryImageID); // query image that is used for this query
			writer.newLine();
			writer.write("" + queryPtID); // query point that is used for this query
			writer.newLine();
			writer.write("" + binIDs.size()); // how many bins are retrieved by this query
			writer.newLine();
			writer.write("" + points.length); // how many points are retrieved by this query
			writer.newLine();
			for(int i = 0; i < bins.length; i++) 
			{
				writer.write(bins[i] + " ");
			}
			writer.newLine();
			for(int i = 0; i < points.length; i++) 
			{
				writer.write(points[i] + " ");
			}
			writer.newLine();
			for(int i = 0; i < images.length; i++) 
			{
				writer.write(images[i] + " ");
			}
			writer.newLine();
			for(int i = 0; i < dists.length; i++)
			{
				writer.write(dists[i] + " ");
			}
			writer.newLine();
		} catch (IOException e) {
			e.printStackTrace();
		}
	       
	}
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the query
	 * point and the number of images.
	 * It does not care if one point votes for multiple points in an image and it weights the votes
	 * according to the distance of the query point to the retrieved bins (where retrieved points belong).
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public void voteForPointWithWeights(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Keypoint queryPoint, int queryPtID, Vector<CountIndexPair> counts, double[] votes)
	{				
		// threshold for max distance between this query pt and points that we consider close enough
		double threshold = 100000000000000000000000000.0;
		// threshold on number of bins we want to look into to find reasonably close points to this query point
		int maxBins = 25;
		// variance for weighting the votes
		double variance = 1 * 1;	 // TO DO: fix	(it was 5000*5000)
		
		// array of squared distances from query point to each of the returned bins
		double[] distances = new double[maxBins];		
		// get all the bins tha contain points 'reasonably' close to the query point		
		Vector<Integer> binIDs = this.bestBinFirst(trainedTree, distances, queryPoint, threshold, maxBins);
		
		// iterate over all returned bins (the bins with 'close' points)				
		for(int i = 0; i < binIDs.size(); i++)
		{
			// for each bin, iterate over its points and vote for 
			// the images where those points belong to
			// the closer bins get more weight
			
			// get a bin
			int binID = binIDs.get(i);
			// get its distance from the query point i.e. set the weights for the points in this bin
			// this weight will apply for all the points in this bin
			double distance = distances[i];	// this is the squared distance
			double weight = Math.exp(-distance/variance);
			
			// get its count and index of its first element
			int count = counts.get(binID).getCount();
			int firstElIndex = counts.get(binID).getIndex();
			
			// go into the populated tree, go to the point at index firstElIndex and 
			// iterate over next count points										
			for(int j = firstElIndex; j < firstElIndex + count; j++)
			{
				// get the imageID of the retrieved point				
				int imageID = (populatedTree.getImageIDs())[j];
				
				// vote for the image imageID				
				//if(lastQueryPtVoted[imageID] != queryPtID)
				//{
					votes[imageID] = votes[imageID] + weight;
				//	lastQueryPtVoted[imageID] = queryPtID;					
				//}						
			}			
		}	
	}
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the query
	 * point and the number of images. The last argument 'stats' contains statistics on distribution of points 
	 * in first 4 bins when voting.
	 * It does not care if one point votes for multiple points in an image and it weights the votes
	 * according to the distance of the query point to the retrieved bins (where retrieved points belong).
	 * It also accounts for density of the bins retrieved by the best bin first search in the following way:
	 * it first calculates the density estimate for the bins retrieved by the best bin first search - it has a 
	 * wide normal distribution and it applies it to each point in the retrieved bins. Then it calculates weight 
	 * according to a narrower normal distribution and divides that weight by the density estimate. That way the 
	 * votes in the dense area are downweighted by the amount of density in that area.
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public void voteForPointWithWeightsAndDensityWithStats(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Keypoint queryPoint, int queryPtID, Vector<CountIndexPair> counts, double[] votes, double[] stats)
	{				
		// threshold for max distance between this query pt and points that we consider close enough
		double threshold = 100000000000000000000000000.0;
		// threshold on number of bins we want to look into to find reasonably close points to this query point
		int maxBins = 200;
		// variance for weighting the votes
		double varianceWide = 2000 * 2000;	 
		double varianceNarrow = 100 * 100;
		
		// increment the value of points that voted
		stats[12]++;
		
		// array of squared distances from query point to each of the returned bins
		double[] distances = new double[maxBins];		
		// get all the bins tha contain points 'reasonably' close to the query point		
		Vector<Integer> binIDs = this.bestBinFirst(trainedTree, distances, queryPoint, threshold, maxBins);
		
		// Calculate the valu for the density estimate
		double densityEstimate = 0;	
		
		// iterate over all returned bins (the bins with 'close' points)
		for(int i = 0; i < binIDs.size(); i++)
		{
			// for each bin, calculate its weight and multiply it
			// by the number of points in the bin.
			// then sum over all bins
			
			// get the bin
			int binID = binIDs.get(i);
			// get its distance from the query point i.e. set the weights for the points in this bin
			// this weight will apply for all the points in this bin
			double distance = distances[i];			
			// get the number of points in this bin
			int numPts = counts.get(binID).getCount();
			double weight = Math.exp(-distance/varianceWide);
			densityEstimate = densityEstimate + weight * numPts;
			
			// get statistics about few closest bins
			// 1st bin
			if(i == 0)
			{	
				stats[0] += weight;
			}
			
			// 2nd bin
			if(i == 1)
			{
				stats[1] += weight;
			}
			
			// 3rd bin
			if(i == 2)
			{
				stats[2] += weight;
				
			}
			
			// 4th bin or 200th bin
			if(i == 199)
			{
				stats[3] += weight;
			}
			
		}				
		
		// iterate over all returned bins (the bins with 'close' points)				
		for(int i = 0; i < binIDs.size(); i++)
		{
			// for each bin, iterate over its points and vote for 
			// the images where those points belong to
			// the closer bins get more weight
			
			// get a bin
			int binID = binIDs.get(i);
			// get its distance from the query point i.e. set the weights for the points in this bin
			// this weight will apply for all the points in this bin
			double distance = distances[i];	// this is the squared distance
			double weight = Math.exp(-distance/varianceNarrow) / Math.sqrt(densityEstimate);
			
			// get its count and index of its first element
			int count = counts.get(binID).getCount();
			int firstElIndex = counts.get(binID).getIndex();
			
			// go into the populated tree, go to the point at index firstElIndex and 
			// iterate over next count points										
			for(int j = firstElIndex; j < firstElIndex + count; j++)
			{
				// get the imageID of the retrieved point				
				int imageID = (populatedTree.getImageIDs())[j];
				
				// vote for the image imageID				
					votes[imageID] = votes[imageID] + weight;					
			}
			
//			 get statistics about few closest bins
			// 1st bin
			if(i == 0)
			{	
				stats[4] += weight * Math.sqrt(densityEstimate);
				stats[8] += weight;
			}
			
			// 2nd bin
			if(i == 1)
			{
				stats[5] += weight * Math.sqrt(densityEstimate);
				stats[9] += weight;
			}
			
			// 3rd bin
			if(i == 2)
			{
				stats[6] += weight * Math.sqrt(densityEstimate);
				stats[10] += weight;
			}
			
			// 4th bin or 200th bin
			if(i == 199)
			{
				stats[7] += weight * Math.sqrt(densityEstimate);
				stats[11] += weight;
			}			
		}	
		
	}
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the image
	 * and the number of images in total.
	 * It makes sure that each query point votes only for one point in each image
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public double[] voteForImageOneVotePerImage(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Vector<CountIndexPair> counts, String imageFileName, int numOfImages)
	{
		double[] votes = new double[numOfImages];
		
		// indeces corespond to imageIDs - this array tells us what is the last query point
		// that voted for this image
		// we don't want to allow one query point to vote for many points in an image so that one image 
		// doesn't get too many votes due to many polled points in that image (this would happen if the
		// image contains a pattern for example)
		// we initialize all elements of lastQueryPtVoted to -1
		int[] lastQueryPtVoted = new int[numOfImages];
		for(int i = 0; i < lastQueryPtVoted.length; i++)
		{
			lastQueryPtVoted[i] = -1;
		}
		
		// read the image that we are voting for
		File objFile = new File(imageFileName);		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(objFile)));
			
			ImageKeypoints imageKeypts;
			ImageKeypointPairs imageKeyptPairs;
			Keypoint[] keyptArr;			
			
			// get the all the keypoints or all the keypoint pairs 
			// (depending on if it is .obj or .pair.obj file) from this image
			if(imageFileName.endsWith(".pair.obj"))
			{				
				imageKeyptPairs = (ImageKeypointPairs)in.readObject();
				keyptArr = imageKeyptPairs.getKeyptPairArray();
			}else
			{				
				imageKeypts = (ImageKeypoints)in.readObject();
				keyptArr = imageKeypts.getKeyptArray();
			}							
			
			// vote for each keypoint in this image
			for(int queryPtID = 0; queryPtID < keyptArr.length; queryPtID++)
			{
				Keypoint queryPoint = keyptArr[queryPtID];
				this.voteForPointOneVotePerImage(trainedTree, populatedTree, queryPoint, queryPtID, counts, votes, lastQueryPtVoted);
			}
			in.close();
			
		}catch (IOException ioe)
		{
		ioe.printStackTrace();
		}
		catch (ClassNotFoundException cnfe)
		{ 
			cnfe.printStackTrace();
		}		
		return votes;
	}
	
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the image
	 * and the number of images in total.
	 * It lets each query point votes for all polled points in an image and it weights those points
	 * according to the distance of the query point to the retrieved bins (where retrieved points belong).
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public double[] voteForImageWithWeights(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Vector<CountIndexPair> counts, String imageFileName, int numOfImages)
	{
		double[] votes = new double[numOfImages];
		
		// indeces corespond to imageIDs - this array tells us what is the last query point
		// that voted for this image
		// we don't want to allow one query point to vote for many points in an image so that one image 
		// doesn't get too many votes due to many polled points in that image (this would happen if the
		// image contains a pattern for example)
		// we initialize all elements of lastQueryPtVoted to -1
		//int[] lastQueryPtVoted = new int[numOfImages];
		//for(int i = 0; i < lastQueryPtVoted.length; i++)
		//{
		//	lastQueryPtVoted[i] = -1;
		//}
		
		// read the image that we are voting for
		File objFile = new File(imageFileName);		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(objFile)));
			
			ImageKeypoints imageKeypts;
			ImageKeypointPairs imageKeyptPairs;
			Keypoint[] keyptArr;			
			
			// get the all the keypoints or all the keypoint pairs 
			// (depending on if it is .obj or .pair.obj file) from this image
			if(imageFileName.endsWith(".pair.obj"))
			{				
				imageKeyptPairs = (ImageKeypointPairs)in.readObject();
				keyptArr = imageKeyptPairs.getKeyptPairArray();
			}else
			{				
				imageKeypts = (ImageKeypoints)in.readObject();
				keyptArr = imageKeypts.getKeyptArray();
			}							
			
			// vote for each keypoint in this image
			for(int queryPtID = 0; queryPtID < keyptArr.length; queryPtID++)
			{
				Keypoint queryPoint = keyptArr[queryPtID];
				this.voteForPointWithWeights(trainedTree, populatedTree, queryPoint, queryPtID, counts, votes);
			}
			in.close();
			
		}catch (IOException ioe)
		{
		ioe.printStackTrace();
		}
		catch (ClassNotFoundException cnfe)
		{ 
			cnfe.printStackTrace();
		}		
		return votes;
	}
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the image
	 * and the number of images in total.
	 * It lets each query point votes for all polled points in an image and it weights those points
	 * according to the distance of the query point to the retrieved bins (where retrieved points belong).
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public double[] voteForImageWithWeightsAndDensityWithStats(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Vector<CountIndexPair> counts, String imageFileName, int numOfImages, double[] stats)
	{
		double[] votes = new double[numOfImages];
		
		// indeces corespond to imageIDs - this array tells us what is the last query point
		// that voted for this image
		// we don't want to allow one query point to vote for many points in an image so that one image 
		// doesn't get too many votes due to many polled points in that image (this would happen if the
		// image contains a pattern for example)
		// we initialize all elements of lastQueryPtVoted to -1
		//int[] lastQueryPtVoted = new int[numOfImages];
		//for(int i = 0; i < lastQueryPtVoted.length; i++)
		//{
		//	lastQueryPtVoted[i] = -1;
		//}
		
		// read the image that we are voting for
		File objFile = new File(imageFileName);		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(objFile)));
			
			ImageKeypoints imageKeypts;
			ImageKeypointPairs imageKeyptPairs;
			Keypoint[] keyptArr;			
			
			// get the all the keypoints or all the keypoint pairs 
			// (depending on if it is .obj or .pair.obj file) from this image
			if(imageFileName.endsWith(".pair.obj"))
			{				
				imageKeyptPairs = (ImageKeypointPairs)in.readObject();
				keyptArr = imageKeyptPairs.getKeyptPairArray();
			}else
			{				
				imageKeypts = (ImageKeypoints)in.readObject();
				keyptArr = imageKeypts.getKeyptArray();
			}							
			
			// vote for each keypoint in this image
			for(int queryPtID = 0; queryPtID < keyptArr.length; queryPtID++)
			{
				Keypoint queryPoint = keyptArr[queryPtID];
				this.voteForPointWithWeightsAndDensityWithStats(trainedTree, populatedTree, queryPoint, queryPtID, counts, votes, stats);
			}
			in.close();
			
		}catch (IOException ioe)
		{
		ioe.printStackTrace();
		}
		catch (ClassNotFoundException cnfe)
		{ 
			cnfe.printStackTrace();
		}		
		return votes;
	}
	
	
	/*
	 * This method does the voting given the populated tree, the trained tree, the image
	 * and the number of images in total.
	 * It lets each query point votes for all polled points in an image
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 */	
	public double[] voteForImage(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Vector<CountIndexPair> counts, String imageFileName, int numOfImages)
	{
		double[] votes = new double[numOfImages];
		
		// indeces corespond to imageIDs - this array tells us what is the last query point
		// that voted for this image
		// we don't want to allow one query point to vote for many points in an image so that one image 
		// doesn't get too many votes due to many polled points in that image (this would happen if the
		// image contains a pattern for example)
		// we initialize all elements of lastQueryPtVoted to -1
		//int[] lastQueryPtVoted = new int[numOfImages];
		//for(int i = 0; i < lastQueryPtVoted.length; i++)
		//{
		//	lastQueryPtVoted[i] = -1;
		//}
		
		// read the image that we are voting for
		File objFile = new File(imageFileName);		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(objFile)));
			
			ImageKeypoints imageKeypts;
			ImageKeypointPairs imageKeyptPairs;
			Keypoint[] keyptArr;			
			
			// get the all the keypoints or all the keypoint pairs 
			// (depending on if it is .obj or .pair.obj file) from this image
			if(imageFileName.endsWith(".pair.obj"))
			{				
				imageKeyptPairs = (ImageKeypointPairs)in.readObject();
				keyptArr = imageKeyptPairs.getKeyptPairArray();
			}else
			{				
				imageKeypts = (ImageKeypoints)in.readObject();
				keyptArr = imageKeypts.getKeyptArray();
			}							
			
			// vote for each keypoint in this image
			for(int queryPtID = 0; queryPtID < keyptArr.length; queryPtID++)
			{
				Keypoint queryPoint = keyptArr[queryPtID];
				this.voteForPoint(trainedTree, populatedTree, queryPoint, queryPtID, counts, votes);
			}
			in.close();
			
		}catch (IOException ioe)
		{
		ioe.printStackTrace();
		}
		catch (ClassNotFoundException cnfe)
		{ 
			cnfe.printStackTrace();
		}		
		return votes;
	}
	

	/*
	 * This method does the voting given the populated tree, the trained tree, the image
	 * and the number of images in total.
	 * It lets each query point votes for all polled points in an image
	 * It returns the vector with that number of elements equal to the number of images where
	 * element at index i represents the image with imageID = i.
	 * This method is used to produce some info for ground truth - it calls voteForPointWithGroundTruth
	 * method that for each point writes some info to a file (one file per image).
	 */	
	public double[] voteForImageWithGroundTruth(Vector<Node> trainedTree, PopulatedKDTree populatedTree,
			Vector<CountIndexPair> counts, String imageFileName, int queryPtID, int numOfImages, BufferedWriter writer)
	{
		double[] votes = new double[numOfImages];
		
		// indeces corespond to imageIDs - this array tells us what is the last query point
		// that voted for this image
		// we don't want to allow one query point to vote for many points in an image so that one image 
		// doesn't get too many votes due to many polled points in that image (this would happen if the
		// image contains a pattern for example)
		// we initialize all elements of lastQueryPtVoted to -1
		//int[] lastQueryPtVoted = new int[numOfImages];
		//for(int i = 0; i < lastQueryPtVoted.length; i++)
		//{
		//	lastQueryPtVoted[i] = -1;
		//}
		
		// read the image that we are voting for
		File objFile = new File(imageFileName);		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(objFile)));
			
			ImageKeypoints imageKeypts;
			ImageKeypointPairs imageKeyptPairs;
			Keypoint[] keyptArr;			
			
			String queryImageID = null;
			
			// get the all the keypoints or all the keypoint pairs 
			// (depending on if it is .obj or .pair.obj file) from this image
			if(imageFileName.endsWith(".pair.obj"))
			{								
				imageKeyptPairs = (ImageKeypointPairs)in.readObject();
				imageKeyptPairs.normalizeKeypointPairs();
				keyptArr = imageKeyptPairs.getKeyptPairArray();
				queryImageID = imageFileName.substring(imageFileName.length() - 1 - 8 - 5, 
						imageFileName.length() - 1 - 8);
			}else
			{								
				imageKeypts = (ImageKeypoints)in.readObject();
				imageKeypts.normalizeKeypoints();
				keyptArr = imageKeypts.getKeyptArray();
				queryImageID = imageFileName.substring(imageFileName.length() - 1 - 3 - 5, 
						imageFileName.length() - 1 - 3);
			}		
			
			//System.out.println("Num points in this image: " + keyptArr.length);
			
			in.close();
			
//			// Info about voting (bins, images, point indeces) for ground truth
//			//String outputFileName = "votingInfoIm" + Integer.parseInt(imageNum) + ".txt";
//			String outputFileName = "votingInfo2" + ".txt";
//			
//			// Open the the file named outputFileName and and pass BufferWrite instance when voting for each point	
//			BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));
			
			
			// vote for each keypoint in this image and write some ground truth info to outputFile
			// write to the file for 3 query points only, count with the couter
//			int counter = 0;
//			for(int queryPtID = 0; queryPtID < keyptArr.length; queryPtID++)
//			{
//				Keypoint queryPoint = keyptArr[queryPtID];				
//				counter++;						
//				this.voteForPointWithGroundTruth(trainedTree, populatedTree, queryPoint, queryPtID, counts, votes, writer);
//				// done voting and writing - close the writer
//				if(counter == 2)
//				{
//					writer.flush();
//					writer.close();
//					System.exit(0);
//				}
//			}
			
			Keypoint queryPoint = keyptArr[queryPtID];									
			this.voteForPointWithGroundTruth(trainedTree, populatedTree, Integer.parseInt(queryImageID), queryPoint, queryPtID, counts, votes, writer);
			
			
		}catch (IOException ioe)
		{
		ioe.printStackTrace();
		}
		catch (ClassNotFoundException cnfe)
		{ 
			cnfe.printStackTrace();
		}		
		return votes;
	}
	
	/*
	 * The method that calculates and returns the dimension bounds for the given data set.
	 * To account for the points lower bound should be lower then the lowest point for 25% 
	 * of the difference between the lowest and highest point and the higher bound should 
	 * be 25% higher then the highest point. 
	 */
	public double[] getDimensionBounds(short[] dataSet)
	{
		// array for dimension bounds contains in its lower/higher part lower/higher bounds
		double[] dimBounds = new double[Keypoint.DESCRIPTOR_LENGTH + Keypoint.DESCRIPTOR_LENGTH];
		
		//array of minimums accross all the points
		double[] mins = new double[Keypoint.DESCRIPTOR_LENGTH];
		double[] maxs = new double[Keypoint.DESCRIPTOR_LENGTH];
		
		for(int i = 0; i < Keypoint.DESCRIPTOR_LENGTH; i++)
		{
			double dimMin = dataSet[i];	//initialize min to the first element of this dimension
			// loop over all points and find min
			for(int j = i; j < dataSet.length; j = j + Keypoint.DESCRIPTOR_LENGTH)
			{											
				if(dataSet[j] < dimMin)
				{
					dimMin = dataSet[j];
				}										
			}
			mins[i] = dimMin; 
						
			double dimMax = dataSet[i];	//initialize max to the first element of this dimension
			// loop over all points and find max
			for(int j = 0; j < dataSet.length; j = j + Keypoint.DESCRIPTOR_LENGTH)
			{				
				if(dataSet[j] > dimMax)
				{
					dimMax = dataSet[j];
				}								 	
			}
			maxs[i] = dimMax;					
		}
		
		Arrays.sort(mins);
		Arrays.sort(maxs);
		double min = mins[0];
		double max = maxs[Keypoint.DESCRIPTOR_LENGTH -1];
		double diff = max - min;
		
		for (int i = 0; i < Keypoint.DESCRIPTOR_LENGTH; i++)
		{
			dimBounds[i] = min - (diff * 0.15);
			dimBounds[i + Keypoint.DESCRIPTOR_LENGTH] = max + (diff * 0.15);
		}
		
		return dimBounds;
	}
	
	
	
	/*
	 * Helping method - used for printing.
	 */
	public void printNodeContent(Vector tree, Node node)
	{
		if(node instanceof KDNode)
		{			
			int leftChildIndex = ((KDNode) node).getLeftChild();
			int rightChildIndex = ((KDNode) node).getRightChild();
			Node leftChild = (Node) tree.get(leftChildIndex);
			Node rightChild = (Node) tree.get(rightChildIndex);
			System.out.println("Dimension: " + ((KDNode) node).getDimension());
			System.out.println("Split point: " + ((KDNode) node).getSplitPt());
			double[] bounds = ((KDNode) node).getDimensionsBounds();
			System.out.println("Bounds for this node: " + bounds[KDNode.LO] + " and "+ bounds[KDNode.HI]);			
			System.out.println("Left child is at index: " + leftChildIndex);
			System.out.println("Right child is at index: " + rightChildIndex);
			System.out.println("\nLeft child info: ");
			printNodeContent(tree, leftChild);
			System.out.print("\nRight child info: \n");
			printNodeContent(tree, rightChild);
		}else
		{
			System.out.println("Node has too few elements - there is no split.");
		}
	}
}