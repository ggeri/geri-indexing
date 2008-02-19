import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Vector;

/**
 * The instance of this class is an object that contains all the
 * keypoint pairs that belong to one image.
 * @author geri
 *
 */
public class ImageKeypointPairs implements Serializable
{
	public static void main(String[] args)
	{
		NumberFormat formatter = new DecimalFormat("00000");
    	String objFileName =KeypointPairsExtraction.imagePath + "ukbench" + formatter.format(0) + ".obj";
    	
    	ImageKeypointPairs pairs = new ImageKeypointPairs(objFileName);
    	
    	String pairobjFileName = KeypointPairsExtraction.imagePathNoBackup + "ukbench" + formatter.format(0) + ".pair.obj";
    	File pairobjFile = new File(pairobjFileName);    		    	 	
        
		try
		{
			pairobjFile.delete();
			pairobjFile.createNewFile();
			ObjectOutputStream out = new ObjectOutputStream(new
					BufferedOutputStream(new FileOutputStream(pairobjFile)));
			out.writeObject(pairs);		        	    		        	    	        	  
    	    out.close();
		}
		catch (SecurityException se)
		{
			se.printStackTrace();
		}
		catch (FileNotFoundException fnfe)
		{
			fnfe.printStackTrace();
		}
		catch (IOException ioe)
		{
			ioe.printStackTrace();
		}    	
	}
	
	
	/*
	 * For serialization.
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * How many pairs do we want - we want a multiple of the number of keypoints, i.e. we want 
	 * NUM_PAIRS * number_of_keypoints
	 */
	private final int NUM_PAIRS = 1;
	
	/**
	 * Number of nearest neighbours we are interested in; i.e. we consider some number x of the closest points to a 
	 * keypoint; we pair that one point with each of its nearest neighbours and calculate the score. Then we 
	 * choose a number NUM_PAIRS of the pairs with the best score for each of the keypoints. The number x is 
	 * given here as a constant.
	 */
	private final int NUM_CLOSEST_PTS = 30;
	
	
	/**
	 * The array that contains all the keypoint pairs of this image
	 */
	private KeypointPair[] keyptPairArray;
	
	
	/**
	 * Helping constructor - doesn't do anything
	 */
	public ImageKeypointPairs()
	{
		this.keyptPairArray = null;
	}
	
	/**
	 * The constructor that creates an instance of this class. The parameter fileName
	 * is the string that represents the name of the pkey or obj image file. The constructor invokes 
	 * the collectKeypointPairs method which reads the contents
	 * of the file, creates the keypoint pairs and stores them in the keyptPairArray attribute of this class.
	 * @param fileName
	 */
	public ImageKeypointPairs(String fileName)
	{
		this.collectKeypointPairs(fileName);		
		
	}
		
	
	
	/**
	 * This method takes the passed file name, opens the file, reads it line by line and 
	 * creates keypoint pairs and populates with them the keyptPairArray attribute of this class.
	 * @param keyFileName
	 */
	public void collectKeypointPairs(String fileName)
	{		
		// first we populate the keyptArray with singletons and then we pair the points
		
		KeypointSingleton[] keyptArray = null;
		if(fileName.endsWith(".pkey"))
		{
			keyptArray = this.getSingletonsFromPkey(fileName);
		}else if(fileName.endsWith(".obj"))
		{
			keyptArray = this.getSingletonsFromObj(fileName);
			
		}else
		{
			System.out.println("Wrong filename extension");
			System.exit(1);
		}

		
		KeypointPair[] keyptPairArrayTemp = new KeypointPair[this.NUM_PAIRS * keyptArray.length];		
		
		// iterate over all the singletons in keyptArray, take one by one and pair it with 30 of its neighbours 			
		
		int lastPosition = 0;
		for(int i = 0; i < keyptArray.length; i++)
		{
			// take one keypt from the array and find the indeces of NUM_CLOSEST_PTS its closest neighbour-points
			KeypointSingleton keypt = keyptArray[i];			
			int[] closestPtsIndeces = this.closestPoints(keypt, i, keyptArray);
			
			// for each pair (keypt at i, one of closest pts) get the score and save it in scores array; 
			// the score at index j corresponds to the score between this point we are looking at 
			// right now (at position i in keyptArray) and the point at index j in closestPts array
			double[] scores = new double[closestPtsIndeces.length];						
			
			for(int j = 0; j < closestPtsIndeces.length; j++)
			{
				int keypt2Indx = closestPtsIndeces[j];
				KeypointSingleton keypt2 = keyptArray[keypt2Indx];
				double score = this.evaluateScore(keypt, keypt2);
				scores[j] = score;			
			}
			
			// now for each point in keyptArray choose NUM_PAIRS keypoints that scored the best, to create final pairs of them
			// we save the indeces of these closest neighbours in pointTwoIndeces (we can forget now closestPtsIndeces)
			// if the number of elements in keyptArray is less then NUM_PAIRS, then we get only the |keyptArray| highest scores
			int[] pointTwoIndeces;
			if(scores.length <= NUM_PAIRS)
			{
				pointTwoIndeces = this.highestScoresIndeces(scores, scores.length);
			}else
			{
				pointTwoIndeces = this.highestScoresIndeces(scores, NUM_PAIRS);
			}			
			
			// now we have the first and second point of each pair - the first point is the keypt at position i in the 
			// keyptArray and the second point is each point from keyptArray that is located in index from pointTwoIndeces
			// we can now create KeypointPair for each 
			// pair and store it into the keyptPairArray
			
			// first point of the pair
			KeypointSingleton ptOne = keypt;
			double[] ptOneLoc = ptOne.getKeyptLocation();
			short[] ptOneDesc = ptOne.getDescriptor();
			int k;
			for(k = 0; k < pointTwoIndeces.length; k++)
			{
				// second point of the pair
				KeypointSingleton ptTwo = keyptArray[pointTwoIndeces[k]];
				double[] ptTwoLoc = ptTwo.getKeyptLocation();
				short[] ptTwoDesc = ptTwo.getDescriptor();
				//create pair descriptor for this pair
				short[] pairDesc = new short[KeypointSingleton.DESCRIPTOR_LENGTH];
				// for now create pairDesc out of 18 first decsription values of each point in
				// that pair - CHANGE it later to 16 from each point + some spacial relations
				int halfWay = (int)(Math.floor(KeypointSingleton.DESCRIPTOR_LENGTH / 2));
				for(int l = 0; l < halfWay; l++)
				{
					pairDesc[l] = ptOneDesc[l];
					pairDesc[l+halfWay] = ptTwoDesc[l];
				}			
				KeypointPair pair = new KeypointPair(ptOneLoc, ptTwoLoc, pairDesc);						
				
				keyptPairArrayTemp[lastPosition + k] = pair;				
			}			
			lastPosition = lastPosition + k;
		}
		// trim the empty spots in keyptPairArray since before we intitialized it to 
		// NUM_PAIRS * keyptArray.length and we might end up having less then that number
		// copy elemnts from keyptPairArray that are not null to this.keyptPairArray			
		
		if(lastPosition == 0)
		{
			this.keyptPairArray = new KeypointPair[0];
		}else
		{
			this.keyptPairArray = new KeypointPair[lastPosition - 1];
		}
		
		System.out.println("Number of singletons is: " + keyptArray.length);
		System.out.println("Number of pairs is: " + this.keyptPairArray.length);
		
		for(int i = 0; i < keyptPairArray.length; i++)
		{
			if(keyptPairArrayTemp != null)
			{
				this.keyptPairArray[i] = keyptPairArrayTemp[i];
			}else
			{
				System.out.println("Keypt at position " + i + " in keypointPairArrayTemp is null!!!");
				System.exit(1);
			}
		}		
	}
	
	/**
	 * This method reads the file named by fileName and gets the keypoint array
	 * from the object read. It them populates with it the passed keyptArray.
	 * @param fileName - the file to read the object from
	 * @param keyptArray - the array of keypoints to populate wtih singleton keypoints
	 */
	
	private KeypointSingleton[] getSingletonsFromPkey(String fileName)
	{
		KeypointSingleton[] keyptArray = null;
		try
		{
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String line = in.readLine();
			String[] parsedLine = line.split("\\s");
			int numOfKeypts = Integer.parseInt(parsedLine[0]);
			// initialize the keypoint array
			keyptArray = new KeypointSingleton[numOfKeypts];
			for(int i = 0; i < numOfKeypts; i++)				
			{
				KeypointSingleton keypt = new KeypointSingleton();
				//these are the location array values
				line = in.readLine();			
				parsedLine = line.split("\\s");
				double[] loc = new double[4];
				for (int l = 0; l < 4; l++)
				{
					loc[l] = Double.parseDouble(parsedLine[l]);
				}
				//set it in this instance of Keypoint class
				keypt.setKeyptLocation(loc);

				// these are the descriptor array values
				String[] wholeParsedLine = new String[KeypointSingleton.DESCRIPTOR_LENGTH];
				
				for (int j = 0; j < 3; j++)
				{
					line = in.readLine();
					//parsedLine array is of length 13 but it has 12 tokens
					parsedLine = line.split("\\s");			
					int off = j * 12;									
					for(int m = 1; m < parsedLine.length; m++)
					{						
						int off2 = off + m -1;
						wholeParsedLine[off2] = parsedLine[m];	
					}												
				}
				short[] descr = new short[KeypointSingleton.DESCRIPTOR_LENGTH];				
				for (int k = 0; k < KeypointSingleton.DESCRIPTOR_LENGTH; k++)
				{										
					descr[k] = Short.parseShort(wholeParsedLine[k]);
				}
				//set the vector in the object				
				keypt.setDescriptor(descr);

				//add the pt to array of points	
				keyptArray[i] = keypt;
			}
			in.close();
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return keyptArray;
	}
	
	
	/**
	 * This method reads the file named by fileName and gets the keypoint array
	 * from the object read. It them populates with it the passed keyptArray.
	 * @param fileName - the file to read the object from
	 * @param keyptArray - the array of keypoints to populate wtih singleton keypoints
	 */
	
	private KeypointSingleton[] getSingletonsFromObj(String fileName)
	{
		KeypointSingleton[] keyptArray = null;
		try
		{
			File file = new File(fileName);
			ObjectInputStream in = new ObjectInputStream(new
			  		BufferedInputStream(new FileInputStream(file)));
		    ImageKeypoints keypts = (ImageKeypoints)in.readObject();			  			  
		    keyptArray = keypts.getKeyptArray();
			in.close();
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		} catch (ClassNotFoundException e) {			
			e.printStackTrace();
		}			
		
		return keyptArray;
	}
	
	/**
	 * This method takes array of scores and an int num_scores. It will find
	 * the num_scores highest scores in this scores array and return their indeces in 
	 * the 'highestScoresIndeces' array.
	 * @param scores - array of double
	 * @param num_scores - number of highest scores we are looking for
	 * @return highestScoresIndeces - indeces of num_scores highest scores in 'scores' array
	 * TESTED - works fine
	 */
	private int[] highestScoresIndeces(double[] scores, int num_scores)
	{
		// this array will contain indeces of highest scores in 'scores' array
		int[] highestScoresIndeces = new int[num_scores];		
		for(int i = 0; i < num_scores; i++)
		{
			//find max score
			double maxScore = -1;			
			int maxScoreIndx = -1;			
			for(int k = 0; k < scores.length; k++)
			{					
				if(scores[k] > maxScore)
				{
					maxScore = scores[k];
					maxScoreIndx = k;
				}
			}
			
			highestScoresIndeces[i] = maxScoreIndx;
			// reassign this score to the min value so that we can find the next max score
			scores[maxScoreIndx] = -1;
		}
		return highestScoresIndeces;
	}
	
	
	/**
	 * Given a keypt, an array of keypoints and a number num, this function returns the indeces
	 * of the num closest keypoints (from the given array) to the keypt.
	 * I.e. the function calculates the distance from keypt to each keypoint in the passed array and returns
	 * the indeces of the num points in that array that have the smallest distance.
	 * TESTED - works fine
	 */
	public int[] closestPoints(KeypointSingleton keypt, int keyptIndx, KeypointSingleton[] keyptArray) // maybe this method should be private?
	{		
		// this array contains the distance from keypt to each point in keyptArray - index i in 
		// dist array corresponds to distance from keypt to Keypoint on index i in keyptArray
		double[] dists = new double[keyptArray.length];
		double dist;
		//get all the distances and put them into dist array
		for(int i = 0; i < keyptArray.length; i++)
		{
			dist = spacialDist(keyptArray[i], keypt);
			dists[i] = dist;
		}
		
		// choose NUM_CLOSEST_PTS smallest distances and put the indeces of those keypoints into closestPtsIndeces array
		// if the length of the keyptArray is less then NUM_CLOSEST_PTS, then we find only |keyptArray| closest points
		int[] closestPtsIndeces;
		if(keyptArray.length <= this.NUM_CLOSEST_PTS)
		{
			closestPtsIndeces = new int[keyptArray.length-1];
			for(int i = 0; i < (keyptArray.length - 1); i++) // need -1 since not taking the distance of the pt to itself
			{
				// find min distance
				double minDist = Double.MAX_VALUE;
				int minDistIndx = -1;
				int j;
				for(j = 0; j < dists.length; j++)
				{
					if(j != keyptIndx)	// we don't want the distance of keypt to itself
					{
						if(dists[j] < minDist)
						{
							minDist = dists[j];
							minDistIndx = j;
						}
					}
				}
				
				closestPtsIndeces[i] = minDistIndx;			
				
				// reassign this dist to the max value so that we can find the next smallest distance
				dists[minDistIndx] = Double.MAX_VALUE;
			}
		}else
		{
			closestPtsIndeces = new int[this.NUM_CLOSEST_PTS];
			for(int i = 0; i < this.NUM_CLOSEST_PTS; i++)
			{
				// find min distance
				double minDist = Double.MAX_VALUE;
				int minDistIndx = -1;
				int j;
				for(j = 0; j < dists.length; j++)
				{
					if(j != keyptIndx)	// we don't want the distance of keypt to itself
					{
						if(dists[j] < minDist)
						{
							minDist = dists[j];
							minDistIndx = j;
						}
					}
				}
				
				closestPtsIndeces[i] = minDistIndx;			
				
				// reassign this dist to the max value so that we can find the next smallest distance
				dists[minDistIndx] = Double.MAX_VALUE;
			}
		}
		return closestPtsIndeces;
	}
	
	
	/**
	 * This method gives the spacial distance between two passed points keypt1 and keypt2.
	 * TESTED - works fine
	 */
	private double spacialDist(KeypointSingleton keypt1, KeypointSingleton keypt2)
	{
		double x1 = (keypt1.getKeyptLocation())[1];
		double y1 = (keypt1.getKeyptLocation())[0];
		double x2 = (keypt2.getKeyptLocation())[1];
		double y2 = (keypt2.getKeyptLocation())[0];
		double distSq = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
		return Math.sqrt(distSq);
	}
	
	/**
	 * This function, given two scales, returns the score in interval [0,1]. pointOne has to be > 0 
	 * and pointOne has to be < pointTwo, otherwise the function returns -1;
	 * It first divides smaller scale by the larger one and gets some value that is always between 0 and 1.
	 * Then for all the values that are <= to pointOne returns 0 and for all the values >= pointTwo
	 * returns 1. For the values that are between pointOne and pointTwo, it returns a value according 
	 * to the given function (it is the line equation).
	 * @param scaleOne - scale of the first keypoint
	 * @param scaleTwo -  scale of the second keypoint
	 * @param pointOne - the point where the score starts increasing from 0 towards 1
	 * @param pointTwo -  after this point, the score is 1
	 * @return - the score calculated according to the given algorithm
	 * TESTED - works
	 */
	public double funScale(double scaleOne, double scaleTwo, double pointOne, double pointTwo)
	{
		if((pointOne <= 0) || (pointOne >= pointTwo))
		{
			return -1;
		}
		double scaleDiv = Math.min(scaleOne, scaleTwo) / Math.max(scaleOne, scaleTwo);				
		
		if(scaleDiv <= pointOne)
		{
			return 0.0;
		}else if(scaleDiv >= pointTwo)
		{
			return 1.0;
		}else
		{
			double pt1x = pointOne;
			double pt1y = 0;
			double pt2x = pointTwo;
			double pt2y = 1;
			double slope = (pt2y - pt1y) / (pt2x - pt1x);
			double y = slope * (scaleDiv - pt1x);
			return y;
		}
	}
	
	/**
	 * This function, given two feature vectors, returns the score in interval [0,1]. point has to be >=0 
	 * otherwise the method returns -1;
	 * It first calculates the distance between the two vectors.
	 * Then for all the values that are >= then point returns 1 and for all the values
	 * < then point, it returns a value according to the given function (it is the line
	 * equation).
	 * @param descOne - feature vector of keypoint one
	 * @param descTwo - feature vector of keypoint two
	 * @param point - if distance between the two keypoints (in feature space) is larger then the value
	 * of point parameter, then the function score is 1, otherwise it calculates the score according 
	 * to the above algorithm
	 * @return -  the score calculated according to the given algorithm
	 * TESTED - works fine
	 */
	public double funFeatureDist(short[] descOne, short[] descTwo, double point)
	{
		if(point <= 0)
		{
			return -1;
		}
		
		double featureDistSq = 0;
		for(int i = 0; i < descOne.length; i++)
		{
			featureDistSq =featureDistSq + ((descOne[i] - descTwo[i]) * (descOne[i] - descTwo[i])); 
		}
		// distance of the two points in feature space
		double featureDist = Math.sqrt(featureDistSq);
		
		if(featureDist >= point)
		{
			return 1;
		}else
		{
			double pt1x = 0;
			double pt1y = 0.25;
			double pt2x = point;
			double pt2y = 1;
			double slope = (pt2y - pt1y) / (pt2x - pt1x);
			double y = slope * (featureDist - pt1x) + pt1y;
			return y;
		}
	}
	
	/**
	 * This function, given two location points, three thresholds and the scales of the
	 * two points, returns the score in interval [0,1]. 
	 * pointOne has to be >=0 and <= then pointTwo, otherwise the function returns -1;
	 * What does the method do:
	 * It first calculates the distance between the two points and then it normalizes it by the min of 
	 * the two scales.
	 * Then for all the values that are >= then pointOne and <= pointTwo returns 1, for all the 
	 * values > pointThree it returns 0 and for points < pointOne and between pointTwo and pointThree
	 * it returns a value according to the given functions (the simple line equations).
	 * @param xOne - x coordinate of the first point
	 * @param yOne - y coordinate of the first point
	 * @param xTwo - x coordinate of the second point
	 * @param yTwo - y coortinate of the second point
	 * @param pointOne - one of the three end points on the graph
	 * @param pointTwo - one of the three end points on the graph
	 * @param pointThree - one of the three end points on the graph
	 * @param scaleOne - scale of the first keypoint
	 * @param scaleTwo - scale of the second keypoint
	 * @return
	 * TESTED - works fine
	 */
	public double funSpacialDist(double xOne, double yOne, double xTwo, double yTwo, double pointOne, double pointTwo, 
			double pointThree, double scaleOne, double scaleTwo)
	{
		if( (pointOne < 0) || (pointOne > pointTwo) )
		{
			return -1;
		}
		
		double spaceDistSq = ((xOne - xTwo)*(xOne - xTwo)) + ((yOne - yTwo)*(yOne - yTwo));
		double scaleMin = Math.min(scaleOne, scaleTwo);
		double spaceDist = Math.sqrt(spaceDistSq) / scaleMin;		
		
		if(spaceDist < pointOne)
		{
			double pt1x = 0;
			double pt1y = 0;
			double pt2x = pointOne;
			double pt2y = 1;
			double slope = (pt2y - pt1y) / (pt2x - pt1x);
			double y = slope * (spaceDist - pt1x);
			return y;		
			
		}else if(spaceDist > pointThree)
		{
			return 0;								
		}else if(spaceDist > pointTwo)
		{
			double pt2x = pointTwo;
			double pt2y = 1;
			double pt3x = pointThree;
			double pt3y = 0;
			double slope = (pt3y - pt2y) / (pt3x - pt2x);
			double y = slope * (spaceDist - pt3x);
			return y;
		}else
		{
			return 1;
		}	
	}
	
	/**
	 * This function takes two keypoints, evaluates for them scale, feature and distance 
	 * scores and combines the three into a final score for this pair of keypoints.
	 * @param keyptOne - the first keypoint of this pair
	 * @param keyptTwo - the second keypoint of this pair
	 * @return - the score for this pair of keypoints
	 */
	public double evaluateScore(KeypointSingleton keyptOne, KeypointSingleton keyptTwo)
	{
		short[] descOne = keyptOne.getDescriptor();
		double[] locOne = keyptOne.getKeyptLocation();
		short[] descTwo = keyptTwo.getDescriptor();
		double[] locTwo = keyptTwo.getKeyptLocation();
		
		//position in the loc vector
		int x = 1, y = 0, scale = 2;
		
		double scaleScore = this.funScale(locOne[scale], locTwo[scale], 0.1, 0.5);
		double featureScore = this.funFeatureDist(descOne, descTwo, 10000);
		double distScore = this.funSpacialDist(locOne[x], locOne[y], locTwo[x], locTwo[y], 8, 16, 30, 
				locOne[scale], locTwo[scale]);
		
		double score = scaleScore * featureScore * distScore;
		return score;
	}
	
	/**
	 * Get method for the attribute keyptPairArray
	 */
	public KeypointPair[] getKeyptPairArray()
	{
		return this.keyptPairArray;
	}
	
}