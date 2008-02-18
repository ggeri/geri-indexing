import java.io.*;
import java.text.*;
import java.util.*;


public class IndexingTest1 {

/**
 * This attribute is used for training the tree. 
 * true = the tree will be trained 
 * false = the tree will be read from the disk - if there is a file with the specified name on the disk.
 */
private static boolean train = false;	
private static boolean populate = false;

private static final int NUM_NODES_IN_TRAINED_TREE = 3 * IndexingTest1.NUM_LEAFS_IN_TRAINED_TREE;
private static final int NUM_LEAFS_IN_TRAINED_TREE = 1000000;		
private static final int NUM_PTS_TRAINING = 10000000;	
private static final int NUM_IMAGES_TRAINING = 10200;
private static final int NUM_IMAGES_POPULATING = 10200; 	// works with 1000 if heap set to 6144max & w/ flag -XX:-UseGCOverheadLimit
private static final int NUM_RUNS = 100;
	
	public static void main(String[] args)
	{	
		String stars = "\n**************************************************************************************\n";
				
		
		/* Instantiate the KeypointsExtraction class - this converts jpg files to pgm, extracts the keypoints and saves them
		 * to .key files and finally reads thesefiles and creates serialized javaobjects on the disk.
		 */
		System.out.println(stars + "\nCreate .obj file for each image and store them to the disk. ");
		//KeypointsExtraction keyptExtract = new KeypointsExtraction();		
		
		/*
		 * Get an instance of the KDTree class
		 */
		KDTree kdTreeInstance = new KDTree();
		Vector<Node> trainedTree = new Vector<Node>(IndexingTest1.NUM_NODES_IN_TRAINED_TREE); 
		
		int dataSetArrLength = IndexingTest1.NUM_PTS_TRAINING * KeypointSingleton.DESCRIPTOR_LENGTH;  
		
		short[] dataSetTraining = new short[dataSetArrLength];


	    if(train)
	    {
	    	/*
			 * Get the set of points for training. We want max 2 000 000 points - we go over images
			 * starting at the first one till we get 2 000 000 points; we take max 1000 images.
			 */
			System.out.println(stars + "\nCreate dataset for training the tree. ");		
			
			long timeStart = System.currentTimeMillis();			
			try
			{
				int position = 0;
				int image = 0;
				for(image = 0; image < IndexingTest1.NUM_IMAGES_TRAINING && position < (dataSetArrLength); image++)
				{					
					short[] dataSet = KeypointsExtraction.getDataSet(image);										
					
					if(dataSet.length > dataSetArrLength - position)
					{
						System.arraycopy(dataSet, 0, dataSetTraining, position, (dataSetArrLength - position));
						position = position + (dataSetArrLength - position);
					}else
					{
						System.arraycopy(dataSet, 0, dataSetTraining, position, dataSet.length);
						position = position + dataSet.length;
					}						
				}
				System.out.println("Number of images used for the dataSetTreining is " + image);
				
			} catch(Exception e)
			{				
				System.err.println( "Chrushed while getting points for training the tree.\nStack trace:" ); 
				e.printStackTrace(System.err);
			}		
			long timeEnd = System.currentTimeMillis();
			System.out.println("Time to get data set for training the tree: " + (((timeEnd - timeStart)/1000)/60) + "min");
			
			
			/*
			 * We will train it on the number of images that is needed in order to have 2 000 000 points.
			 */
			
			System.out.println(stars + "\nGetting the dimension bounds. ");
			double[] dimBounds = kdTreeInstance.getDimensionBounds(dataSetTraining);
			
			System.out.println(stars + "\nTraining the tree. ");
			try
			{
				timeStart = System.currentTimeMillis();
				kdTreeInstance.trainKDTree(dataSetTraining, 0, ((dataSetTraining.length / KeypointSingleton.DESCRIPTOR_LENGTH) - 1), trainedTree, dimBounds);
				timeEnd = System.currentTimeMillis();
				
				System.out.println("Training time for the tree: " +  (((timeEnd - timeStart)/1000)/60) + "min");						
	
				// write the trained tree to the disk under the name trainedTree.obj
				String trainedTreeName = "trainedTree.obj";
				File trainedTreeFile = new File(trainedTreeName);

				ObjectOutputStream out = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(trainedTreeFile)));
				out.writeObject(trainedTree);
				out.close();
				
			}catch(Exception e)
			{
				System.err.println(e.getMessage());
				System.err.println( "Crushed while training the tree and sving it to the disk.\nStack trace:" ); 
				e.printStackTrace(System.err);
			}
	      }else
	      {
	    	  System.out.println("in else statement");
	    	  
	    	  System.out.println(stars + "\nReading the trained tree from the disk. ");
	    	  String trainedTreeName = "trainedTree.obj";
			  File trainedTreeFile = new File(trainedTreeName);
			  try{
			  ObjectInputStream in = new ObjectInputStream(new
			  		BufferedInputStream(new FileInputStream(trainedTreeFile)));
			  trainedTree = (Vector<Node>)in.readObject();
			  
			  System.out.println("Read the trainedTree - its size is " + trainedTree.size());			  
			  
			  in.close();
			  }catch(IOException ioe)
			  {
				  System.err.println( "IOException stack trace:\n");
				  ioe.printStackTrace(System.err);
			  }catch(ClassNotFoundException cnfe)
			  {
				  System.err.println( "ClassNotFoundException stack trace:" ); 
				  cnfe.printStackTrace(System.err);
			  }
		}		
		
		kdTreeInstance.mapLeafsToBins(trainedTree);
	     	
		
		System.out.println(stars + "\nPopulating the tree. ");
		
		int dataSetPopulArrLength = IndexingTest1.NUM_IMAGES_POPULATING * 2000 * KeypointSingleton.DESCRIPTOR_LENGTH; 
		short[] dataSetPopulatingTemp = new short[dataSetPopulArrLength];
		short[] dataSetPopulating = new short[0];
		// imageIDs
		//int arraySize = dataSetPopulating.length / KeypointSingleton.DESCRIPTOR_LENGTH;
		int[] imageIDs = new int[IndexingTest1.NUM_IMAGES_POPULATING];
		// imageIndeces - used in populating the tree - each index coresponds to an image number and
		// the element at that position is the index of the first keypt of that image in the dataSetpopulating array
		int[] imageIndeces = new int[IndexingTest1.NUM_IMAGES_POPULATING];
		//maping of image names to IDs
		KeypointsExtraction keyptsExt = new KeypointsExtraction("helping constructor");
		Vector<String> nameIdMap = KeypointsExtraction.getNameIdMap();	
		
		/*
		 * Get the set of points for populating. We want NUM_IMAGES_POPULATING.
		 */
		if(populate)
		{
			System.out.println(stars + "\nCreate dataset for populating the tree. ");										
			
			long timeStart = System.currentTimeMillis();
			try
			{
				int position = 0;
				int image = 0;
				for(image = 0; image < IndexingTest1.NUM_IMAGES_POPULATING; image++)
				{					
					// imageIDs   	
					NumberFormat formatter = new DecimalFormat("00000");
					String objFileName = KeypointsExtraction.imagePath + "ukbench" + formatter.format(image) + ".obj";
					int imageID = nameIdMap.indexOf(objFileName);
					imageIDs[image] = imageID;
								
					
					short[] dataSet = KeypointsExtraction.getDataSet(image);
								
					System.arraycopy(dataSet, 0, dataSetPopulatingTemp, position, dataSet.length);
					imageIndeces[image] = position;
					position = position + dataSet.length;										
				}
				
				// copy dataSetPopulatingTemp to dataSetPopulating which is of the right size
				dataSetPopulating = new short[position];
				System.arraycopy(dataSetPopulatingTemp, 0, dataSetPopulating, 0, dataSetPopulating.length);
				dataSetPopulatingTemp = null;
				
				System.out.println("Number of images used for the dataSetPopulating is " + image);
				
			} catch(Exception e)
			{
				System.err.println( "Chrushed while getting points for populating the tree.\nStack trace:" ); 
				e.printStackTrace(System.err);
			}		
			long timeEnd = System.currentTimeMillis();
			System.out.println("Time to get data set for populating the tree: " + (((timeEnd - timeStart)/1000)/60) + "min");
			
			try
			{
				// write the data set for populating to the disk under the name opopulateDataSet.obj
				String populateDataSetName = "populateDataSet.obj";
				File populateDataSetFile = new File(populateDataSetName);
				ObjectOutputStream out = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(populateDataSetFile)));
				out.writeObject(dataSetPopulating);
				out.close();
				
	            // write the imageIndeces to the disk under the name imageIndeces.obj
				String imageIndecesName = "imageIndeces.obj";
				File imageIndecesFile = new File(imageIndecesName);	
				ObjectOutputStream out2 = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(imageIndecesFile)));
				out2.writeObject(imageIndeces);
				out2.close();
			}catch(IOException ioe)
		    {
			    System.err.println( "IOException stack trace:\n");
			    ioe.printStackTrace(System.err);
		    }
			
		}else
		{
			System.out.println(stars + "\nReading the data set for populating the tree from the disk. ");
	    	  String populateDataSetName = "populateDataSet.obj";			  
			  String imageIndecesName = "imageIndeces.obj";
			  
			  try{
				  ObjectInputStream in = new ObjectInputStream(new
				  		BufferedInputStream(new FileInputStream(populateDataSetName)));
				  dataSetPopulating = (short[])in.readObject();
				  
				  ObjectInputStream in2 = new ObjectInputStream(new
					  		BufferedInputStream(new FileInputStream(imageIndecesName)));			   
				  imageIndeces = (int[])in2.readObject();
				  
				  System.out.println("Read the data set for populating - its size is " + dataSetPopulating.length);
				  System.out.println("Read the imageIndeces used for populating - imageIndeces.size is " + imageIndeces.length);
				  
				  in.close();
				  }catch(IOException ioe)
				  {
					  System.err.println( "IOException stack trace:\n");
					  ioe.printStackTrace(System.err);
					  
			  }catch(ClassNotFoundException cnfe)
			  {
				  System.err.println( "ClassNotFoundException stack trace:" ); 
				  cnfe.printStackTrace(System.err);
			  }
			 
			  int image = 0;
			  for(image = 0; image < IndexingTest1.NUM_IMAGES_POPULATING; image++)
			  {					
					// imageIDs   	
					NumberFormat formatter = new DecimalFormat("00000");
					String objFileName = KeypointsExtraction.imagePath + "ukbench" + formatter.format(image) + ".obj";
					int imageID = nameIdMap.indexOf(objFileName);
					imageIDs[image] = imageID;																							
										
			  }
				System.out.println("Number of images used for the dataSetPopulating is " + image);

		}
		
		
		
		Vector<CountIndexPair> counts = new Vector<CountIndexPair>(IndexingTest1.NUM_LEAFS_IN_TRAINED_TREE);																
		
		PopulatedKDTree populatedTree = kdTreeInstance.populateKDTree(trainedTree, dataSetPopulating, imageIDs, imageIndeces, counts);
		
//		Do NUM_RUNS run, vote and record ranks...		
		double[] reciprocals = new double[IndexingTest1.NUM_RUNS * 3];
		int index = 0;
		int queryImage = 3;
		for(int j = 0; j < NUM_RUNS; j++)
		{					
			//take one image and query all the points in it and vote...			
			NumberFormat formatter = new DecimalFormat("00000");
			String imageFileName = KeypointsExtraction.imagePath + "ukbench" + formatter.format(queryImage) + ".obj";
		
		
			int[] votes = kdTreeInstance.voteForImage(trainedTree, populatedTree, counts, imageFileName, IndexingTest1.NUM_IMAGES_POPULATING); 			
			
				
			// we put all the image-votes pairs in a list and sort it
			LinkedList<ImageScorePair> scoresList = new LinkedList<ImageScorePair>();
			
			ImageScorePair pair1 = null, pair2 = null, pair3 = null, queryPair = null;
			for(int i = 0; i < votes.length; i++)
			{
				if(i != queryImage) // we take the query image out of the set
				{
					ImageScorePair pair = new ImageScorePair(i, votes[i]);
					scoresList.add(pair);
					
					if(i == queryImage - 3)
					{
						pair1 = pair;						
					}
					if(i == queryImage - 2)
					{						
						pair2 = pair;
					}
					if(i == queryImage - 1)
					{						
						pair3 = pair;
					}					
				}else
				{
					ImageScorePair pair = new ImageScorePair(i, votes[i]);
					
					System.out.println("queryPair has " + votes[i] + " scores and is located at position " + i);					
					
					scoresList.add(pair);
					queryPair = pair;
				}		
			}
			
			// create instance of Scores class with above created scoresList
			Scores scores = new Scores(scoresList);
			
			// sort scores - create instance of Scores class with the same scoresList, but sorted
			Scores scoresSorted = scores.sortScores(scores);
			// get the sorted list of scores		
			LinkedList<ImageScorePair> scoresListSorted = scoresSorted.getScoresList();
			
			// get and print the rank of queryPair - should be first
			int queryPairRank = scoresListSorted.indexOf(queryPair);
			System.out.println("\nRank of the queryPair - it should be 0: " + queryPairRank);
			
			//get and print out the number of votes for the queryPair - it should be equal to the number of keypoints in the query image
			int queryImageScore = queryPair.getScore();
			
			// check how many votes query image gets
			short[] dataSet = KeypointPairsExtraction.getDataSet(queryImage);
			System.out.println("The score of the queryImage: " + queryImageScore + 
					" and the number of keypoints in queryImage is: " + (dataSet.length / Keypoint.DESCRIPTOR_LENGTH));
						
			// now take the query image out of the scoresList so that it doesn't influence the results
			scoresListSorted.remove(queryPair);
			
			// Now for each of the remaining three images, we find the reciprocals and add them to the reciprocal set
			// For now, for simplicity we always query with the 4th of the 4 images in a block				
			
			int pair1Rank = scoresListSorted.indexOf(pair1) + 1;  // we don't need rank 0 bcuz that'd give us 1/0
			int pair2Rank = scoresListSorted.indexOf(pair2) + 1;
			int pair3Rank = scoresListSorted.indexOf(pair3) + 1;
			
			//sort the ranks
			int[] ranks = new int[3];
			ranks[0] = pair1Rank;
			ranks[1] = pair2Rank;
			ranks[2] = pair3Rank;
			Arrays.sort(ranks);
			pair1Rank = ranks[0];
			pair2Rank = ranks[1];
			pair3Rank = ranks[2];		
			
			double pair1Precision = 1.0/pair1Rank;
			double pair2Precision = 2.0/pair2Rank;
			double pair3Precision = 3.0/pair3Rank;
			
			System.out.println("pair1rank = " + pair1Rank + ", pair2rank = " + pair2Rank + ", pair3rank = " + pair3Rank);			
			
			reciprocals[index] = pair1Precision;
			index++;
			reciprocals[index] = pair2Precision;
			index++;
			reciprocals[index] = pair3Precision;
			index++;
			
			// take the next query image
			queryImage = queryImage + 4;
		}
		
		double sum1 = 0, sum2 = 0, sum3 = 0;
		
		for(int i = 0; i < reciprocals.length; i = i + 3)
		{
			sum1 = sum1 + reciprocals[i];
			sum2 = sum2 + reciprocals[i+1];
			sum3 = sum3 + reciprocals[i+2];
		}		
		
		double avrg1 = sum1 / IndexingTest1.NUM_RUNS;
		double avrg2 = sum2 / IndexingTest1.NUM_RUNS;
		double avrg3 = sum3 / IndexingTest1.NUM_RUNS;
		
		System.out.println("Average for the image 1: " + avrg1 + "\nAverage for the image 2: " + avrg2 + 
				"\nAverage for the image 3: " + avrg3);
		
		double avrg = (avrg1 + avrg2 + avrg3) / 3;
		
		System.out.println("\n\nOverall average is " + avrg + "\n");				
	}	
}
