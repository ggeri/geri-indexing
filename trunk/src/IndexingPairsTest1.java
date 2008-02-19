import java.io.*;
import java.text.*;
import java.util.*;


public class IndexingPairsTest1 {

/**
 * This attribute is used for training the tree. 
 * true = the tree will be trained 
 * false = the tree will be read from the disk - if there is a file with the specified name on the disk.
 */
private static boolean train = false;	
private static boolean populate = false;

//if there is apporximatelly 2000 points per image, then the number of bins is (num_images*2000)/bin_size
private static final int NUM_LEAFS_IN_TRAINED_TREE = (KeypointPairsExtraction.NUM_IMAGES * 2000)/KDTree.THRESHOLD; //1000000;
// tree has three times more nodes then leafs
private static final int NUM_NODES_IN_TRAINED_TREE = 3 * IndexingPairsTest1.NUM_LEAFS_IN_TRAINED_TREE;
private static final int NUM_PTS_TRAINING = 10000000;	
private static final int NUM_IMAGES_TRAINING = 10200;
private static final int NUM_IMAGES_POPULATING = 10200; 	// set heap to 6144max & w/ flag -XX:-UseGCOverheadLimit
private static final int NUM_RUNS = 2549;
	
	public static void main(String[] args) 
	{	
		String stars = "\n**************************************************************************************\n";
				
		
		/* Instantiate the KeypointsExtraction class - this converts jpg files to pgm, extracts the keypoints and saves them
		 * to .key files and finally reads thesefiles and creates serialized javaobjects on the disk.
		 */
		System.out.println(stars + "\nCreate .pair.obj file for each image and store them to the disk. ");
		//KeypointPairsExtraction keyptPairExtract = new KeypointPairsExtraction();


		//Get an instance of the KDTree class
		KDTree kdTreeInstance = new KDTree();
		Vector<Node> trainedTree = new Vector<Node>(IndexingPairsTest1.NUM_NODES_IN_TRAINED_TREE); 				
		
		int dataSetArrLength = IndexingPairsTest1.NUM_PTS_TRAINING * Keypoint.DESCRIPTOR_LENGTH; 
		
		
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
				for(image = 0; image < IndexingPairsTest1.NUM_IMAGES_TRAINING && position < dataSetArrLength; image++)
				{					
					short[] dataSet = KeypointPairsExtraction.getDataSet(image);					
					
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
				System.out.println("Number of images used for the dataSetTraining is " + image);
				System.out.println("Number of points used for the dataSetTraining is " + dataSetTraining.length);
				
			} catch(Exception e)
			{
				System.err.println( "Chrushed while getting points for training the tree.\nStack trace:" ); 
				e.printStackTrace(System.err);
			}		
			long timeEnd = System.currentTimeMillis();
			System.out.println("Time to get data set for training the tree: " + (((timeEnd - timeStart)/1000)/60) + "min");				
			
			System.out.println(stars + "\nGetting the dimension bounds. ");
			double[] dimBounds = kdTreeInstance.getDimensionBounds(dataSetTraining);
			
			System.out.println(stars + "\nTraining the tree. ");
			try
			{
				timeStart = System.currentTimeMillis();
				kdTreeInstance.trainKDTree(dataSetTraining, 0, ((dataSetTraining.length / KeypointSingleton.DESCRIPTOR_LENGTH) - 1), trainedTree, dimBounds);

				//trim the trainedTree to it's actual size
				trainedTree.trimToSize();
				
				timeEnd = System.currentTimeMillis();
				
				System.out.println("Training time for the tree: " +  (((timeEnd - timeStart)/1000)/60) + "min");
				
	
				// write the trained tree to the disk under the name trainedTree.obj
				// changed from trainedTreePairs.obj to trainedTreePairs2.obj in order to create save a new trained tree 
				String trainedTreeName = "trainedTreePairs.obj"; 
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
	    	  String trainedTreeName = "trainedTreePairs.obj";
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
	     			
		
		/*
		 * 
		 * 
		 */
		System.out.println(stars + "\nPopulating the tree. ");
		
		int dataSetPopulArrLength = IndexingPairsTest1.NUM_IMAGES_POPULATING * 2000 * KeypointSingleton.DESCRIPTOR_LENGTH; 
		short[] dataSetPopulatingTemp = new short[dataSetPopulArrLength];
		short[] dataSetPopulating = new short[0];
		// imageIDs
		//int arraySize = dataSetPopulating.length / KeypointSingleton.DESCRIPTOR_LENGTH;
		int[] imageIDs = new int[IndexingPairsTest1.NUM_IMAGES_POPULATING];
		// imageIndeces - used in populating the tree - each index coresponds to an image number and
		// the element at that position is the index of the first keypt of that image in the dataSetpopulating array
		int[] imageIndeces = new int[IndexingPairsTest1.NUM_IMAGES_POPULATING];
		//maping of image names to IDs
		KeypointPairsExtraction keyptsExt = new KeypointPairsExtraction("helping constructor");
		Vector<String> nameIdMap = KeypointPairsExtraction.getNameIdMap();	
		
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
				for(image = 0; image < IndexingPairsTest1.NUM_IMAGES_POPULATING; image++)
				{					
					// imageIDs   	
					NumberFormat formatter = new DecimalFormat("00000");
					String objFileName = KeypointsExtraction.imagePath + "ukbench" + formatter.format(image) + ".pair.obj";
					int imageID = nameIdMap.indexOf(objFileName);
					imageIDs[image] = imageID;
					
					short[] dataSet = KeypointPairsExtraction.getDataSet(image);									
								
					System.arraycopy(dataSet, 0, dataSetPopulatingTemp, position, dataSet.length);
					imageIndeces[image] = position;
					position = position + dataSet.length;														
				}
				
				// copy dataSetPopulatingTemp to dataSetPopulating which is of the right size
				dataSetPopulating = new short[position];
				System.arraycopy(dataSetPopulatingTemp, 0, dataSetPopulating, 0, dataSetPopulating.length);
				dataSetPopulatingTemp = null;
				
				System.out.println("Number of keypoints in dataSetPopulating " + dataSetPopulating.length/Keypoint.DESCRIPTOR_LENGTH);
				
			} catch(Exception e)
			{				 
				e.printStackTrace(System.err);
			}		
			long timeEnd = System.currentTimeMillis();
			System.out.println("Time to get data set for populating the tree: " + (((timeEnd - timeStart)/1000)/60) + "min");
			
			try
			{
				// write the data set for populating to the disk under the name populateDataSet.obj
				String populateDataSetName = "populateDataSetPairs.obj";
				File populateDataSetFile = new File(populateDataSetName);
				ObjectOutputStream out = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(populateDataSetFile)));
				out.writeObject(dataSetPopulating);
				out.close();
				
	            // write the imageIndeces to the disk under the name imageIndeces.obj
				String imageIndecesName = "imageIndecesPairs.obj";
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
	    	  String populateDataSetName = "populateDataSetPairs.obj";			  
			  String imageIndecesName = "imageIndecesPairs.obj";
			  
			  try{
				  ObjectInputStream in = new ObjectInputStream(new
				  		BufferedInputStream(new FileInputStream(populateDataSetName)));
				  dataSetPopulating = (short[])in.readObject();
				  
				  ObjectInputStream in2 = new ObjectInputStream(new
					  		BufferedInputStream(new FileInputStream(imageIndecesName)));			   
				  imageIndeces = (int[])in2.readObject();
				  
				  System.out.println("Read dataSetPopulating for populating - its size is " + dataSetPopulating.length);
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
			  for(image = 0; image < IndexingPairsTest1.NUM_IMAGES_POPULATING; image++)
			  {					
					// imageIDs   	
					NumberFormat formatter = new DecimalFormat("00000");
					String objFileName = KeypointsExtraction.imagePath + "ukbench" + formatter.format(image) + ".pair.obj";
					// imageID is a unique id of each image that coresponds to it name
					int imageID = nameIdMap.indexOf(objFileName);
					imageIDs[image] = imageID;																							
													
			  }			  			  
		}								
		
		// contains a CountIndexPair per image for all the images in populated tree - count gives number of points 
		//in that specific image and index gives the index of the first point of that image
		Vector<CountIndexPair> counts = new Vector<CountIndexPair>(IndexingPairsTest1.NUM_LEAFS_IN_TRAINED_TREE);																
		
		// imageIndeces tells us where each new image starts in this lond dataSetPopulating array
		PopulatedKDTree populatedTree = kdTreeInstance.populateKDTree(trainedTree, dataSetPopulating, imageIDs, imageIndeces, counts);																			
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		
        // Do NUM_RUNS run, vote and record ranks...
		double[] reciprocals = new double[IndexingPairsTest1.NUM_RUNS * 3];
		int index = 0;
		int queryImage = 3; 
		
		// To get average number of votes for the 3 matching images and the average for the 
		// remaining (NUM_IMAGES - 4) images
		int im1VotesSum = 0, im2VotesSum = 0, im3VotesSum = 0;
		int imRestSum = 0;
		
		for(int j = 0; j < NUM_RUNS; j++)
		{			
			//take one image and query with all the points in it and vote...
			
			NumberFormat formatter = new DecimalFormat("00000");
			String imageFileName = KeypointPairsExtraction.imagePath + "ukbench" + formatter.format(queryImage) + ".q.pair.obj";														
		
			int[] votes = kdTreeInstance.voteForImage(trainedTree, populatedTree, counts, imageFileName, IndexingPairsTest1.NUM_IMAGES_POPULATING); 
				
//			System.out.println("\nImage " + queryImage + " has " + votes[queryImage] + " votes.");
//			System.out.println("In votes array, queryImage is located at position " + queryImage + " and has " + votes[queryImage] + " votes.");
			
			// we put all the image-votes pairs in a list and sort it
			LinkedList<ImageScorePair> scoresList = new LinkedList<ImageScorePair>();
			
			ImageScorePair pair1 = null, pair2 = null, pair3 = null, queryPair = null;
			for(int i = 0; i < votes.length; i++)
			{
				// add change here
				if(i != queryImage) // we take the query image out of the set
				{
					ImageScorePair pair = new ImageScorePair(i, votes[i]);
					scoresList.add(pair);
					
					if(i == queryImage - 3)
					{
						pair1 = pair;	
						im1VotesSum = im1VotesSum + pair.getScore();
					}
					else if(i == queryImage - 2)
					{						
						pair2 = pair;
						im2VotesSum = im2VotesSum + pair.getScore();
					}
					else if(i == queryImage - 1)
					{						
						pair3 = pair;
						im3VotesSum = im3VotesSum + pair.getScore();
					}	
					else
					{
						imRestSum = imRestSum + pair.getScore();
					}
				}else
				{
					ImageScorePair pair = new ImageScorePair(i, votes[i]);									
					
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
//			int queryPairRank = scoresListSorted.indexOf(queryPair);
//			System.out.println("Rank of the queryImage in array of sorted scores - it should be 0: " + queryPairRank);
			
			//get and print out the number of votes for the queryPair - it should be equal to the number of keypoints in the query image
//			int queryImageScore = queryPair.getScore();
//			
//			// check how many votes query image gets
//			short[] dataSet = KeypointPairsExtraction.getDataSet(queryImage);
//			System.out.println("queryImage score & number of pairs in queryImage should be the same: " + 
//					queryImageScore + " = " + (dataSet.length / Keypoint.DESCRIPTOR_LENGTH));
			
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
						
			//System.out.println("pair1rank = " + pair1Rank + "  pair2rank = " + pair2Rank + "  pair3rank = " + pair3Rank);
			
			reciprocals[index] = pair1Precision;
			index++;
			reciprocals[index] = pair2Precision;
			index++;
			reciprocals[index] = pair3Precision;
			index++;
						
			// take the next queryImage
			queryImage = queryImage + 4;
		}
		
		// Average number of votes for matching images 1, 2 and 3
		int im1AvrgVotes = im1VotesSum / IndexingPairsTest1.NUM_RUNS;
		int im2AvrgVotes = im2VotesSum / IndexingPairsTest1.NUM_RUNS;
		int im3AvrgVotes = im3VotesSum / IndexingPairsTest1.NUM_RUNS;
		System.out.println("\nAvrg votes im1 = " + im1AvrgVotes + ", avrg votes im2 = " + im2AvrgVotes + ", avrg votes im3 = " + im3AvrgVotes);
		
		// Average for all 3 mathing images
		int imMatchAvrgVotes = (im1AvrgVotes + im2AvrgVotes + im3AvrgVotes) / 3;
		System.out.println("Average votes for matching im (im1,im2,im3) = " + imMatchAvrgVotes);
		
		// Average number of votes for all other non matching images
		int imRestAvrgVotes = imRestSum / (IndexingPairsTest1.NUM_RUNS * (KeypointsExtraction.NUM_IMAGES-4));
		System.out.println("Average votes for non-matching im (all except im1,im2,im3,im4) = " + imRestAvrgVotes);
		
		double sum1 = 0, sum2 = 0, sum3 = 0;
		
		for(int i = 0; i < reciprocals.length; i = i + 3)
		{
			sum1 = sum1 + reciprocals[i];
			sum2 = sum2 + reciprocals[i+1];
			sum3 = sum3 + reciprocals[i+2];
		}
		
		double avrg1 = sum1 / IndexingPairsTest1.NUM_RUNS;
		double avrg2 = sum2 / IndexingPairsTest1.NUM_RUNS;
		double avrg3 = sum3 / IndexingPairsTest1.NUM_RUNS;
		
		System.out.println("\nAverage prec for the image 1: " + avrg1 + "\nAverage prec for the image 2: " + avrg2 + 
				"\nAverage prec for the image 3: " + avrg3);
		
		double avrg = (avrg1 + avrg2 + avrg3) / 3;
		
		System.out.println("\n\nOverall average precisio is " + avrg + "\n");	
	}	
}
