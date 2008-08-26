package ca.utoronto.cs.kdtree.groundtruth;
import java.io.*;
import java.text.*;
import java.util.*;

import ca.utoronto.cs.kdtree.basic.CountIndexPair;
import ca.utoronto.cs.kdtree.basic.ImageScorePair;
import ca.utoronto.cs.kdtree.basic.KDTree;
import ca.utoronto.cs.kdtree.basic.Keypoint;
import ca.utoronto.cs.kdtree.basic.KeypointSingleton;
import ca.utoronto.cs.kdtree.basic.KeypointsExtraction;
import ca.utoronto.cs.kdtree.basic.Node;
import ca.utoronto.cs.kdtree.basic.Path;
import ca.utoronto.cs.kdtree.basic.PopulatedKDTree;
import ca.utoronto.cs.kdtree.basic.Scores;


public class IndexingGroundTruthTest2 {

/**
 * This attribute is used for training the tree. 
 * true = the tree will be trained 
 * false = the tree will be read from the disk - if there is a file with the specified name on the disk.
 */
private static boolean train = true;	
private static boolean populate = true;

// if there is apporximatelly 2000 points per image, then the number of bins is (num_images*2000)/bin_size
private static final int NUM_LEAFS_IN_TRAINED_TREE = (KeypointsExtraction.NUM_IMAGES * 2000)/KDTree.THRESHOLD; //1000000;
// tree has three times more nodes then leafs
private static final int NUM_NODES_IN_TRAINED_TREE = 3 * IndexingGroundTruthTest2.NUM_LEAFS_IN_TRAINED_TREE;
private static final int NUM_PTS_TRAINING = 10000000;	
private static final int NUM_IMAGES_TRAINING = 10200;
private static final int NUM_IMAGES_POPULATING = 10200; 	// works with 1000 if heap set to 6144max & w/ flag -XX:-UseGCOverheadLimit
private static final int NUM_RUNS = 2550;
	
	public static void main(String[] args)
	{	
		long beginTime = System.currentTimeMillis();
		
		String stars = "\n**************************************************************************************\n";
				
		
		/* Instantiate the KeypointsExtraction class - this converts jpg files to pgm, extracts the keypoints and saves them
		 * to .key files and finally reads thesefiles and creates serialized javaobjects on the disk.
		 */
		//System.out.println(stars + "\nCreate .obj file for each image and store them to the disk. ");
		//KeypointsExtraction keyptExtract = new KeypointsExtraction();
		
		//System.out.println("Done creating obj files");
		//System.exit(1);
		
		
		//Get an instance of the KDTree class		 
		KDTree kdTreeInstance = new KDTree();
		Vector<Node> trainedTree = new Vector<Node>(IndexingGroundTruthTest2.NUM_NODES_IN_TRAINED_TREE); 
		
		int dataSetArrLength = IndexingGroundTruthTest2.NUM_PTS_TRAINING * Keypoint.DESCRIPTOR_LENGTH;  
		
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
				for(image = 0; image < IndexingGroundTruthTest2.NUM_IMAGES_TRAINING && position < dataSetArrLength; image++)
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
				System.out.println("Number of points used for the dataSetTreining is " + dataSetTraining.length);
				
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
	
				// write the trained tree to the disk under the name trainedTreeSingl.obj
				String trainedTreeName = "trainedTreeSinglFullDim.obj";
				//String trainedTreeName = "trainedTreeSingl.obj";
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
	    	  String trainedTreeName = "trainedTreeSinglFullDim.obj";	
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
		
		int dataSetPopulArrLength = IndexingGroundTruthTest2.NUM_IMAGES_POPULATING * 2000 * KeypointSingleton.DESCRIPTOR_LENGTH; 
		short[] dataSetPopulatingTemp = new short[dataSetPopulArrLength];
		short[] dataSetPopulating = new short[0];
		// imageIDs
		//int arraySize = dataSetPopulating.length / KeypointSingleton.DESCRIPTOR_LENGTH;
		int[] imageIDs = new int[IndexingGroundTruthTest2.NUM_IMAGES_POPULATING];
		// imageIndeces - used in populating the tree - each index coresponds to an image number and
		// the element at that position is the index of the first keypt of that image in the dataSetpopulating array
		int[] imageIndeces = new int[IndexingGroundTruthTest2.NUM_IMAGES_POPULATING];
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
				for(image = 0; image < IndexingGroundTruthTest2.NUM_IMAGES_POPULATING; image++)
				{					
					// imageIDs   	
					//NumberFormat formatter = new DecimalFormat("00000");
					//String objFileName = KeypointsExtraction.imagePath + "ukbench" + formatter.format(image) + ".obj";
					//int imageID = nameIdMap.indexOf(objFileName);
					//imageIDs[image] = imageID;
					imageIDs[image] = image;
								
					
					short[] dataSet = KeypointsExtraction.getDataSet(image);
								
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
				System.err.println( "Chrushed while getting points for populating the tree.\nStack trace:" ); 
				e.printStackTrace(System.err);
			}		
			long timeEnd = System.currentTimeMillis();
			System.out.println("Time to get data set for populating the tree: " + (((timeEnd - timeStart)/1000)/60) + "min");
			
			try
			{
				// write the data set for populating to the disk under the name opopulateDataSet.obj
				String populateDataSetName = "dataForTreePopulatingSinglFullDim.obj";
				//String populateDataSetName = "dataForTreePopulatingSingl.obj";
				File populateDataSetFile = new File(populateDataSetName);
				ObjectOutputStream out = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(populateDataSetFile)));
				out.writeObject(dataSetPopulating);
				out.close();
				
	            // write the imageIndeces to the disk under the name imageIndeces.obj
				String imageIndecesName = "imageIndecesSinglFullDim.obj";
				//String imageIndecesName = "imageIndecesSingl.obj";
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
	    	  String populateDataSetName = "dataForTreePopulatingSinglFullDim.obj";				  
			  String imageIndecesName = "imageIndecesSinglFullDim.obj";	
			  
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
			  for(image = 0; image < IndexingGroundTruthTest2.NUM_IMAGES_POPULATING; image++)
			  {					
					// imageIDs   	
					//NumberFormat formatter = new DecimalFormat("00000");
					//String objFileName = KeypointsExtraction.imagePath + "ukbench" + formatter.format(image) + ".obj";
					//int imageID = nameIdMap.indexOf(objFileName);
					//imageIDs[image] = imageID;																							
				  imageIDs[image] = image;
										
			  }
				System.out.println("Number of images used for the dataForTreePopulating is " + image);

		}
		
		
		// contains a CountIndexPair per image for all the images in populated tree - count gives number of points 
		//in that specific image and index gives the index of the first point of that image
		Vector<CountIndexPair> counts = new Vector<CountIndexPair>(IndexingGroundTruthTest2.NUM_LEAFS_IN_TRAINED_TREE);																
		
		PopulatedKDTree populatedTree = kdTreeInstance.populateKDTree(trainedTree, dataSetPopulating, imageIDs, imageIndeces, counts);
		
//		Do NUM_RUNS run, vote and record ranks...		
		double[] reciprocals = new double[IndexingGroundTruthTest2.NUM_RUNS * 3];
		int index = 0;
		int queryImage = 3;
		
		// To get average number of votes for the 3 matching images and the average for the 
		// remaining (NUM_IMAGES - 4) images
		double im1VotesSum = 0, im2VotesSum = 0, im3VotesSum = 0;
		double imRestSum = 0;
		
		// stats is array of double of length 13 - it contains 4 averages for varianceWide, 4 averages for 
		// varianceNarrow and 4 averages for varianceNarrow where values are normalized by densityEstimate
		// and the last value is the number of points from this image that voted
		// This information gathered in vote-for-point function
		double[] stats = {0.0, 0.0, 0.0, 0.0,	// average for bins 1-4 for varianceWide
						  0.0, 0.0, 0.0, 0.0,	// average for bins 1-4 for varianceNarrow
						  0.0, 0.0, 0.0, 0.0,	// average for bins 1-4 for varianceNarrow normalized
						  0.0 };			// number of points that voted
		
		for(int j = 0; j < NUM_RUNS; j++)
		{					
			//take one image and query all the points in it and vote...			
			NumberFormat formatter = new DecimalFormat("00000");			
			// get the path names from the Paths class
			//String imageFileName = Path.imageSinglObj + "ukbench" + formatter.format(queryImage) + ".obj";
			
			String imageFileName = Path.fileSinglOBJ_4 + "ukbench" + formatter.format(queryImage) + ".obj";
			//String imageFileName = Path.fileSinglOBJ_3 + "ukbench" + formatter.format(queryImage) + ".obj";
			
			double[] votes = kdTreeInstance.voteForImageWithWeightsAndDensityWithStats(trainedTree, 
					populatedTree, counts, imageFileName, IndexingGroundTruthTest2.NUM_IMAGES_POPULATING, stats); 									
			
//			System.out.println("\nImage " + queryImage + " has " + votes[queryImage] + " votes.");
//			System.out.println("In votes array, queryImage is located at position " + queryImage + " and has " + votes[queryImage] + " votes.");
				
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
			//Scores scoresSorted = scores.sortScores(scores);
			// get the sorted list of scores		
			//LinkedList<ImageScorePair> scoresListSorted = scoresSorted.getScoresList();
			
			// sort the scoresList
			scores.sortDescending();
			
			// get the sorted scoreList
			LinkedList<ImageScorePair> scoresListSorted = scores.getScoreList();
			
			// get and print the rank of queryPair - should be first
			//int queryPairRank = scoresListSorted.indexOf(queryPair);
			//System.out.println("Rank of the queryImage in array of sorted scores - it should be 0: " + queryPairRank);
			
			//get and print out the number of votes for the queryPair - it should be equal to the number of keypoints in the query image
			//int queryImageScore = queryPair.getScore();
			
			// check how many votes query image gets
			//short[] dataSet = KeypointsExtraction.getDataSet(queryImage);
			//System.out.println("queryImage score & number of keypoints in queryImage should be the same: " + 
			//		queryImageScore + " = " + (dataSet.length / Keypoint.DESCRIPTOR_LENGTH));
						
			// now take the query image out of the scoresList so that it doesn't influence the results
			scoresListSorted.remove(queryPair);
			
/////////////////////////////////			
			
			// Now for each of the remaining three images, we find the reciprocals of their ranks
			// and add them to the reciprocal set
			// For now, for simplicity we always query with the 4th of the 4 images in a block				
			
			int pair1Rank = scoresListSorted.indexOf(pair1) + 1;  // we don't need rank 0 bcuz that'd give us 1/0
			int pair2Rank = scoresListSorted.indexOf(pair2) + 1;
			int pair3Rank = scoresListSorted.indexOf(pair3) + 1;
			
			// Sort the ranks
			int[] ranks = new int[3];
			ranks[0] = pair1Rank;	// rank of first image in the set of 4
			ranks[1] = pair2Rank;	// rank of second image in the set of 4
			ranks[2] = pair3Rank;	// rank of third image in the set of 4
									// fourth image in the set of 4 is the query image
			Arrays.sort(ranks);
			int bestRank = ranks[0];	// best rank of the three
			int midRank = ranks[1];		// mid rank of the three
			int worstRank = ranks[2];	// worst rank of the three
			
			double precision1 = 1.0/bestRank;
			double precision2 = 2.0/midRank;
			double precision3 = 3.0/worstRank;
						
			//System.out.println("Best rank = " + bestRank + ",\tmidle rank = " + midRank + ",\tworstRank3 = " + worstRank);			
			
			// index i - best scored image of 3 (4-query) images
			// index i+1 - mid scored image of 3 (4-query) images
			// index i+2 - worst scored image of 3 (4-query) images
			reciprocals[index] = precision1;
			index++;
			reciprocals[index] = precision2;
			index++;
			reciprocals[index] = precision3;
			index++;
			
			// take the next query image
			queryImage = queryImage + 4;
		}
		
		// Print out the statistics for varianceWide
		System.out.println("############################################################");
		System.out.println("Statistics on average weight for bins 1-4, for varianceWide");
		System.out.println("Bin 1:\t average = " + stats[0] / stats[12]);
		System.out.println("Bin 2:\t average = " + stats[1] / stats[12]);
		System.out.println("Bin 3:\t average = " + stats[2] / stats[12]);
		System.out.println("Bin 200:\t average = " + stats[3] / stats[12]);
		System.out.println("############################################################");
		System.out.println("Statistics on average weight for bins 1-4, for varianceNarrow");
		System.out.println("Bin 1:\t average = " + stats[4] / stats[12]);
		System.out.println("Bin 2:\t average = " + stats[5] / stats[12]);
		System.out.println("Bin 3:\t average = " + stats[6] / stats[12]);
		System.out.println("Bin 200:\t average = " + stats[7] / stats[12]);
		System.out.println("############################################################");
		System.out.println("Statistics on average weight for bins 1-4, for varianceNarrow with normalization");
		System.out.println("Bin 1:\t average = " + stats[8] / stats[12]);
		System.out.println("Bin 2:\t average = " + stats[9] / stats[12]);
		System.out.println("Bin 3:\t average = " + stats[10] / stats[12]);
		System.out.println("Bin 200:\t average = " + stats[11] / stats[12]);
		System.out.println("############################################################");
		
		//System.out.println("Last image used is image " + (queryImage - 4));
		
		// Average number of votes for matching images 1, 2 and 3 (not sorted)
		double im1AvrgVotes = im1VotesSum / IndexingGroundTruthTest2.NUM_RUNS;
		double im2AvrgVotes = im2VotesSum / IndexingGroundTruthTest2.NUM_RUNS;
		double im3AvrgVotes = im3VotesSum / IndexingGroundTruthTest2.NUM_RUNS;
		System.out.println("\nAverage number of votes for matching images - not sorted by number of votes:");
		System.out.println("\nAvrg votes im1 = " + im1AvrgVotes + ", avrg votes im2 = " + im2AvrgVotes + ", avrg votes im3 = " + im3AvrgVotes);
		
		// Average for all 3 matching images
		double imMatchAvrgVotes = (im1AvrgVotes + im2AvrgVotes + im3AvrgVotes) / 3;
		System.out.println("Average votes for all 3 matching imgs (im1,im2,im3) = " + imMatchAvrgVotes);
		
		// Average number of votes for all other non matching images
		double imRestAvrgVotes = imRestSum / (IndexingGroundTruthTest2.NUM_RUNS * (KeypointsExtraction.NUM_IMAGES-4));
		System.out.println("Average votes for all non-matching imgs (all except im1,im2,im3,im4) = " + imRestAvrgVotes);
		
		// sum1 - sum of reciprocals for best scored images
		// sum2 - sum of reciprocals for mid scored images
		// sum3 - sum of reciprocals for worst scored images
		double sum1 = 0, sum2 = 0, sum3 = 0;
		
		for(int i = 0; i < reciprocals.length; i = i + 3)
		{
			sum1 = sum1 + reciprocals[i];
			sum2 = sum2 + reciprocals[i+1];
			sum3 = sum3 + reciprocals[i+2];
		}		
		
		double avrg1 = sum1 / IndexingGroundTruthTest2.NUM_RUNS;
		double avrg2 = sum2 / IndexingGroundTruthTest2.NUM_RUNS;
		double avrg3 = sum3 / IndexingGroundTruthTest2.NUM_RUNS;
		
		System.out.println("\nAverage prec for the best scored image: " + avrg1 + "\nAverage prec for the mid scored image: " + avrg2 + 
				"\nAverage prec for the worst scored image: " + avrg3);
		
		double avrg = (avrg1 + avrg2 + avrg3) / 3;
		
		System.out.println("\n\nOverall average precision is " + avrg + "\n");		
		
		long endTime = System.currentTimeMillis();
		System.out.println("\nTime needed to index with all images: " + (((endTime - beginTime)/1000)/60) + "min");
		
	}	
}
