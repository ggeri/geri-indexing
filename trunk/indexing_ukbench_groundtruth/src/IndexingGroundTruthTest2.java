import java.io.*;
import java.text.*;
import java.util.*;

public class IndexingGroundTruthTest2 {

/**
 * This attribute is used for training the tree. 
 * true = the tree will be trained 
 * false = the tree will be read from the disk - if there is a file with the specified name on the disk.
 */
private static boolean train = false;	
private static boolean populate = false;

// if there is apporximatelly 2000 points per image, then the number of bins is (num_images*2000)/bin_size
public static final int NUM_LEAFS_IN_TRAINED_TREE = (KeypointsExtraction.NUM_IMAGES * 2000)/KDTree.THRESHOLD; //1000000;
// tree has three times more nodes then leafs
private static final int NUM_NODES_IN_TRAINED_TREE = 3 * IndexingGroundTruthTest2.NUM_LEAFS_IN_TRAINED_TREE;
private static final int NUM_PTS_TRAINING = 10000000;	
private static final int NUM_IMAGES_TRAINING = 10200;
private static final int NUM_IMAGES_POPULATING = 10200; 	// works with 1000 if heap set to 6144max & w/ flag -XX:-UseGCOverheadLimit
private static final int NUM_RUNS = 100;
	
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
				//String trainedTreeName = "trainedTreeSinglFullDim.obj";
				//String trainedTreeName = "trainedTreeSingl.obj";
				String trainedTreeName = "trainedTreeSinglNormalized.obj";
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
	    	  //System.out.println("in else statement");
	    	  
	    	  System.out.println(stars + "\nReading the trained tree from the disk. ");
	    	  //String trainedTreeName = "trainedTreeSinglFullDim.obj";
	    	 // String trainedTreeName = "trainedTreeSingl.obj";
	    	  String trainedTreeName = "trainedTreeSinglNormalized.obj";
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
	     	
		
		System.out.println(stars + "\nGetting the data for populating the tree and populating the tree. ");
		
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
			System.out.println(stars + "\nCreate dataset for populating the tree. \n");										
			
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
				//String populateDataSetName = "dataForTreePopulatingSinglFullDim.obj";
				//String populateDataSetName = "dataForTreePopulatingSingl.obj";
				String populateDataSetName = "dataForTreePopulatingSinglNormalized.obj";
				File populateDataSetFile = new File(populateDataSetName);
				ObjectOutputStream out = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(populateDataSetFile)));
				out.writeObject(dataSetPopulating);
				out.close();
				
	            // write the imageIndeces to the disk under the name imageIndeces.obj
				//String imageIndecesName = "imageIndecesSinglFullDim.obj";
				//String imageIndecesName = "imageIndecesSingl.obj";
				String imageIndecesName = "imageIndecesSinglNormalized.obj";
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
			System.out.println(stars + "\nReading the data set for populating the tree from the disk. \n");
	    	  //String populateDataSetName = "dataForTreePopulatingSinglFullDim.obj";				  
			  //String imageIndecesName = "imageIndecesSinglFullDim.obj";	
			  //String populateDataSetName = "dataForTreePopulatingSingl.obj";
			  //String imageIndecesName = "imageIndecesSingl.obj";
			  String populateDataSetName = "dataForTreePopulatingSinglNormalized.obj";
			  String imageIndecesName = "imageIndecesSinglNormalized.obj";
			  
			  try{
				  ObjectInputStream in = new ObjectInputStream(new
				  		BufferedInputStream(new FileInputStream(populateDataSetName)));
				  dataSetPopulating = (short[])in.readObject();
				  //System.out.println("Total number of points in the tree: " + dataSetPopulating.length / Keypoint.DESCRIPTOR_LENGTH);
				  
				  ObjectInputStream in2 = new ObjectInputStream(new
					  		BufferedInputStream(new FileInputStream(imageIndecesName)));			   
				  imageIndeces = (int[])in2.readObject();
				  
				  System.out.println("Finished reading dataSetPopulating for populating - its size is " + dataSetPopulating.length);
				  System.out.println("Finished reading imageIndeces used for populating - imageIndeces size is " + imageIndeces.length);
				  
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
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		System.out.println(stars + "\nPopulating the tree with the read data.\n");
		
		// contains a CountIndexPair per image for all the images in populated tree - count gives number of points 
		//in that specific image and index gives the index of the first point of that image
		Vector<CountIndexPair> counts = new Vector<CountIndexPair>(IndexingGroundTruthTest2.NUM_LEAFS_IN_TRAINED_TREE);																
		
		PopulatedKDTree populatedTree = kdTreeInstance.populateKDTree(trainedTree, dataSetPopulating, imageIDs, imageIndeces, counts);
		
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		System.out.println(stars + "\nVoting and writing out the analysis information.\n");
		
		
		// Vote with one image and write out a text file named votingInfoIm0.txt (ex. for im0) with the following info:
		// what point in order is this point, in this image - one number, alone on line
		// info about retrieved bins: bin_num, point_num, im_num - bin_nums on one line, point_nums on next line, im_nums on third line etc.
		// ex. 5 - point that is voting
		// 20 - number of bins that the point that is voting retrieves
		// 300 - number of points that the point that is voting retrieves
		// 3 9 8 9 7 - numbers of the retrieved bins
		// 980 456 342 980 459 - point numbers in the below image
		// 45 9 32 239 9987 - image numbers
		// 23.4 55.6 4.5 67.5 34.2 - distnaces to nearest edge of a bin (for each retrieved pt)
		
		
		//double[] reciprocals = new double[IndexingGroundTruthTest2.NUM_RUNS * 3];
		int index = 0;
		
		// Read the file with query images and queyr points; file name = randomPoints.txt
		// and store read image and point IDs to two arrays
		int[] queryImages = new int[100];
		int[] queryPts = new int[100];
		
		String fileName = Path.randomQueryPtsFileSingl;
		BufferedReader inputStream = null;		
		
		try {
			inputStream = 
                new BufferedReader(new FileReader(fileName));
			String imagesLine = inputStream.readLine();
			String pointsLine = inputStream.readLine();
			
			String[] imagesTokens = imagesLine.split(" ");
			String[] pointsTokens = pointsLine.split(" ");

			for(int i = 0; i < imagesTokens.length; i++)
			{
				int imageID = Integer.parseInt(imagesTokens[i]); 
				queryImages[i] = imageID;
			}
			
			for(int i = 0; i < pointsTokens.length; i++)
			{
				int pointID = Integer.parseInt(pointsTokens[i]); 
				queryPts[i] = pointID - 1; // -1 becasue the info thta comes from matlab starts at 1, not 0
			}
			
		} catch (FileNotFoundException e1) 
		{
			e1.printStackTrace();
		} catch (IOException e) 
		{
			e.printStackTrace();
		}
		
		try {
			inputStream.close();
		} catch (IOException e1) 
		{
			e1.printStackTrace();
		}
		
//		for(int g = 0; g < queryImages.length; g++)
//		{
//			System.out.print(queryImages[g] + " ");
//		}
//		System.out.println();
//		for(int g = 0; g < queryPts.length; g++)
//		{
//			System.out.print(queryPts[g] + " ");
//		}
//		System.exit(0);
		
		// To get average number of votes for the 3 matching images and the average for the 
		// remaining (NUM_IMAGES - 4) images
		//double im1VotesSum = 0, im2VotesSum = 0, im3VotesSum = 0;
		//double imRestSum = 0;
		
		// stats is array of double of length 13 - it contains 4 averages for varianceWide, 4 averages for 
		// varianceNarrow and 4 averages for varianceNarrow where values are normalized by densityEstimate
		// and the last value is the number of points from this image that voted
		// This information gathered in vote-for-point function
		/*double[] stats = {0.0, 0.0, 0.0, 0.0,	// average for bins 1-4 for varianceWide
						  0.0, 0.0, 0.0, 0.0,	// average for bins 1-4 for varianceNarrow
						  0.0, 0.0, 0.0, 0.0,	// average for bins 1-4 for varianceNarrow normalized
						  0.0 };	*/		// number of points that voted
		
		// Info about voting (bins, images, point indeces) for ground truth
		//String outputFileName = "votingInfoIm" + Integer.parseInt(imageNum) + ".txt";
		String outputFileName = "votingInfoSingletons" + ".txt";
		
		// Open the the file named outputFileName and and pass BufferWrite instance when voting for each point
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(outputFileName));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("Writting to the file " + outputFileName);
		
		int queryImage;
		int queryPt;
		
		for(int j = 0; j < NUM_RUNS; j++)
		{	
			//Take one by one query image and query point
			queryImage = queryImages[j];
			queryPt = queryPts[j];
			
			System.out.println("Query image: " + queryImage);
			System.out.println("Query point: " + queryPt);
			
			//take one image and query all the points in it and vote...			
			NumberFormat formatter = new DecimalFormat("00000");			
			// get the path names from the Paths class
			//String imageFileName = Path.imageSinglObj + "ukbench" + formatter.format(queryImage) + ".obj";
			
			String imageFileName = Path.fileSinglOBJ_2 + "ukbench" + formatter.format(queryImage) + ".obj";
			//String imageFileName = Path.fileSinglOBJ_3 + "ukbench" + formatter.format(queryImage) + ".obj";
			
			double[] votes = kdTreeInstance.voteForImageWithGroundTruth(trainedTree, 
					populatedTree, counts, imageFileName, queryPt, IndexingGroundTruthTest2.NUM_IMAGES_POPULATING, writer); 				

		}
		
		System.out.println("Flushing and closing the writer.");
		
		try {
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}	
}
