package ca.utoronto.cs.kdtree.basic;
import java.io.*;
import java.text.*;
import java.util.*;

/**
 * This class represents the engine of the system. It has a main method
 * where the keypoints are extracted from the images, paired and sets of pairs
 * are written on the disk. The .key and .pgm files
 * are deleted from the system.
 * @author geri
 *
 */
public class KeypointPairsExtraction 
{
	/**
	 * The number of images we are working with.
	 */
	public static final int NUM_IMAGES = 10200;
	
	/**
	 * Maps that holds all the image file names and gives the mapping
	 * from these names to image ids (index of the file name in this vector)
	 */
	private static Vector<String> nameIdMap = new Vector<String>(KeypointPairsExtraction.NUM_IMAGES);
	
	
	
	public KeypointPairsExtraction()
	{
		for(int i = 3; i < KeypointPairsExtraction.NUM_IMAGES; i=i+4)
		{	
			//progress
			if((i % 100) == 0 || i == 1) 
			{
				System.out.println("Image number " + i);
			}
			
			NumberFormat formatter = new DecimalFormat("00000");
		
			String imageJPG = Path.imageJPG + "ukbench" + formatter.format(i) + ".jpg ";	
			String imagePGM = Path.imagePGM + "ukbench" + formatter.format(i) + ".pgm";		
			String fileNameKEY = Path.filePairOBJ_1 + "ukbench" + formatter.format(i) + ".key";			
			String fileNameOBJ = Path.filePairOBJ_1 + "ukbench" + formatter.format(i) + ".obj";
			//String fileNamePAIROBJ = KeypointPairsExtraction.imagePath + "ukbench" + formatter.format(i) + ".pair.obj";
			//String fileNameQueryPAIROBJ = KeypointPairsExtraction.imagePath + "ukbench" + formatter.format(i) + ".q.pair.obj";
			
			//String fileNamePAIRPCAKEY = KeypointPairsExtraction.imagePathNoBackupMeaningPair10PcaKey + "ukbench" + formatter.format(i) + ".pair5.pca.key";
			//String fileNamePAIROBJ = KeypointPairsExtraction.imagePathNoBackupMeaningObj + "ukbench" + formatter.format(i) + ".pair.obj";
			
			String fileNamePAIRPCAKEY = Path.filePairMix05KEY_4 + "imagemix" + formatter.format(i) + ".pair.key";
			String fileNamePAIROBJ = Path.filePairMix05OBJ_4 + "imagemix" + formatter.format(i) + ".pair.obj";
			
			// extracts keypoints from the images
			// this.generateKEYFile(imageJPG, imagePGM, fileNameKEY);		
			
			//extracting points using PCA			
			// this.generatePKEYfile(imagePGM, fileNameKEY, fileNamePKEY);						
			
			// creates the files that represent the sets of keypoints (one file per image).
			// this.generateOBJFile(fileNamePKEY, fileNameOBJ);
			
			// create the files that represent the sets of keypoint pairs (one file per image)
			// this.generatePAIROBJFile(fileNamePKEY, fileNamePAIROBJ);
			
			// create the files that represent the sets of keypoint pairs (one file per image)
			//this.generatePAIROBJFile(fileNameOBJ, fileNamePAIROBJ);
			this.generatePAIROBJFile(fileNamePAIRPCAKEY, fileNamePAIROBJ);
			
			
			// remove .key files		
			// this.deleteFile(fileNameKEY);				
			
			// add .obj file to the nameIdMap - we map image file names to image ids			
			//KeypointPairsExtraction.nameIdMap.add(fileNamePAIROBJ);
			KeypointPairsExtraction.nameIdMap.add(fileNamePAIROBJ);
		
		}
	}
	
	/**
	 * This constructor can be used to instantiate the class without running all the conversions but 
	 * with creating the name-id map for all the images
	 * @param helpingConstructor - anyhting can be passed since this argument is used only to distinguish
	 * 				this constructor from the other one
	 */
	public KeypointPairsExtraction(String helpingConstructor)
	{
		// won't use the passed parametter at all
		for(int i = 0; i < KeypointPairsExtraction.NUM_IMAGES; i++)
		{
			NumberFormat formatter = new DecimalFormat("00000");
			String fileNamePAIROBJ = Path.filePairOBJ_1 + "ukbench" + formatter.format(i) + ".pair.obj";
			KeypointPairsExtraction.nameIdMap.add(fileNamePAIROBJ);
		}
	}
	
	/**
	 * This constructor is just like the first one but it deals with the image number imageNum 
	 * instead of all of them
	 * @param numImages - which image to process
	 */
	public KeypointPairsExtraction(int imageNum)
	{		
		NumberFormat formatter = new DecimalFormat("00000");
		
		String imageJPG = Path.imageJPG + "ukbench" + formatter.format(imageNum) + ".jpg ";	
		String imagePGM = Path.imagePGM + "ukbench" + formatter.format(imageNum) + ".pgm";
		String fileNamePKEY = Path.fileSinglPKEY_1 + "ukbench" + formatter.format(imageNum) + ".pkey";
		String fileNameOBJ = Path.fileSinglOBJ_1 + "ukbench" + formatter.format(imageNum) + ".obj";
		String fileNamePAIROBJ = Path.filePairOBJ_1 + "ukbench" + formatter.format(imageNum) + ".pair.obj";
		//String fileNameQueryPAIROBJ = KeypointPairsExtraction.imagePath + "ukbench" + formatter.format(imageNum) + ".q.pair.obj";
		
		// extracts keypoints from the images
		//this.generateKEYFile(imageJPG, imagePGM, fileNameKEY);		
		
		//extracting points using PCA			
		//this.generatePKEYfile(imagePGM, fileNameKEY, fileNamePKEY);						
		
		// creates the files that represent the sets of keypoints (one file per image).
		//this.generateOBJFile(fileNamePKEY, fileNameOBJ);
		
		// create the files that represent the sets of keypoint pairs (one file per image)
		this.generatePAIROBJFile(fileNameOBJ, fileNamePAIROBJ);
		
		// remove .key files		
		//this.deleteFile(fileNameKEY);				
		
		// add .obj file to the nameIdMap - we map image file names to image ids			
		//KeypointPairsExtraction.nameIdMap.add(fileNamePAIROBJ);
	
	}
	
	
	/**
	 * The main method.
	 * @param args
	 */		
	public static void main(String[] args) 	
	{				
		long timeStart = System.currentTimeMillis();
		
		KeypointPairsExtraction kPts = new KeypointPairsExtraction(9305);
		
		long timeEnd = System.currentTimeMillis();
		System.out.println("Time to get data set for training the tree: " + ((timeEnd - timeStart)/1000) + "sec");
		System.out.println("Time to get data set for training the tree: " + ((timeEnd - timeStart)) + " milisec");
		
		//kPts.test();		
	}
	
	/**
	 * This method converts jpg to pgm format.
	 * The images are in jpg format. They are converted to pgm format using the
	 * 'convert' system call.
	 * Then the keypoints are extracted using the extractSift script; this script
	 * takes as an argument the number of images we are dealing with.
	 */
	
	
	public void generateKEYFile(String imageJPG, String imagePGM, String imageKEY)
	{
		try
		{
			Runtime runtime = Runtime.getRuntime();								
												
			// Convert images from .jpg to .pgm format
			Process convert = runtime.exec("convert " + imageJPG + imagePGM);
			convert.waitFor();										
			
			//Extract keypoints from the .pgm images using the shell script exractSift
			//extractSift extracts keypoints from one image at a time				
			String extractSift = "/h/291/geri/Workspace/indexing_ukbench_singletons/extractSift ";
			Process extract = runtime.exec(extractSift + imagePGM + " " + imageKEY);									
			extract.waitFor();							
			
		} catch(IOException ioe)
		{
			ioe.printStackTrace();
			
		} catch (InterruptedException ie) 
		{
			ie.printStackTrace();
		}			
	}
	
	
	/**
	 * 
	 * @param imagePGM
	 * @param fileNameKEY
	 * @param fileNamePKEY
	 */
	public void generatePKEYfile(String imagePGM, String fileNameKEY, String fileNamePKEY)	
	{
		try
		{
			Runtime runtime = Runtime.getRuntime();								
												
			// run the PCA extractor
			// recalckeys gpcavects.txt image_name.pgm image_name.key image_name.pkey
			Process generate = runtime.exec("recalckeys gpcavects.txt " + imagePGM + 
					" " + fileNameKEY + " " + fileNamePKEY);
			generate.waitFor();					
			
			System.out.println("Recalc for image " + imagePGM + " done!");
						
		} catch(IOException ioe)
		{
			ioe.printStackTrace();
			
		} catch (InterruptedException ie) 
		{
			ie.printStackTrace();
		}			
	}	
	
	
	/**
	 * This method is used to create the .obj files that represent the keypoints.
	 * Each file is written on the disk as a serialized object and it represents
	 * all the keypoints found in an image. The file has extension .obj.	 
	 */
	public void generateOBJFile(String pkeyFileName, String objFileName)
	{
			// Create one file per image						
			ImageKeypoints keypts = new ImageKeypoints(pkeyFileName);
			// Store the object on the disk under a name with extension .obj							
    		File objFile = new File(objFileName);    		    	 	
    
			try
			{
				objFile.delete();
				objFile.createNewFile();
				ObjectOutputStream out = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(objFile)));
				out.writeObject(keypts);		        	    		        	    	        	  
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
	
	/**
	 * This method is used to create the .pair.obj files that represent the keypoint pairs.
	 * Each file is written on the disk as a serialized object and it represents
	 * all the keypoint pairs found in an image. The file has extension .pair.obj.	 
	 */
	public ImageKeypointPairs generatePAIROBJFile(String fileName, String pairobjFileName)
	{
			// Create one file per image		
		    // constructor of ImageKeypointPairs takes pkeyFileName, extracts singletons and 
		    // creates pairs that are now in ImageKeypointPairs instance
			ImageKeypointPairs keyptPairs = new ImageKeypointPairs(fileName);
			
			// Store the object on the disk under a name with extension .pair.obj							
    		File pairobjFile = new File(pairobjFileName);    		    	 	
			try
			{
				pairobjFile.delete();
				pairobjFile.createNewFile();
				ObjectOutputStream out = new ObjectOutputStream(new
						BufferedOutputStream(new FileOutputStream(pairobjFile)));
				out.writeObject(keyptPairs);		        	    		        	    	        	  
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
			return null;
	}
	
	
	/**
	 * This method is used to delete the files from the disk that we do
	 * not need anymore. It takes the parameter extension that distinguishes which file
	 * to delete (since all the files have the same name and deiffer only by extension).
	 * @param extension
	 */
	public void deleteFile(String fileName)
	{
			File file = new File(fileName);
			file.delete();				
	}
	
	/**
	 * Getter for the imageNameMap
	 * @return
	 */
	public static Vector<String> getNameIdMap()
	{
		return KeypointPairsExtraction.nameIdMap;
	}
	
	/**
	 * Given the image number, returns the short[][] points set containing descriptors of that images.	
	 * It is used to get the argument for the trainKDTree method, for voting etc.
	 */	
	public static short[] getDataSet(int imageNum)
	{		
		NumberFormat formatter = new DecimalFormat("00000");
				
		short[] dataSet = new short[0];
		
		// get the path names from the path class
		String pairObjFileName = Path.filePairOBJ_3 + "ukbench" + formatter.format(imageNum) + ".pair.obj";
		
		File pairObjFile = new File(pairObjFileName);		
		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(pairObjFile)));
			ImageKeypointPairs imageKeyptPairs = (ImageKeypointPairs)in.readObject();	
			
			KeypointPair[] keyptPairArr = imageKeyptPairs.getKeyptPairArray();	
						
			dataSet = new short[Keypoint.DESCRIPTOR_LENGTH * keyptPairArr.length];	
			
			int position = 0;			
			for(int k = 0; k < keyptPairArr.length; k++)
			{	
				short[] desc = keyptPairArr[k].getDescriptor();	
				for(int j = 0; j < desc.length; j++, position++)
				{
					dataSet[position] = desc[j]; 						
				}					
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
				
	return dataSet;
	}
	
	/**
	 * Used for testing only.
	 */
	private void test()
	{	
				
	}
}


