import java.io.*;
import java.text.*;
import java.util.*;

/**
 * This class represents the engine of the system. It has a main method
 * where the keypoints get extracted from the images and the ImageKeypoints 
 * objects are created and serialized on the disk. The .key and .pgm files
 * are deleted from the system.
 * @author geri
 *
 */
public class KeypointsExtraction 
{
	/**
	 * The number of images we are working with.
	 */
	public static final int NUM_IMAGES = 10200;
	
	/**
	 * Maps that holds all the image file names and gives the mapping
	 * from these names to image ids (index of the file name in this vector)
	 */
	private static Vector<String> nameIdMap = new Vector<String>(KeypointsExtraction.NUM_IMAGES);
	
	public KeypointsExtraction()
	{
		for(int i = 3; i < KeypointsExtraction.NUM_IMAGES; i=i+4)
		{	
			//progress
			if((i % 100) == 0) 
			{
				System.out.println("Image number " + i);
			}
			
			NumberFormat formatter = new DecimalFormat("00000");
		
			//String imageJPG = Path.imageJPG + "ukbench" + formatter.format(i) + ".jpg ";	
			//String imagePGM = Path.imagePGM + "ukbench" + formatter.format(i) + ".pgm";		
			//String fileNameKEY = Path.fileNameKEY + "ukbench" + formatter.format(i) + ".key";
			//String fileNamePKEY = Path.fileNamePKEY + "ukbench" + formatter.format(i) + ".pkey";
			//String fileNamePCAKEY = Path.fileNamePCAKEY + "ukbench" + formatter.format(i) + ".pca.key";
					
			String fileNamePAIRPCAKEY = Path.fileSinglMix05PCAKEY_3 + "imagemix" + formatter.format(i) + ".pca.key";
			String fileNameOBJ = Path.fileSinglMix05OBJ_3 + "imagemix" + formatter.format(i) + ".obj";
			
			// extracts keypoints from the images
			//this.generateKEYFile(imageJPG, imagePGM, fileNameKEY);		
			
			//extracting points using PCA-SIFT			
			//this.generatePKEYfile(imagePGM, fileNameKEY, fileNamePKEY);						
			
			// creates the files that represent the sets of keypoints (one file per image).
			this.generateOBJFile(fileNamePAIRPCAKEY, fileNameOBJ);
			
			// remove .key files		
			//this.deleteFile(fileNameKEY);				
			
			// add .obj file to the nameIdMap - we map image file names to image ids	
			KeypointsExtraction.nameIdMap.add(fileNameOBJ);
		}
	}
	
	/**
	 * This constructor can be used to instantiate the class without running all the conversions but 
	 * with creating the name-id map for all the images
	 * @param helpingConstructor - anyhting can be passed since this argument is used only to distinguish
	 * 				this constructor from the other one
	 */
	public KeypointsExtraction(String helpingConstructor)
	{
		// we don't use the passed parametter at all
		for(int i = 0; i < KeypointsExtraction.NUM_IMAGES; i++)
		{
			NumberFormat formatter = new DecimalFormat("00000");
			String fileNameOBJ = Path.fileSinglOBJ_1 + "ukbench" + formatter.format(i) + ".obj";
			KeypointsExtraction.nameIdMap.add(fileNameOBJ);
		}
	}
	
	/**
	 * This constructor is just like the first one but it deals with the image number imageNum 
	 * instead of all of them
	 * @param numImages - which image to process
	 */
	public KeypointsExtraction(int imageNum)
	{		
		NumberFormat formatter = new DecimalFormat("00000");
		
		String imageJPG = Path.imageJPG + "ukbench" + formatter.format(imageNum) + ".jpg ";	
		String imagePGM = Path.imagePGM + "ukbench" + formatter.format(imageNum) + ".pgm";
		String fileNameKEY = Path.fileSinglKEY_2 + "ukbench" + formatter.format(imageNum) + ".key";
		String fileNamePKEY = Path.fileSinglPKEY_1 + "ukbench" + formatter.format(imageNum) + ".pkey";
		String fileNameOBJ = Path.fileSinglOBJ_1 + "ukbench" + formatter.format(imageNum) + ".obj";
		
		// extracts keypoints from the images
		//this.generateKEYFile(imageJPG, imagePGM, fileNameKEY);			
		
		//extracting points using PCA			
		//this.generatePKEYfile(imagePGM, fileNameKEY, fileNamePKEY);		
							
		// creates the files that represent the sets of keypoints (one file per image).
		this.generateOBJFile(fileNamePKEY, fileNameOBJ);
		
		// remove .key files		
		//this.deleteFile(fileNameKEY);				
		
		// add .obj file to the nameIdMap - we map image file names to image ids	
		//KeypointsExtraction.nameIdMap.add(fileNameOBJ);
	}
	
	
	/**
	 * The main method.
	 * @param args
	 */		
	public static void main(String[] args) 	
	{		
		//KeypointsExtraction kPts = new KeypointsExtraction("");
		//KeypointsExtraction kPts = new KeypointsExtraction(10048);

		//kPts.test();		
	}
	
	/**
	 * This method converts jpg to pgm format.
	 * The images are in jpg format. They are converted to pgm format using the
	 * 'convert' system call.
	 * Then they keypoints are extracted using the extractSift script; this script
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
	 * This method is used to create the .obj files that represent the keypoints.
	 * It takes as a parameter the .pkey or .pcakey file -  and it acts the same on both.
	 * Each file is written on the disk as a serialized object and it represents
	 * all the keypoints found in an image. The file has extension .ojb.	 
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
		return KeypointsExtraction.nameIdMap;
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
		String objFileName = Path.fileSinglOBJ_3 + "ukbench" + formatter.format(imageNum) + ".obj";
		
		// System.out.println(objFileName);
					
		File objFile = new File(objFileName);		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(objFile)));				
			
			ImageKeypoints imageKeypts = (ImageKeypoints)in.readObject();	
							
			KeypointSingleton[] keyptArr = imageKeypts.getKeyptArray();			
						
			dataSet = new short[Keypoint.DESCRIPTOR_LENGTH * keyptArr.length];						
			
			int position = 0;
			for(int k = 0; k < keyptArr.length; k++)
			{					
				short[] desc = keyptArr[k].getDescriptor();					
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
		for(int k = 0; k < KeypointsExtraction.NUM_IMAGES; k++)
		{
			NumberFormat formatter = new DecimalFormat("00000");
			String objFileName = Path.fileSinglOBJ_1 + "ukbench" + formatter.format(k) + ".obj";			
			
			File objFile = new File(objFileName);		
			try
			{						
				ObjectInputStream in = new ObjectInputStream(new
						BufferedInputStream(new FileInputStream(objFile)));
				ImageKeypoints keypoints = (ImageKeypoints)in.readObject();	
				KeypointSingleton[] keyptArr = keypoints.getKeyptArray();
				
				if(keyptArr.length == 0)
				{
					System.out.println("Image has no points - bug in pca recalc!");
				}else
				{
					System.out.println("length of keyptArr is " + keyptArr.length);
                    // first keypoint in this image			
					KeypointSingleton pt1 = keyptArr[0];

					short[] desc1 = pt1.getDescriptor();
					System.out.println("The length of the descriptor is " + desc1.length);
					System.out.print("The description values of the first point are: ");
					for(int i = 0; i < desc1.length; i++)
					{
						System.out.print(desc1[i] + " ");
					}
				}									
				
			    in.close();
			}
			catch (IOException ioe)
			{
				ioe.printStackTrace();
			}
			catch (ClassNotFoundException cnfe)
			{ 
				cnfe.printStackTrace();
			}	
		}
	}
}


