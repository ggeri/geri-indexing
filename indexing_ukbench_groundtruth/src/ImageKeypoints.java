import java.io.*;

/**
 * The instance of this class is an object that contains all the
 * keypoints that belong to one image.
 * @author geri
 *
 */
public class ImageKeypoints implements Serializable
{

	/**
	 * For serialization
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * The array that contains all the keypoints of this image.
	 */
	private KeypointSingleton[] keyptArray;
	
	/**
	 * The constructor that creates an instance of this class. The parameter keyFileName
	 * is the string that represents the name of the image file. The constructor invokes 
	 * the method collectKeypoints which opens the ASCII file, reads it line by line, 
	 * creates the keypoint and stores them in the keyptArray.
	 * @param keyFileName
	 */
	public ImageKeypoints(String keyFileName)
	{
		this.collectKeypoints128(keyFileName);
	}
		
	/**
	 * This method takes the passed file name, opens the pkey file, reads it line by line and 
	 * creates keypoints and populates with them the keyptArray attribute of this class.
	 * @param pkeyFileName
	 */
	public void collectKeypoints(String pkeyFileName)
	{		
		try
		{
			// if file doesn't exist, create empty keyptArray and return it
			if(!new File(pkeyFileName).exists())
			{
				System.out.println("File didn't exist - we generate an empty keyptArray");
				this.keyptArray = new KeypointSingleton[0];
			}else
			{
				BufferedReader in = new BufferedReader(new FileReader(pkeyFileName));
				String line = in.readLine();
				String[] parsedLine = line.split("\\s");
				int numOfKeypts = Integer.parseInt(parsedLine[0]);

				// initialize the keypoint array
				this.keyptArray = new KeypointSingleton[numOfKeypts];
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
					String[] wholeParsedLine = new String[Keypoint.DESCRIPTOR_LENGTH];
					
					for (int j = 0; j < 3; j++)
					{
						line = in.readLine();
						//parsedLine array is of length 13 but it has 12 tokens
						parsedLine = line.split("\\s");			

						int off = j * 12;									
						for(int m = 0; m < parsedLine.length; m++)  // for .pkey file, change m to m=1
						{						
							int off2 = off + m; // use m-1 fir .pkey files
							wholeParsedLine[off2] = parsedLine[m];	
						}												
					}
					short[] descr = new short[Keypoint.DESCRIPTOR_LENGTH];									
					
					for (int k = 0; k < Keypoint.DESCRIPTOR_LENGTH; k++)
					{																
						descr[k] = Short.parseShort(wholeParsedLine[k]);
					}
					//set the vector in the object				
					keypt.setDescriptor(descr);

					//add the pt to array of points	
					this.keyptArray[i] = keypt;
				}
				in.close();
			}
			
		}
		catch (FileNotFoundException e)
		{
			System.out.println("This file was not found; we create an empty keyptArray");
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * This method takes the passed file name, opens the key file that contains 128 dimensional features, 
	 * reads it line by line and creates keypoints and populates with them the keyptArray attribute of 
	 * this class.
	 * @param pkeyFileName
	 */
	public void collectKeypoints128(String keyFileName)
	{		
		try
		{
			// if file doesn't exist, create empty keyptArray and return it
			if(!new File(keyFileName).exists())
			{
				System.out.println("File didn't exist - we generate an empty keyptArray");
				this.keyptArray = new KeypointSingleton[0];
			}else
			{
				BufferedReader in = new BufferedReader(new FileReader(keyFileName));
				String line = in.readLine();
				String[] parsedLine = line.split("\\s");
				int numOfKeypts = Integer.parseInt(parsedLine[0]);

				// initialize the keypoint array
				this.keyptArray = new KeypointSingleton[numOfKeypts];
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
					String[] wholeParsedLine = new String[Keypoint.DESCRIPTOR_LENGTH];
					
					for (int j = 0; j < 6; j++)
					{
						line = in.readLine();
						//parsedLine array is of length 20 but it has 21 tokens
						parsedLine = line.split("\\s");			

						int off = j * 20;									
						for(int m = 1; m < parsedLine.length; m++)  // for .pkey and .key file, change m to m=1
						{						
							int off2 = off + m-1; // use m-1 fir .pkey and .key files
							wholeParsedLine[off2] = parsedLine[m];	
						}												
					}
					line = in.readLine();
					//parsedLine array is of length 8 but it has 9 tokens
					parsedLine = line.split("\\s");			

					int off = 6 * 20;									
					for(int m = 1; m < parsedLine.length; m++)  // for .pkey and .key file, change m to m=1
					{						
						int off2 = off + m-1; // use m-1 fir .pkey and .key files
						wholeParsedLine[off2] = parsedLine[m];	
					}				
				
					short[] descr = new short[Keypoint.DESCRIPTOR_LENGTH];									
					
					for (int k = 0; k < Keypoint.DESCRIPTOR_LENGTH; k++)
					{																
						descr[k] = Short.parseShort(wholeParsedLine[k]);
					}
					//set the vector in the object				
					keypt.setDescriptor(descr);

					//add the pt to array of points	
					this.keyptArray[i] = keypt;
				}
				in.close();
			}
			
		}
		catch (FileNotFoundException e)
		{
			System.out.println("This file was not found; we create an empty keyptArray");
			e.printStackTrace();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Get method for the attribute keyptArray
	 */
	public KeypointSingleton[] getKeyptArray()
	{
		return this.keyptArray;
	}
	
	/**
	 * This method takes all the points of this image and normalizes their
	 * feature vectors.
	 */
	public void normalizeKeypoints() 
	{
		for(int i = 0; i < keyptArray.length; i++ )
		{
			KeypointSingleton keypt = keyptArray[i];
			short[] desc = keypt.getDescriptor();
			int magn = 0;
			for(int j = 0; j < desc.length; j++)
			{
				magn = magn + desc[j];
			}
			if(magn == 0)
			{
				System.out.println("magn: " + magn);
				System.out.print("desc: ");
				for(int j = 0; j < desc.length; j++)
				{
					System.out.print(desc[j] + " ");
				}			
				System.out.println();
			}			
			
			for(int j = 0; j < desc.length; j++)
			{
				desc[j] = (short) Math.floor(desc[j] / magn);
			}
			keypt.setDescriptor(desc);
		}
	}
	
}
