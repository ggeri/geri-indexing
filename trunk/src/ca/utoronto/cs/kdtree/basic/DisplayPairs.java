package ca.utoronto.cs.kdtree.basic;
import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * This class is responsible for displaying the created pairs. 
 * It gets the spatial information of the pairs in the image of interest and writes it to
 * a file. Then a matlab script is called what reads that file and displays the image with the pairs on it.
 * @author geri
 *
 */

public class DisplayPairs 
{
    public static void main(String[] args)
    {
    	DisplayPairs dispPairs = new DisplayPairs();
    	dispPairs.writeSpatialInfo(0, "src\\matlabData\\spatialInfoIm0.txt");
    	
    }
	
	
    /**
     * This method reads the .pair.obj file and gets the spatial info of all the pairs
     * in that file. It, then writes them out to a file (so that matlab can read them from there).
     */
	public void writeSpatialInfo(int imageNum, String fileName)
	{	
		int count = 0;
		
		System.out.println("matlab file to write data for matlab to: " + fileName);
		
		NumberFormat formatter = new DecimalFormat("00000");
		
		String pairObjFileName = Path.filePairOBJ_2 + "ukbench" + formatter.format(imageNum) + ".pair.obj";
			
		File pairObjFile1 = new File(pairObjFileName);
		
		// file to write to
		File file1 = new File(fileName);		
		
		try
		{						
			ObjectInputStream in = new ObjectInputStream(new
					BufferedInputStream(new FileInputStream(pairObjFile1)));
			
			ImageKeypointPairs imageKeyptPairs = (ImageKeypointPairs)in.readObject();
							
			KeypointPair[] keyptPairArr = imageKeyptPairs.getKeyptPairArray();
			
			Writer out = new BufferedWriter(new FileWriter(file1));			
			
			System.out.println("Length of the keypoint pair array is " + keyptPairArr.length);
			
			for(int i = 0; i < keyptPairArr.length; i++)
			{
				KeypointPair pair = keyptPairArr[i];
				double[] loc1 = pair.getKeyptLocationOne();
				double[] loc2 = pair.getKeyptLocationTwo();
				
				for(int j = 0; j < loc1.length; j++)
				{													    	
					out.write((new Double(loc1[j])).toString());
					out.write('\t');
					count++;
				}
				for(int j = 0; j < loc2.length; j++)
				{													    	
						out.write((new Double(loc2[j])).toString());
						out.write('\t');
						count++;					
				}
				
			}		
			
			out.close();
			in.close();
			
			// test
			
			File file = new File("src\\matlabData\\spatialInfoIm0.txt");			
			BufferedReader in2 = new BufferedReader(new FileReader(file));
			String line = in2.readLine();
			String[] parsedLine = line.split(" ");
			System.out.println("Length of the parsed line is " + parsedLine.length);
			
			in.close();
			
		}catch (IOException ioe)
		{
		ioe.printStackTrace();
		}
		catch (ClassNotFoundException cnfe)
		{ 
			cnfe.printStackTrace();
		}
		
		System.out.println("Number of spatial elements is " + count);
	}
}
