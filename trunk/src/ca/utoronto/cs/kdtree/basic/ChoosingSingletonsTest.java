package ca.utoronto.cs.kdtree.basic;
import java.util.Random;

/**
 * This class is used to test if we are choosing the right pairs. We do the test on singletons 
 * and on pairs to be able to compare the numbers at the end.
 * First we use 36D singletons that we have already created from the software. We do two things:
 * 1. we use two images that belong to the same set of 4 and look for the fraction of the points 
 * in image_1 that have at least one match in the image_2.
 * 2. we use two images that belong to two different sets of 4 and look for the fraction of the points 
 * in image_1 that have at least one match in the image_2.
 * Then we use 36D pairs that we already have from the software (here we have 1 pair per singleton). 
 * We now do four different things:
 * 1. we use two images that belong to the same set of 4 and look for the fraction of the points 
 * in image_1 that have at least one match in the image_2.
 * 2. we use two images that belong to two different sets of 4 and look for the fraction of the points 
 * in image_1 that have at least one match in the image_2.
 * Now we extract 5 or 10 pairs per singleton and repeat the same experiment as in 1 and 2 above.
 * @author geri
 *
 */
public class ChoosingSingletonsTest {

	// max euclidean distance that we can tolerate between two points to still consider them the same
	private static double maxDistance = 4500; 
	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		
		long timeStart = System.currentTimeMillis();					
		
		ChoosingSingletonsTest cpTest = new ChoosingSingletonsTest();
		double fractionOfMatchPts = 0;
		double fractionSum = 0;
		int numImages = 1;

		// the pairs of images from the same set and average the results	
		for(int i = 0; i < KeypointsExtraction.NUM_IMAGES-1; i = i + 4)
		{
			int image1 = i;
			int image2 = i + 1;
			fractionOfMatchPts = cpTest.numSamePoints(image1, image2);				
			
			if(fractionOfMatchPts != -1)
			{
				fractionSum = fractionSum + fractionOfMatchPts;
				numImages++;
			}			
		}
		
		System.out.println("Two images from the same set");
		System.out.format("On average, fraction of points in the first image that "
				+ "have at least one match in the second image is: %.3f", (fractionSum / numImages));
		System.out.println(" %");	
		
		
		fractionOfMatchPts = 0;
		fractionSum = 0;
		numImages = 1;
		Random rand = new Random(1000);
		
		// take pairs of images from different sets and average the resutls
		for(int i = 0; i < KeypointsExtraction.NUM_IMAGES-1; i = i + 4)
		{			
			int image3 = i;
			int randNum = rand.nextInt(KeypointsExtraction.NUM_IMAGES);
			while(randNum == i || randNum == i + 1 || randNum == i + 2 || randNum == i + 3)
			{
				randNum = rand.nextInt(KeypointsExtraction.NUM_IMAGES);
			}
			int image4 = randNum;			
			fractionOfMatchPts = cpTest.numSamePoints(image3, image4);
			if(fractionOfMatchPts != -1)
			{
				fractionSum = fractionSum + fractionOfMatchPts;
				numImages++;
			}			
		}
		
		System.out.println("\nTwo images from different sets");
		System.out.format("On average, fraction of points in the first image that "
				+ "have at least one match in the second image is: %.3f", (fractionSum / numImages));
		System.out.println(" %");

		long timeEnd = System.currentTimeMillis();
		System.out.println("Time to get statistics for all images: " + (((timeEnd - timeStart)/1000)/60) + "min");
	}
	
	// returns the fraction of points in image 1 that had at least one match in image 2
	private double numSamePoints(int image1, int image2)
	{		
		// number of same points		
		double numSamePoints = 0;

		// get the data set for both images
		short[] image1DataSet = KeypointsExtraction.getDataSet(image1);
		short[] image2DataSet = KeypointsExtraction.getDataSet(image2);
		
		// if one of the images has no points in itself, return -1 as a sign not to count that one in
		if(image1DataSet.length == 0 || image2DataSet.length == 0)
		{
			return -1;
		}
		// compare the points in the two images
		// increment i by 36 since all the points are in the same vector - a new one starts every 36 cells
		
		for(int i = 0; i < image1DataSet.length; i = i + KeypointSingleton.DESCRIPTOR_LENGTH)
		{
			short[] point1 = new short[KeypointSingleton.DESCRIPTOR_LENGTH];			
			short[] point2 = new short[KeypointSingleton.DESCRIPTOR_LENGTH];
			
			//get point 1
			for(int j = 0; j < KeypointSingleton.DESCRIPTOR_LENGTH; j++)
			{				
				point1[j] = image1DataSet[i+j];				
			}
			
			// we want fraction of points in image 1 that have at least one match in image 2 so
			// we set 'match' to true if we find at least one match and the later if 'match' is 
			// true, we increment numSamePoints by one
			boolean match = false;
			// compare point 1 to each point in image 2
			for(int k = 0; k < image2DataSet.length; k = k + KeypointSingleton.DESCRIPTOR_LENGTH)
			{			
				// get point 2
				for(int l = 0; l < KeypointSingleton.DESCRIPTOR_LENGTH; l++)
				{				
					point2[l] = image2DataSet[k+l];		
				}
				// compare point 1 and point 2 and if the same, if they are the same, set 'match' to true
				if(this.distance(point1, point2) < ChoosingSingletonsTest.maxDistance)
				{					
					match = true;
				}
			}		
			// if we found at least one match, increment numSamePoints by one
			if(match)
			{
				numSamePoints++;
			}
		}
		
		int im1NumPts = image1DataSet.length / 36;
		double fraction = (numSamePoints / im1NumPts) * 100; 
		
		return fraction;
	}
	
	private double distance(short[] point1, short[] point2)
	{
		double distSqr = 0;
		for(int i = 0; i < point1.length; i++)
		{
			double temp = (point1[i] - point2[i])*(point1[i] - point2[i]);
			distSqr = distSqr + temp;
		}
		return Math.sqrt(distSqr);
	}
}
