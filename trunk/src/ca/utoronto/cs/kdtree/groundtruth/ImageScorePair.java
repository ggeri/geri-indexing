package ca.utoronto.cs.kdtree.groundtruth;
/**
 * The class ImageScorePair represents a pair that consists of an image and its score.
 * It is used for the purpose of calculating scores.
 * 
 * @author geri
 *
 */
public class ImageScorePair 
{
	int imageNo;
	double score;
	
	public ImageScorePair(int no, double score)
	{
		this.imageNo = no;
		this.score = score;
	}
	
	public double getScore()
	{
		return this.score;
	}
	
	public int getImageNo()
	{
		return this.imageNo;
	}
	
}
