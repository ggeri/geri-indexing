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
	int score;
	
	public ImageScorePair(int no, int score)
	{
		this.imageNo = no;
		this.score = score;
	}
	
	public int getScore()
	{
		return this.score;
	}
	
	public int getImageNo()
	{
		return this.imageNo;
	}
	
}
