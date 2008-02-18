/*
 * The class Scores represents a collection of ImageScorePairs. It has a private attribute 'scores' 
 * which is a linked list that contains all the ImageScorePairs passed to the constructor of this class.
 * It also has a method used to sorting ImageScorePairs based on their scores in descending order.
 */

import java.util.*;

public class Scores 
{
	private LinkedList<ImageScorePair> scores;
	
	public Scores(LinkedList<ImageScorePair> scoresList)
	{
		this.scores = scoresList;
	}

	public LinkedList<ImageScorePair> getScoresList()
	{
		return this.scores;
	}
	
	public Scores sortScores(Scores scores)
	{
		LinkedList<ImageScorePair> sortedList = this.mergeSort(scores.scores);
		return new Scores(sortedList);
	}
	
	public LinkedList<ImageScorePair> mergeSort(LinkedList<ImageScorePair> scores)
	{
		LinkedList<ImageScorePair> left = new LinkedList<ImageScorePair>();
		LinkedList<ImageScorePair> right = new LinkedList<ImageScorePair>();
		LinkedList<ImageScorePair> scoresSorted;
		
		if(scores.size() <= 1)
		{
			scoresSorted = scores;
		}else
		{
			int middle = (int)Math.floor(scores.size() / 2);
			
			for(int i = 0; i < middle; i++)
			{
				left.add(scores.get(i));
			}
			
			for(int i = middle; i < scores.size(); i++)
			{
				right.add(scores.get(i));
			}
			
			left = mergeSort(left);
			right = mergeSort(right);
			
			scoresSorted = merge(left, right);
		}
		
		return scoresSorted;
	}
    
	public LinkedList<ImageScorePair> merge(LinkedList<ImageScorePair> left, LinkedList<ImageScorePair> right)
	{
		LinkedList<ImageScorePair> mergedList = new LinkedList<ImageScorePair>();
		
		int i = 0, j = 0;
		while(i < left.size() && j < right.size())
		{
			if(left.get(i).score >= right.get(j).score)
			{
				mergedList.add(left.get(i));
				i++;
			}else
			{
				mergedList.add(right.get(j));
				j++;
			}
		}
		
		if(i < left.size())
		{
			for(; i < left.size(); i++)
			{
				mergedList.add(left.get(i));
			}
		}
		
		if(j < right.size())
		{
			for(; j < right.size(); j++)
			{
				mergedList.add(right.get(j));
			}
		}
		
		return mergedList;
	}
}
