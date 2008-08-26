package ca.utoronto.cs.kdtree.groundtruth;
/*
 * The class Scores represents a collection of ImageScorePairs. It has a private attribute 'scores' 
 * which is a linked list that contains all the ImageScorePairs passed to the constructor of this class.
 * It also has a method used to sorting ImageScorePairs based on their scores in descending order.
 */

import java.util.*;

public class Scores 
{
	private LinkedList<ImageScorePair> scoreList;		
	
	
	// constructor
	public Scores(LinkedList<ImageScorePair> scoresList)
	{
		this.scoreList = scoresList;
	}

	// getter for 'scores' list
	public LinkedList<ImageScorePair> getScoreList()
	{
		return this.scoreList;
	}
	
	// sorts the 'scores' list
	public void sortDescending()
	{
		this.scoreList = this.mergeSort(this.scoreList);
		
	}
	
	// takes the Scores object and sorts it's list of scores
	public Scores sortScores(Scores scores)
	{
		LinkedList<ImageScorePair> sortedList = this.mergeSort(scores.scoreList);
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
	
	
	// for testing
	public static void main(String[] args)
	{
		test();
	}
	
	private static void test()
	{
		LinkedList<ImageScorePair> scoreList = new LinkedList<ImageScorePair>();
		ImageScorePair p1 = new ImageScorePair(0,3);
		ImageScorePair p2 = new ImageScorePair(1,1);
		ImageScorePair p3 = new ImageScorePair(2,5);		
		scoreList.add(p1);
		scoreList.add(p2);
		scoreList.add(p3);
		System.out.println("Unsorted");
		for(int i = 0; i < 3; i++)
		{
			System.out.println(scoreList.get(i).getImageNo() + " " + scoreList.get(i).getScore());
		}
		
		Scores scores = new Scores(scoreList);
		scores.sortDescending();
		LinkedList<ImageScorePair> scoreListSorted = scores.getScoreList();
		System.out.println("Sorted");
		for(int i = 0; i < 3; i++)
		{
			System.out.println(scoreListSorted.get(i).getImageNo() + " " + scoreListSorted.get(i).getScore());
		}
		scoreListSorted.remove(p3);
		
		System.out.println("Sorted without first element");
		for(int i = 0; i < 3; i++)
		{
			System.out.println(scoreListSorted.get(i).getImageNo() + " " + scoreListSorted.get(i).getScore());
		}
		
	}
}
