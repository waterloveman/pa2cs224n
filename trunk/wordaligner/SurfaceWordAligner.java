package cs224n.wordaligner;  

import cs224n.util.*;
import java.util.List;

/**
* Simple alignment baseline which maps french positions to english positions.
* If the french sentence is longer, all final word map to null.
*/
public class SurfaceWordAligner extends WordAligner
{
	private CounterMap<String,String> dummy;
	private Counter<String> english;
	private Counter<String> french;
	private double englishTotal;
	private double englishCount;
	private double frenchTotal;
	private double frenchCount;
	private double total;
	private double count;
	
	public Alignment alignSentencePair(SentencePair sentencePair)
	{
		Alignment alignment = new Alignment();
		int numFrenchWords = sentencePair.getFrenchWords().size();
		int numEnglishWords = sentencePair.getEnglishWords().size();
		for (int frenchPosition = 0; frenchPosition < numFrenchWords; frenchPosition++)
		{
			String f = sentencePair.getFrenchWords().get(frenchPosition);
			int bestEngPos = -1;
			double bestProb = -1;
			for (int englishPosition = 0; englishPosition < numEnglishWords; englishPosition++)
			{
				String e = sentencePair.getEnglishWords().get(englishPosition);
				double prob = getWordAlignmentProb(f, e);
				if(prob > bestProb) {
					bestProb = prob;
					bestEngPos = englishPosition;
				}
			}
			alignment.addAlignment(bestEngPos, frenchPosition, true);
		}
		return alignment;
	}
	
	public double getWordAlignmentProb(String targetWord, String sourceWord)
	{
		double count = english.getCount(targetWord);
		if(count == 0) {
			return 0;
		}
		double p = dummy.getCount(sourceWord, targetWord) / count;
		return p;
	}
	
	public double getAlignmentProb(List<String> targetSentence, List<String> sourceSentence, Alignment alignment)
	{
		return 0;
	}
	
	public CounterMap<String,String> getProbSourceGivenTarget()
	{
		return dummy;
	}
	
	public void train(List<SentencePair> trainingPairs)
	{
		dummy = new CounterMap<String,String>();
		english = new Counter<String>();
		french = new Counter<String>();
		for(SentencePair pair : trainingPairs){
			List<String> targetWords = pair.getEnglishWords();
			List<String> sourceWords = pair.getFrenchWords();
			for(String source : sourceWords)
			{
				for(String target : targetWords)
				{
					if(dummy.getCount(source, target) == 0)
					{
						dummy.setCount(source, target, 1.0);
					}
					else
					{
						dummy.incrementCount(source, target, 1.0);
					}
				}
				if(french.getCount(source) == 0)
				{
					french.setCount(source, 1.0);
				}
				else
				{
					french.incrementCount(source, 1.0);
				}
			}
			for(String target : targetWords)
			{
				if(english.getCount(target) == 0)
				{
					english.setCount(target, 1.0);
				}
				else
				{
					english.incrementCount(target, 1.0);
				}
			}
		}
		englishTotal = english.totalCount();
		englishCount = english.size();
		frenchTotal = french.totalCount();
		frenchCount = french.size();
		total = dummy.totalCount();
		count = dummy.size();
	}	
}
