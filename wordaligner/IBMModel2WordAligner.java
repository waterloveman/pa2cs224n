package cs224n.wordaligner;  

import cs224n.util.*;
import java.util.List;

/**
 * Simple alignment baseline which maps french positions to english positions.
 * If the french sentence is longer, all final word map to null.
 */
public class IBMModel2WordAligner extends WordAligner {
	
	private CounterMap<String,String> alignmentProbs, alignmentCounts;
	//private double[][][][] dProbs, dCounts;
	private double[] dProbs, dCounts;
	public static int MAX_SENTENCE_LENGTH = 3;
	
	public Alignment alignSentencePair(SentencePair sentencePair) {
		Alignment alignment = new Alignment();
		List<String> sourceWords = sentencePair.getFrenchWords();
		List<String> targetWords = sentencePair.getEnglishWords();
		for (int i = 0; i < sourceWords.size(); i++) {
			String curSourceWord = sourceWords.get(i);
			int maxAlignPos = 0;
			double maxAlignProb = alignmentProbs.getCount(curSourceWord, targetWords.get(0));
			//System.out.println("Cur prob for " + curSourceWord + " " + targetWords.get(0) + ":" + maxAlignProb);
			for (int j = 1; j < targetWords.size(); j++) {
				double curAlignProb = alignmentProbs.getCount(curSourceWord, targetWords.get(j));
				//System.out.println("Cur prob for " + curSourceWord + " " + targetWords.get(j) + ":" + curAlignProb);
				if (curAlignProb > maxAlignProb) {
					maxAlignPos = j;
					maxAlignProb = curAlignProb;
				}
			}
			if (alignmentProbs.getCount(curSourceWord, NULL_WORD) > maxAlignProb) {
				//System.out.println("Max align: -1");
				alignment.addAlignment(-1, i, true);
			}
			else {
				//System.out.println("Max align: " + maxAlignPos);
				alignment.addAlignment(maxAlignPos, i, true);
			}
		}
		return alignment;
	}
	
	
	public double getAlignmentProb(List<String> targetSentence, List<String> sourceSentence, Alignment alignment) {
		//Can be thought of as returning P(S, A|T) where T = target, S = source and A is the alignment
		//Simplifies to (sigma/(length_f +1)^(length_e)) * (product for each word in alignment->) * P(S_i | T_ai)
		double normalFactor = 1 / Math.pow(targetSentence.size() + 1, sourceSentence.size());
		double pairProduct = 1.0;
		for (int i = 0; i < sourceSentence.size(); i++) {
			String sourceWord = sourceSentence.get(i);
			System.out.println("ALIGNMENT IS: " + alignment.getAlignedTarget(i));
			int curAlignTarget = alignment.getAlignedTarget(i);
			if (curAlignTarget == -1)
				pairProduct *= alignmentProbs.getCount(sourceWord, NULL_WORD);
			else
				pairProduct *= alignmentProbs.getCount(sourceWord, targetSentence.get(curAlignTarget));
		}
		double thingy = normalFactor * pairProduct;
		System.out.println("Prob for:\n" + targetSentence + "\n" + sourceSentence + "\np: " + thingy);
		return normalFactor * pairProduct;
	}

    
	public CounterMap<String,String> getProbSourceGivenTarget(){
		//Return the counter map storing P(s|t) for each french/english word pair
		return alignmentProbs;
	}
	
	private int bucket(double sourcePos, double targetPos, double sourceLen, double targetLen)
	{
		
		//TODO: create bucketed value
		//should be a static function (i.e. not learned)
		if(sourceLen != targetLen) {
			System.out.print("");
		}
		double bval = sourcePos - targetPos * (sourceLen / targetLen);
		int i = (int) Math.abs(Math.round(bval));
		if(i <= 3) {
			return 0;
		}
		else if(i <= 6) {
			return 1;
		}
		else if(i <= 9) {
			return 2;
		}
		return 3;
	}

	public void train(List<SentencePair> trainingPairs) {
		alignmentProbs = new CounterMap<String,String>();
		alignmentCounts = new CounterMap<String,String>();
		dProbs = new double[MAX_SENTENCE_LENGTH+2];
		dCounts = new double[MAX_SENTENCE_LENGTH+2];
		double totalD = 0.0;
		
		Counter<String> targetWordCounts = new Counter<String>();
		
		//Get total count of word pairs
		for (SentencePair pair : trainingPairs) {
			List<String> targetWords = pair.getEnglishWords();
			List<String> sourceWords = pair.getFrenchWords();
			for (String source : sourceWords) {
				for (String target : targetWords) {
					alignmentCounts.incrementCount(source, target, 1.0);
				}
				alignmentCounts.incrementCount(source, NULL_WORD, 1.0);
			}
		}
		
		//Initialize P(s|t) to uniform values
		double wordPairCount = alignmentCounts.totalSize();
		double uniformInitValue = 1 / wordPairCount;
		for (String source : alignmentCounts.keySet()) {
			for (String target : alignmentCounts.getCounter(source).keySet()) {
				alignmentProbs.setCount(source, target, uniformInitValue);
			}
		}
		
		//Initialize distortion
		//dProbs[0] = 1.0;
		for(int i = 0; i <= MAX_SENTENCE_LENGTH+1; i++) {
			dProbs[i] = 1.0 / (double) (i+1);
			dCounts[i] = 0.0;
		}
		
		//Continue until convergence
		for (int i = 0; i < 100; i++) {
			//Set count(s|t) to 0 for all s,t and total(t) to 0 for all t
			for (String source : alignmentCounts.keySet()) {
				for (String target : alignmentCounts.getCounter(source).keySet()) {
					alignmentCounts.setCount(source, target, 0.0);
					targetWordCounts.setCount(target, 0.0);
				}
			}

			for (SentencePair pair : trainingPairs) {
				List<String> targetWords = pair.getEnglishWords();
				List<String> sourceWords = pair.getFrenchWords();
				int targetLen = targetWords.size();
				int sourceLen = sourceWords.size();
				Counter<String> sourceWordCounts = new Counter<String>();
				int sourcePos = 0;
				for (String source : sourceWords) {
					sourcePos++;
					sourceWordCounts.setCount(source, 0.0);
					int targetPos = 0;
					for (String target : targetWords) {
						targetPos++;
						int b = bucket(targetPos, sourcePos, targetLen, sourceLen);
						double bval = dProbs[b];
						double pairProb = alignmentProbs.getCount(source, target) * bval;
						sourceWordCounts.incrementCount(source, pairProb);
					}
					double bval = dProbs[MAX_SENTENCE_LENGTH+1];
					double nullPairProb = alignmentProbs.getCount(source, NULL_WORD) * bval;
					sourceWordCounts.incrementCount(source, nullPairProb);
					double curSourceCount = sourceWordCounts.getCount(source);
					targetPos = 0;
					for (String target : targetWords) {
						targetPos++;
						double curAlignmentProb = alignmentProbs.getCount(source, target);
						int b = bucket(targetPos, sourcePos, targetLen, sourceLen);
						bval = dProbs[b];
						double incVal = curAlignmentProb / curSourceCount * bval;
						alignmentCounts.incrementCount(source, target, incVal);
						targetWordCounts.incrementCount(target, incVal);
						dCounts[b] += incVal;
						totalD += incVal;
					}
					double nullAlignmentProb = alignmentProbs.getCount(source, NULL_WORD);
					bval = dProbs[MAX_SENTENCE_LENGTH+1];
					double incVal = nullAlignmentProb / curSourceCount * bval;
					alignmentCounts.incrementCount(source, NULL_WORD, incVal);
					targetWordCounts.incrementCount(NULL_WORD, incVal);
					dCounts[MAX_SENTENCE_LENGTH+1] += incVal;
					totalD += incVal;
				}
			}
			System.out.print("");
			//Compute new word pair probabilities
			for (String source : alignmentCounts.keySet()) {
				for (String target : alignmentCounts.getCounter(source).keySet()) {
					double curAlignmentCount = alignmentCounts.getCount(source, target);
					double curTargetTotal = targetWordCounts.getCount(target);
					double thing = curAlignmentCount / curTargetTotal;
					alignmentProbs.setCount(source, target, curAlignmentCount / curTargetTotal);
					//System.out.println("PROBS: " + source + " " + target + " " + thing);
				}
			}
			
			//Smooth distortion counts
			double l = 1.0;
			for(int a = 0; a <= MAX_SENTENCE_LENGTH+1; a++) {
				double e1 = dCounts[a];
				if(e1 > 0 && e1 < l) {
					l = e1;
				}
			}
			l = 0.5 * l;
			for(int a = 0; a <= MAX_SENTENCE_LENGTH+1; a++) {
				dCounts[a] += l;
			}
			totalD += l*(MAX_SENTENCE_LENGTH+1);
			
			//Compute new distortion probabilities
			for(int a = 0; a <= MAX_SENTENCE_LENGTH+1; a++) {
				dProbs[a] = dCounts[a] / totalD;
				dCounts[a] = 0.0;
			}
			totalD = 0.0;
		}
	}	
}
