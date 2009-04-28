package cs224n.wordaligner;  

import cs224n.util.*;
import java.util.List;

/**
   * Simple alignment baseline which maps french positions to english positions.
   * If the french sentence is longer, all final word map to null.
   */
  public class IBMModel1WordAligner extends WordAligner {

    private CounterMap<String,String> alignmentProbs, alignmentCounts;

    public Alignment alignSentencePair(SentencePair sentencePair) {
      Alignment alignment = new Alignment();
      List<String> sourceWords = sentencePair.getFrenchWords();
      List<String> targetWords = sentencePair.getEnglishWords();
      for (int i = 0; i < sourceWords.size(); i++) {
	String curSourceWord = sourceWords.get(i);
	int maxAlignPos = 0;
	double maxAlignProb = alignmentProbs.getCount(curSourceWord, targetWords.get(0));
	for (int j = 1; j < targetWords.size(); j++) {
	  double curAlignProb = alignmentProbs.getCount(curSourceWord, targetWords.get(j));
	  if (curAlignProb > maxAlignProb) {
	    maxAlignPos = j;
	    maxAlignProb = curAlignProb;
	  }
	}
	//Deal with NULL_WORD case
	if (alignmentProbs.getCount(curSourceWord, NULL_WORD) > maxAlignProb)
	  maxAlignPos = -1;
	alignment.addAlignment(maxAlignPos, i, true);
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
	int curAlignTarget = alignment.getAlignedTarget(i);
	String targetWord = (curAlignTarget == -1 ? NULL_WORD : targetSentence.get(curAlignTarget));
	pairProduct *= alignmentProbs.getCount(sourceWord, targetWord);
      }
      double thingy = normalFactor * pairProduct;
      return normalFactor * pairProduct;
    }

    
    public CounterMap<String,String> getProbSourceGivenTarget(){
      //Return the counter map storing P(s|t) for each french/english word pair
      return alignmentProbs;
    }

    public void train(List<SentencePair> trainingPairs) {
      alignmentProbs = new CounterMap<String,String>();
      alignmentCounts = new CounterMap<String,String>();

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

      //Continue until convergence
      for (int i = 0; i < 500; i++) {
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
	  Counter<String> sourceWordCounts = new Counter<String>();
	  for (String source : sourceWords) {
	    sourceWordCounts.setCount(source, 0.0);
	    for (String target : targetWords) {
	      double pairProb = alignmentProbs.getCount(source, target);
	      sourceWordCounts.incrementCount(source, pairProb);
	    }
	    //Populate NULL_WORD entries
	    double nullPairProb = alignmentProbs.getCount(source, NULL_WORD);
	    sourceWordCounts.incrementCount(source, nullPairProb);
	  }
	  for (String source : sourceWords) {
	    double curSourceCount = sourceWordCounts.getCount(source);
	    for (String target : targetWords) {
	      double curAlignmentProb = alignmentProbs.getCount(source, target);
	      alignmentCounts.incrementCount(source, target, (curAlignmentProb / curSourceCount));
	      targetWordCounts.incrementCount(target, (curAlignmentProb / curSourceCount));
	    }
	    //Populate NULL_WORD entries
	    double nullAlignmentProb = alignmentProbs.getCount(source, NULL_WORD);
	    alignmentCounts.incrementCount(source, NULL_WORD, (nullAlignmentProb / curSourceCount));
	    targetWordCounts.incrementCount(NULL_WORD, (nullAlignmentProb / curSourceCount));
	  }
	}

	//Compute new word pair probabilities
	for (String source : alignmentCounts.keySet()) {
	  for (String target : alignmentCounts.getCounter(source).keySet()) {
	    double curAlignmentCount = alignmentCounts.getCount(source, target);
	    double curTargetTotal = targetWordCounts.getCount(target);
	    double thing = curAlignmentCount / curTargetTotal;
	    alignmentProbs.setCount(source, target, curAlignmentCount / curTargetTotal);
	  }
	}
      }	
    }
  }
