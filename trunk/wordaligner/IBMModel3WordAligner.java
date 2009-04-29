package cs224n.wordaligner;  

import cs224n.util.*;
import java.util.*;

/**
 * Simple alignment baseline which maps french positions to english positions.
 * If the french sentence is longer, all final word map to null.
 */
public class IBMModel3WordAligner extends WordAligner {
  
  private CounterMap<String,String> alignmentProbs, alignmentCounts;
  private CounterMap<String, Integer> fertilityProbs, fertilityCounts;
  //private double[][][][] dProbs, dCounts;
  private double[] dProbs, dCounts;
  public static int MAX_SENTENCE_LENGTH = 3;
  public static int MAX_FERTILITY = 5;
  public double NULL_PROB = 0.2;
  
  public Alignment alignSentencePair(SentencePair sentencePair) {
    return model2AlignSentencePair(sentencePair);
  }

  private Alignment model2AlignSentencePair(SentencePair sentencePair) {
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
    return 0.0;

  }

  private double model2GetAlignmentProb(List<String> targetSentence, List<String> sourceSentence, Alignment alignment) {
    //Can be thought of as returning P(S, A|T) where T = target, S = source and A is the alignment
    //Simplifies to (sigma/(length_f +1)^(length_e)) * (product for each word in alignment->) * P(S_i | T_ai)
    double normalFactor = 1 / Math.pow(targetSentence.size() + 1, sourceSentence.size());
    double pairProduct = 1.0;
    for (int i = 0; i < sourceSentence.size(); i++) {
      String sourceWord = sourceSentence.get(i);
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
    trainModel2(trainingPairs);
    trainModel3(trainingPairs);
  }

  private void trainModel3(List<SentencePair> trainingPairs) {
    fertilityCounts = new CounterMap<String, Integer>();
    fertilityProbs = new CounterMap<String, Integer>();
    //Initialize fertility counts and pi count
    for (SentencePair pair : trainingPairs) {
      List<String> targetWords = pair.getEnglishWords();
      List<String> sourceWords = pair.getFrenchWords();
      double alpha[][] = new double[targetWords.size()][MAX_FERTILITY];
      for (int i  = 0; i < targetWords.size(); i++) {
	for (int k = 1; k < MAX_FERTILITY; k++) {
	  double beta = 0.0;
	  for (int j = 0; j < sourceWords.size(); j++) {
	    double aProbCount = alignmentProbs.getCount(sourceWords.get(j), targetWords.get(i));
	    beta += Math.pow(aProbCount / (1-aProbCount),k);
	  }
	  alpha[i][k] = Math.pow(-1, k+1) * beta / k;
	}
      }
      for (int i = 0; i < targetWords.size(); i++) {
	double r = 1.0;
	for (int j = 0; j < sourceWords.size(); j++) {
	  r *= (1 - alignmentProbs.getCount(sourceWords.get(j), targetWords.get(i)));
	}
	for (int phi = 0; phi < MAX_FERTILITY; phi++) {
	  double sum = 0.0;
	  ArrayList<ArrayList<Integer>> partitions = createPartitions(phi);
	  for (ArrayList<Integer> curPartition : partitions) {
	    double prod = 1.0;
	    for (int k = 0; k < phi; k++) {
	      int timesC = times(curPartition, k);
	      prod *= Math.pow(alpha[i][k], timesC) / factorial(timesC);
	    }
	    sum += prod;
	  }
	  fertilityCounts.incrementCount(targetWords.get(i), phi, r*sum);
	}
      }
    }

    for (int count = 0; count < 50; count++) {
    for (SentencePair pair : trainingPairs) {
      List<String> targetWords = pair.getEnglishWords();
      List<String> sourceWords = pair.getFrenchWords();
      Alignment viterbiAlign = model2AlignSentencePair(pair);
      Alignment bestAlign = viterbiAlign;
      double maxScore = getAlignmentProb(targetWords, sourceWords, bestAlign);
      while (true) {
	boolean hadBetter = false;
	for (int fr = 0; fr < sourceWords.size(); fr++) {
	  int initEnPos = viterbiAlign.getAlignedTarget(fr);
	  int initFrPos = fr;
	  
	  for (int d1 = -1; d1 <= 1; d1++) {
	    for (int d2 = -1; d2 <= 1; d2++) {
		if (d1 == 0 && d2 == 0) continue;
		int newEnPos = initEnPos+d1;
		int newFrPos = initFrPos+d2;
		if (newEnPos >= -1 && newEnPos <= (sourceWords.size() - 1) && newFrPos >= -1 && newFrPos <= (sourceWords.size() - 1)) {
		  viterbiAlign.removeAlignment(initEnPos, initFrPos);
		  viterbiAlign.addAlignment(newEnPos, newFrPos, true);
		  double curScore = getAlignmentProb(targetWords, sourceWords, viterbiAlign);
		  if (curScore > maxScore) {
		    hadBetter = true;
		    maxScore = curScore;
		    bestAlign = viterbiAlign;
		  }	
		  else {
		    viterbiAlign.removeAlignment(newEnPos, newFrPos);
		    viterbiAlign.addAlignment(initEnPos, initFrPos, true);
		  }
		}
	    }
	  }
	}
	if (!hadBetter) break;
      }

      for (int fr = 0; fr < sourceWords.size(); fr++) {
	int initEnPos = viterbiAlign.getAlignedTarget(fr);
	int initFrPos = fr;
	  
	for (int d1 = -2; d1 <= 2; d1++) {
	  for (int d2 = -2; d2 <= 2; d2++) {
	    if (d1 == 0 && d2 == 0) continue;
	    int newEnPos = initEnPos+d1;
	    int newFrPos = initFrPos+d2;
	    if (newEnPos >= -1 && newEnPos <= (sourceWords.size() - 1) && newFrPos >= -1 && newFrPos <= (sourceWords.size() - 1)) {
	      bestAlign.removeAlignment(initEnPos, initFrPos);
	      bestAlign.addAlignment(newEnPos, newFrPos, true);
	      int numNullWords = 0;
	      int[] freqArr = new int[targetWords.size()];
	      for (int a = 0; a < sourceWords.size(); a++) {
		int curAlignTarget = bestAlign.getAlignedTarget(a);
		System.out.println("source size " + sourceWords.size() + " targ size " + targetWords.size());
		System.out.println(a + " aligned with " + curAlignTarget);
		if(curAlignTarget == -1) {
		  numNullWords++;
		}
		else {
		  if(freqArr[curAlignTarget] < 5) freqArr[curAlignTarget]++;
		}
	      }
	      for (int fInd = 0; fInd < freqArr.length; fInd++) {
		fertilityCounts.incrementCount(targetWords.get(fInd), freqArr[fInd], 1.0);
	      }
	      NULL_PROB += (numNullWords/sourceWords.size() - 0.2) / trainingPairs.size();
	      bestAlign.removeAlignment(newEnPos, newFrPos);
	      bestAlign.addAlignment(initEnPos, initFrPos, true);
	    }
	  }
	}			
      }
    }
    }
  }
  
  private ArrayList<ArrayList<Integer>> createPartitions(int phi) {
    ArrayList<ArrayList<Integer>> partitions = new ArrayList<ArrayList<Integer>>();
    createPartitionsRecur(phi, 0, new ArrayList<Integer>(), partitions);
    return partitions;
  }

  private void createPartitionsRecur(int phi, int start, ArrayList<Integer> curPart, ArrayList<ArrayList<Integer>> partitions) {
    if (phi <= 0) {
      for (int i = 0; i < partitions.size(); i++) {
	if (partitions.get(i).size() == curPart.size()) return;
      }
      partitions.add(new ArrayList<Integer>(curPart));
    }

    for (int i = 1; i <= phi; i++) {
      curPart.add(start, new Integer(i));
      createPartitionsRecur(phi - i, start + 1, curPart, partitions);
      curPart.remove(start);
    }
  }

  private int times(ArrayList<Integer> partition, int num) {
    int count = 0;
    for (int i = 0; i < partition.size(); i++) {
      if(partition.get(i).equals(new Integer(i)))
	count++;
    }
    return count;
  }
  
  private int factorial(int input) {
    int result = 1;
    for (int i = 2; i <= input; i++) {
      result *= input;
    }
    return result;
  }

  private void trainModel2(List<SentencePair> trainingPairs) {
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
