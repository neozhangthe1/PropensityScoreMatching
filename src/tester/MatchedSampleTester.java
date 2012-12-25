package tester;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;


import weka.classifiers.functions.Logistic;
import weka.core.Attribute;
import weka.core.Instances;
import weka.core.SparseInstance;

public class MatchedSampleTester {
	
	Instances dataset = null;
	ArrayList<Subject> subjects = new ArrayList<Subject>();
	int treatedNum;
	int untreatedNum;
	public MatchedSampleTester() {
		treatedNum = 0;
		untreatedNum = 0;
	}
	
	
	public void AddSubject(ArrayList<Double> features, boolean adopted, boolean treated) {
		if (dataset == null) {
			ArrayList<Attribute> attrv = new ArrayList<Attribute>();
			for (int i = 0; i < features.size(); ++i) attrv.add(new Attribute(Integer.toString(i)));
			ArrayList<String> classList = new ArrayList<String>();
			classList.add("untreated");
			classList.add("treated");
			attrv.add(new Attribute("class", classList));
			dataset = new Instances("set", attrv, 100);
			dataset.setClassIndex(features.size());
		}
		SparseInstance ins = new SparseInstance(features.size() + 1);
		ins.setDataset(dataset);
		for (int i = 0; i < features.size(); ++i) {
			ins.setValue(i, features.get(i));
		}
		if (treated == true) {
			ins.setValue(features.size(), "treated");
			++treatedNum;
		}
		else {
			ins.setValue(features.size(), "untreated");
			++untreatedNum;
		}
		//dataset.add(ins);
		Subject subject = new Subject(adopted, treated, ins);
		subjects.add(subject);
	}
	
	public void MatchedSampleTest(ArrayList<Integer> treatedResult, ArrayList<Integer> untreatedResult) {
		long ts0, ts1;
		//if (treatedNum < 100) return;
		treatedResult.clear();
		untreatedResult.clear();
		
		double r = 0.0;
		ts0 = System.currentTimeMillis();
		
		for (int i = 0; i < subjects.size(); ++i) {
			Subject sub = subjects.get(i);
			if (sub.treated) {
				dataset.add(sub.instance);
			}
			else {
				//r += 1.0;
				r += (double)  treatedNum / untreatedNum;
				if (r >= 1.0) {
					dataset.add(sub.instance);
					r = 0.0;
				}
			}
		}
		
		if (dataset.numInstances() < 2) {
			return;
		}
		
		ts1 = System.currentTimeMillis();
		System.out.println("Sample Timecost: " + (ts1-ts0) + "ms.");
		
		System.out.println("treated:" + treatedNum + ",  untreated:" + untreatedNum);
		System.out.println("dataset size:" + dataset.numInstances());
		ts0 = System.currentTimeMillis();

		try{
			Logistic classifier = new Logistic();
			classifier.buildClassifier(dataset);
			@SuppressWarnings("unused")
			int tn = 0, tp = 0, fn = 0, fp = 0;
			for (Subject sub : subjects) {
				sub.propensity = classifier.distributionForInstance(sub.instance)[1];
				if (sub.propensity < 0.5) {
					if (sub.treated == false){
						++tn;
					}
					else{
						++fn;
					}
				}
				else {
					if (sub.treated == false) ++fp;
					else {
						++tp;
					}
				}
			}
			double prec = (double)tp/(tp+fp);
			double rec = (double)tp/(tp+fn);
			double f1 = 2*prec*rec/(prec+rec);
			System.out.println("Prec:" + prec);
			System.out.println("Rec: " + rec);
			System.out.println("F1:  " + f1);
			System.out.println("Coefficents:");
			for (int i = 0; i < 1; ++i) {
				for (int j = 0; j < dataset.numAttributes() - 1; ++j) {
					System.out.print(classifier.coefficients()[j][i] + "\t");
				}
				System.out.println();
			}
		} catch(Exception e) {
			e.printStackTrace();
			return;
		}
		Collections.sort(subjects, new Comparator<Subject>() {
			public int compare(Subject s1, Subject s2) {
				return Double.compare(s1.propensity, s2.propensity);
			} 
		});
		ts1 = System.currentTimeMillis();
		System.out.println("Calculate Propensity Timecost: " + (ts1-ts0) + "ms.");
		ts0 = System.currentTimeMillis();
		Subject[] treatedList = new Subject[treatedNum];
		Subject[] untreatedList = new Subject[untreatedNum];
		boolean[] untreatedMark = new boolean[untreatedNum];
		int t1 = 0, t2 = 0;
		for (int i = 0; i < subjects.size(); ++i) {
			if (subjects.get(i).treated == true) {
				treatedList[t1++] = subjects.get(i);
			}
			else {
				untreatedMark[t2] = false;
				untreatedList[t2++] = subjects.get(i);
			}
		}
		int[] treatedCnt = new int[2];
		int[] untreatedCnt = new int[2];
		treatedCnt[0] = treatedCnt[1] = untreatedCnt[0] = untreatedCnt[1] = 0;
		int p2 = 0;
		int[] matchedId = new int[treatedList.length];
		double diffSum = 0.0;
		double diffSquareSum = 0.0;
		int N = 0;
		for (int p1 = 0; p1 < treatedList.length; ++p1) {
			double q1 = treatedList[p1].propensity;
			double q2;
			while (p2 < untreatedList.length && (untreatedList[p2].propensity < q1 || untreatedMark[p2] == true)) ++p2;
			int tmpp2;
			if (p2 >= untreatedList.length) tmpp2 = untreatedList.length - 1;
			else tmpp2 = p2;
			while (tmpp2 >= 0 && (untreatedList[tmpp2].propensity >= q1 || untreatedMark[tmpp2] == true)) --tmpp2;
			int fp2;
			if (p2 >= untreatedList.length) {
				if (tmpp2 < 0) {
					p2 = tmpp2 = 0;
					continue;
				}
				else {
					fp2 = tmpp2;
				}
			}
			else {
				if (tmpp2 < 0) {
					tmpp2 = 0;
					fp2 = p2;
				}
				else {
					if (Math.abs(untreatedList[tmpp2].propensity - q1) < Math.abs(untreatedList[p2].propensity - q1)) {
						fp2 = tmpp2;
					}
					else {
						fp2 = p2;
					}
				}
			}
			q2 = untreatedList[fp2].propensity;
			double diff = Math.abs(q2 - q1);
			diffSum += diff;
			diffSquareSum += diff*diff;
			matchedId[p1] = fp2;
			if (untreatedMark[fp2] == true) {
				System.err.println("MARKED!");
				matchedId[p1] = -1;
				continue;
			}
			++N;
			untreatedMark[fp2] = true;
			p2 = tmpp2;
		}
		double diffSigma = Math.sqrt((double)diffSquareSum/N - ((double)diffSum/N) * ((double)diffSum/N));
		System.out.println("DiffSigma: " + diffSigma);
		for (int p1 = 0; p1 < treatedList.length; ++p1) {
			if (matchedId[p1] == -1) continue;
			Subject s1 = treatedList[p1];
			Subject s2 = untreatedList[matchedId[p1]];
			if (Math.abs(s1.propensity - s2.propensity) >= 2 * diffSigma) continue;
			if (s1.adopted == false) ++treatedCnt[0];
			else ++treatedCnt[1];
			if (s2.adopted == false) ++untreatedCnt[0];
			else ++untreatedCnt[1];
			/*
			if (s2.adopted && !s1.adopted) 
			{
				System.out.println("Treated:  " + s1.instance + " propensity:" + s1.propensity);
				System.out.println("Untreated:" + s2.instance + " propensity:" + s2.propensity);
			}*/
		}
		
		ts1 = System.currentTimeMillis();
		System.out.println("Matching Timecost: " + (ts1-ts0) + "ms.");
		treatedResult.add(treatedCnt[0]);
		treatedResult.add(treatedCnt[1]);
		untreatedResult.add(untreatedCnt[0]);
		untreatedResult.add(untreatedCnt[1]);
	}
}
