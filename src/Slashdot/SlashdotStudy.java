package Slashdot;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.sql.Date;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import tester.MatchedSampleTester;
import weka.classifiers.functions.Logistic;
import weka.core.Attribute;
import weka.core.Instances;
import weka.core.SparseInstance;


import Network.NetworkData;
import Network.NetworkEdge;




public class SlashdotStudy {
	NetworkData nwdata;
	SlashdotCommentData cmdata;
	final long MAXTIME = Long.MAX_VALUE;
	
	
	private long[] checkAdoptedVector;
	
	
	public SlashdotStudy(String nwdataFileName, String cmdataFileName) {
		nwdata = new NetworkData();
		nwdata.ReadNetwork(nwdataFileName);
		cmdata = new SlashdotCommentData(nwdata);
		cmdata.readSlashdotComments(cmdataFileName);
		checkAdoptedVector = new long[nwdata.getNodeNum()];
		System.out.println("Nodenum:" + nwdata.getNodeNum());
		System.out.println("Edgenum:" + nwdata.getEdgeNum());
		System.out.println("Newsnum:" + cmdata.getNewsNum());
	}
	
	public void UpdateCheckAdoptedVector(int newsid) {
		for (int i = 0; i < nwdata.getNodeNum(); ++i) {
			checkAdoptedVector[i] = CheckAdopted(newsid, i);
		}
	}
	
	public long CheckAdopted(int newsid, int nodeid) {
		if (!cmdata.userComment.get(nodeid).containsKey(newsid)) return MAXTIME;
		else {
			return cmdata.userComment.get(nodeid).get(newsid).get(0).t.getTime();
		}
	}
	
	public void StudyNegativeCorrByAssortativeMixing(String outputFileName) {
		final int MAXFOE = 200;
		int nnum = cmdata.getNewsNum();
		int nodenum = nwdata.getNodeNum();
		int[][][] stat = new int[2][nnum][MAXFOE];
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < nnum; ++j) {
				for (int k = 0; k < MAXFOE; ++k) {
					stat[i][j][k] = 0;
				}
			}
		}
		for (int nid = 0; nid < nnum; ++nid) {
			for (int id = 0; id < nodenum; ++id) {
				//Check Adopted
				int adopted;
				long adtime = CheckAdopted(nid, id);
				if (adtime == MAXTIME) {
					adopted = 0;
				}
				else {
					adopted = 1;
				}
				
				//Check # of Adopter Foes
				int adfoes = 0;
				ArrayList<NetworkEdge> elist = nwdata.edgeFromList.get(id);
				for (NetworkEdge e : elist) {
					if (e.weight > 0) continue;
					int v2id = e.v2id;
					long check = CheckAdopted(nid, v2id);
					if (check >= adtime) continue;
					++adfoes;
				}
				
				stat[adopted][nid][adfoes]++;
			}
			if (nid % 10 == 0 || nid == nnum - 1) {
				System.out.println("Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
			}
		}

		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			for (int k = 0; k < MAXFOE; ++k) {
				long macropad = 0, macrononad = 0;
				for (int i = 0; i < nnum; ++i) {
					macropad += stat[1][i][k];
					macrononad += stat[0][i][k];
				}
				bw.write(k + "\t" + macropad + "\t" + macrononad + "\t" + (macropad + macrononad) + "\n");
			}
			bw.flush();
			bw.close();
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void StudyCorrByAssortativeMixing(String outputFileName) {
		final int MAXFOE = 200;
		final int MAXFRI = 200;
		int nnum = cmdata.getNewsNum();
		int nodenum = nwdata.getNodeNum();
		int[][][][] stat = new int[2][nnum][MAXFOE][MAXFRI];
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < nnum; ++j) {
				for (int k = 0; k < MAXFOE; ++k) {
					for (int l = 0; l < MAXFRI; ++l) {
						stat[i][j][k][l] = 0;
					}
				}
			}
		}
		for (int nid = 0; nid < nnum; ++nid) {
			for (int id = 0; id < nodenum; ++id) {
				//Check Adopted
				int adopted;
				long adtime = CheckAdopted(nid, id);
				if (adtime == MAXTIME) {
					adopted = 0;
				}
				else {
					adopted = 1;
				}
				
				//Check # of Adopter Foes
				int adfoes = 0;
				int adfris = 0;
				ArrayList<NetworkEdge> elist = nwdata.edgeFromList.get(id);
				for (NetworkEdge e : elist) {
					int v2id = e.v2id;
					long check = CheckAdopted(nid, v2id);
					if (check >= adtime) continue;
					if (e.weight > 0) {
						++adfris;
					}
					else {
						++adfoes;
					}
				}
				
				stat[adopted][nid][adfoes][adfris]++;
			}
			if (nid % 10 == 0 || nid == nnum - 1) {
				System.out.println("Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
			}
		}

		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			for (int k = 0; k < MAXFOE; ++k) {
				for (int l = 0; l < MAXFRI; ++l) {
					long macropad = 0, macrononad = 0;
					for (int i = 0; i < nnum; ++i) {
						macropad += stat[1][i][k][l];
						macrononad += stat[0][i][k][l];
					}
					bw.write(k + "\t" + macropad + "\t" + macrononad + "\t" + (macropad + macrononad) + "\t");
				}
				bw.write("\n");
			}
			bw.flush();
			bw.close();
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	
	
	
	private SparseInstance makePositiveInstance(int userId, int newsId, Instances dataset, int threshold) {
		SparseInstance ins = new SparseInstance(5);
		ins.setDataset(dataset);

		int attrId = 0;
		int fri = 0, foe = 0, fan = 0, frk = 0;
		int adfri = 0, adfoe = 0, adfan = 0, adfrk = 0;
		ArrayList<NetworkEdge> elist;
		long adtime = checkAdoptedVector[userId];
		elist = nwdata.edgeFromList.get(userId);
		for (NetworkEdge e: elist) {
			if (e.weight > 0) {
				if (checkAdoptedVector[e.v2id] < adtime) ++adfri;
				++fri;
			}
			else {
				if (checkAdoptedVector[e.v2id] < adtime) ++adfoe;
				++foe;
			}
		}
		

		elist = nwdata.edgeToList.get(userId);
		for (NetworkEdge e: elist) {
			if (e.weight > 0) {
				if (checkAdoptedVector[e.v1id] < adtime) ++adfan;
				++fan;
			}
			else {
				if (checkAdoptedVector[e.v1id] < adtime) ++adfrk;
				++frk;
			}
		}		
		//ins.setValue(attrId++, adfri);
		//ins.setValue(attrId++, adfoe);
		//ins.setValue(attrId++, adfan);
		//ins.setValue(attrId++, adfrk);


		ins.setValue(attrId++, fri);
		ins.setValue(attrId++, foe);
		ins.setValue(attrId++, fan);
		ins.setValue(attrId++, frk);

		//ins.setValue(attrId++, adfoe);

		
		if (adfri >= threshold) {
			ins.setValue(attrId++, "YES");
			return ins;
		}
		else if (adfri == 0){
			ins.setValue(attrId++, "NO");
			return ins;
		}
		else return null;
		
		
	}

	private SparseInstance makeInstance(int userId, int newsId, Instances dataset, int threshold) {
		//System.out.println("Making instance...");
		SparseInstance ins = new SparseInstance(6);
		ins.setDataset(dataset);
		
		int attrId = 0;
		int pos, neg;
		int adfri = 0;
		int adfoe = 0;
		ArrayList<NetworkEdge> elist;
		long adtime = checkAdoptedVector[userId];
		pos = 0; neg = 0;
		elist = nwdata.edgeFromList.get(userId);
		for (NetworkEdge e: elist) {
			if (e.weight > 0) {
				if (checkAdoptedVector[e.v2id] < adtime) ++adfri;
				++pos;
			}
			else {
				if (checkAdoptedVector[e.v2id] < adtime) ++adfoe;
				++neg;
			}
		}		
		ins.setValue(attrId++, pos);
		ins.setValue(attrId++, neg);
		
		pos = 0; neg = 0;
		elist = nwdata.edgeToList.get(userId);
		for (NetworkEdge e: elist) {
			if (e.weight > 0) ++pos;
			else ++neg;
		}		
		ins.setValue(attrId++, pos);
		ins.setValue(attrId++, neg);
		ins.setValue(attrId++, adfri);
		
		if (adfoe >= threshold) {
			ins.setValue(attrId++, "YES");
			return ins;
		}
		else if (adfoe == 0){
			ins.setValue(attrId++, "NO");
			return ins;
		}
		else return null;
		
		
	}
	
	public void StudyCorrByAssortativeMixingMatchedSample(String outputFileName, int threshold) {
		
		//Extract attributes and add to dataset.
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			System.out.println("Adding...");
			int nnum = cmdata.getNewsNum();
			//ArrayList<SlashdotUserNewsPair> ctrlGroup = new ArrayList<SlashdotUserNewsPair>();
			//ArrayList<SlashdotUserNewsPair> testGroup = new ArrayList<SlashdotUserNewsPair>();
			long[] timecounter = new long[6];
			for (int nid = 0; nid < nnum; ++nid) {
				System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
				long t0 = System.currentTimeMillis();
				ArrayList<Attribute> attrv = new ArrayList<Attribute>();
				attrv.add(new Attribute("friendsnum"));
				attrv.add(new Attribute("foesnum"));
				attrv.add(new Attribute("fansnum"));
				attrv.add(new Attribute("freaksnum"));
				attrv.add(new Attribute("adoptedfoe"));
				ArrayList<String> classList = new ArrayList<String>();
				classList.add("YES");
				classList.add("NO");
				attrv.add(new Attribute("class", classList));
				Instances dataset = new Instances("trainset", attrv, nwdata.getNodeNum());
				dataset.setClassIndex(dataset.numAttributes() - 1);
				long t1 = System.currentTimeMillis();
				UpdateCheckAdoptedVector(nid);
				long t2 = System.currentTimeMillis();
				ArrayList<SlashdotUserNewsPair> userNewsPairList = new ArrayList<SlashdotUserNewsPair>();
				int posnum = 0, negnum = 0;
				for (int uid = 0; uid < nwdata.getNodeNum(); ++uid) {
					SlashdotUserNewsPair userNewsPair = new SlashdotUserNewsPair();
					//long adptime = CheckAdopted(nid, uid);
					long adptime = checkAdoptedVector[uid];
					if (adptime == MAXTIME) {
						userNewsPair.adopted = 0;
					}
					else {
						userNewsPair.adopted = 1;
					}
					
					SparseInstance ins = makeInstance(uid, nid, dataset, threshold);
					if (ins == null) continue;
					//dataset.add(ins);
					//inslist.add(ins);
					userNewsPair.ins = ins;
					userNewsPair.treated = (int) (1 - ins.classValue());
					if (userNewsPair.treated == 1) {
						++posnum;
					}
					else {
						++negnum;
					}
					userNewsPairList.add(userNewsPair);
				}
				double r = 0.0;
				for (int j = 0; j < userNewsPairList.size(); ++j) {
					if (userNewsPairList.get(j).treated == 1) {
						dataset.add(userNewsPairList.get(j).ins);
					}
					else {
						r += (double) posnum / negnum;
						if (r >= 1) {
							dataset.add(userNewsPairList.get(j).ins);
							r = 0;
						}
					}
				}
				long t3 = System.currentTimeMillis();
				double propensitySum = 0.0;
				double propensitySquareSum = 0.0;
				try {
					Logistic classifier = new Logistic();
					classifier.buildClassifier(dataset);
					
					for (SlashdotUserNewsPair pair : userNewsPairList) {
						double[] res = classifier.distributionForInstance(pair.ins);
						pair.treatedPropensity = res[0];
						propensitySum += pair.treatedPropensity;
						propensitySquareSum += pair.treatedPropensity * pair.treatedPropensity;
					}
				}catch(Exception e) {
					e.printStackTrace();
					continue;
				}
				double propensitySigma = Math.sqrt((double)propensitySquareSum/userNewsPairList.size() - ((double)propensitySum/userNewsPairList.size()) * ((double)propensitySum/userNewsPairList.size()));
				long t4 = System.currentTimeMillis();
				Collections.sort(userNewsPairList, new Comparator<SlashdotUserNewsPair>() {
					public int compare(SlashdotUserNewsPair c1, SlashdotUserNewsPair c2) {
						return Double.compare(c1.treatedPropensity, c2.treatedPropensity);
					} 
				});
				
				ArrayList<SlashdotUserNewsPair> treatedList = new ArrayList<SlashdotUserNewsPair>();
				ArrayList<SlashdotUserNewsPair> untreatedList = new ArrayList<SlashdotUserNewsPair>();
				
				for (int p = 0; p < userNewsPairList.size(); ++p) {
					if (userNewsPairList.get(p).treated == 1) {
						treatedList.add(userNewsPairList.get(p));
					}
					else {
						untreatedList.add(userNewsPairList.get(p));
					}
				}
				boolean[] treatedMark = new boolean[treatedList.size()];
				boolean[] untreatedMark = new boolean[untreatedList.size()];
				int p2 = 0;
				int[] treatedCnt = new int[2];
				int[] untreatedCnt = new int[2];

				for (int p1 = 0; p1 < treatedList.size(); ++p1) {
					double q1 = treatedList.get(p1).treatedPropensity;
					double q2;
					while (p2 < untreatedList.size() && untreatedList.get(p2).treatedPropensity < q1) ++p2;
					if (p2 < untreatedList.size()) q2 = untreatedList.get(p2).treatedPropensity;
					else break;
					double diff = q2 - q1;
					if (p2 != 0 && untreatedMark[p2 - 1] == false) {
						if (diff > q1 - untreatedList.get(p2-1).treatedPropensity) {
							--p2;
							diff = q1 - untreatedList.get(p2).treatedPropensity;
						}
					}
					if (diff < propensitySigma * 2) {
						treatedMark[p1] = true;
						untreatedMark[p2] = true;
						SlashdotUserNewsPair pair1 = untreatedList.get(p2);
						SlashdotUserNewsPair pair2 = treatedList.get(p1);
						if (pair1.adopted == 0) ++untreatedCnt[0];
						else ++untreatedCnt[1];
						if (pair2.adopted == 0) ++treatedCnt[0];
						else ++treatedCnt[1];
						//ctrlGroup.add(untreatedList.get(p2));
						//testGroup.add(treatedList.get(p1));

						//============DEBUG==============
						/*
						SlashdotUserNewsPair pair1 = untreatedList.get(p2);
						SlashdotUserNewsPair pair2 = treatedList.get(p1);
						System.out.println(pair1.treatedPropensity + "  vs  " + pair2.treatedPropensity);
						System.out.println(pair1.ins);
						System.out.println(pair2.ins);
						System.out.println("===");
						*/
						//================================
						++p2;
					}
				}		
				long t5 = System.currentTimeMillis();
				bw.write(treatedCnt[1] + "\t" + treatedCnt[0] + "\t");
				bw.write(untreatedCnt[1] + "\t" + untreatedCnt[0] + "\t");
				bw.write("\n");
				bw.flush();
				long t6 = System.currentTimeMillis();
				timecounter[0] += t1 - t0;
				timecounter[1] += t2 - t1;
				timecounter[2] += t3 - t2;
				timecounter[3] += t4 - t3;
				timecounter[4] += t5 - t4;
				timecounter[5] += t6 - t5;
				for (int i = 0; i < 6; ++i) {
					System.out.println("Timecost " + i + ":  " + timecounter[i] + " ms.");
				}
			}
			bw.close();

		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void StudyCorrByAssortativeMixingMatchedSamplePos(String outputFileName, int threshold) {
		
		//Extract attributes and add to dataset.
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			System.out.println("Adding...");
			int nnum = cmdata.getNewsNum();
			long[] timecounter = new long[6];
			for (int nid = 0; nid < nnum; ++nid) {
				long t0 = System.currentTimeMillis();
				ArrayList<Attribute> attrv = new ArrayList<Attribute>();
				//attrv.add(new Attribute("adfoe"));
				//attrv.add(new Attribute("adfan"));
				//attrv.add(new Attribute("adfrk"));
				attrv.add(new Attribute("fri"));
				attrv.add(new Attribute("foe"));
				attrv.add(new Attribute("fan"));
				attrv.add(new Attribute("frk"));
				//attrv.add(new Attribute("adoptedfri"));
				ArrayList<String> classList = new ArrayList<String>();
				classList.add("YES");
				classList.add("NO");
				attrv.add(new Attribute("class", classList));
				Instances dataset = new Instances("trainset", attrv, nwdata.getNodeNum());
				dataset.setClassIndex(dataset.numAttributes() - 1);
				long t1 = System.currentTimeMillis();
				UpdateCheckAdoptedVector(nid);
				long t2 = System.currentTimeMillis();
				ArrayList<SlashdotUserNewsPair> userNewsPairList = new ArrayList<SlashdotUserNewsPair>();
				int posnum = 0, negnum = 0;
				for (int uid = 0; uid < nwdata.getNodeNum(); ++uid) {
					SlashdotUserNewsPair userNewsPair = new SlashdotUserNewsPair();
					userNewsPair.uid = uid;
					userNewsPair.nid = nid;
					long adptime = checkAdoptedVector[uid];
					if (adptime == MAXTIME) {
						userNewsPair.adopted = 0;
					}
					else {
						userNewsPair.adopted = 1;
					}
					
					SparseInstance ins = makePositiveInstance(uid, nid, dataset, threshold);
					if (ins == null) continue;
					//dataset.add(ins);
					userNewsPair.ins = ins;
					userNewsPair.treated = (int) (1 - ins.classValue());
					if (userNewsPair.treated == 1) {
						++posnum;
					}
					else {
						++negnum;
					}
					userNewsPairList.add(userNewsPair);
				}
				
				double r = 0.0;
				for (int j = 0; j < userNewsPairList.size(); ++j) {
					if (userNewsPairList.get(j).treated == 1) {
						dataset.add(userNewsPairList.get(j).ins);
					}
					else {
						r += (double) posnum / negnum;
						if (r >= 1) {
							dataset.add(userNewsPairList.get(j).ins);
							r = 0;
						}
					}
				}
				long t3 = System.currentTimeMillis();
				double propensitySum = 0.0;
				double propensitySquareSum = 0.0;
				try {
					Logistic classifier = new Logistic();
					classifier.buildClassifier(dataset);
					
					for (SlashdotUserNewsPair pair : userNewsPairList) {
						double[] res = classifier.distributionForInstance(pair.ins);
						pair.treatedPropensity = res[0];
						propensitySum += pair.treatedPropensity;
						propensitySquareSum += pair.treatedPropensity * pair.treatedPropensity;
					}
				}catch(Exception e) {
					e.printStackTrace();
					continue;
				}
				double propensitySigma = Math.sqrt((double)propensitySquareSum/userNewsPairList.size() - ((double)propensitySum/userNewsPairList.size()) * ((double)propensitySum/userNewsPairList.size()));
				long t4 = System.currentTimeMillis();
				Collections.sort(userNewsPairList, new Comparator<SlashdotUserNewsPair>() {
					public int compare(SlashdotUserNewsPair c1, SlashdotUserNewsPair c2) {
						return Double.compare(c1.treatedPropensity, c2.treatedPropensity);
					} 
				});
				
				ArrayList<SlashdotUserNewsPair> treatedList = new ArrayList<SlashdotUserNewsPair>();
				ArrayList<SlashdotUserNewsPair> untreatedList = new ArrayList<SlashdotUserNewsPair>();
				
				for (int p = 0; p < userNewsPairList.size(); ++p) {
					if (userNewsPairList.get(p).treated == 1) {
						treatedList.add(userNewsPairList.get(p));
					}
					else {
						untreatedList.add(userNewsPairList.get(p));
					}
				}
				boolean[] treatedMark = new boolean[treatedList.size()];
				boolean[] untreatedMark = new boolean[untreatedList.size()];
				for (int i = 0; i < treatedList.size(); ++i) treatedMark[i] = false;
				for (int i = 0; i < untreatedList.size(); ++i) untreatedMark[i] = false;
				
				int[] treatedCnt = new int[2];
				int[] untreatedCnt = new int[2];
				
				/*
				for (int p1 = 0; p1 < treatedList.size(); ++p1) {
					double q1 = treatedList.get(p1).treatedPropensity;
					double mindiff = 1.0;
					int minp2 = -1;
					for (int p2 = 0; p2 < untreatedList.size(); ++p2) {
						if (untreatedMark[p2] == true) continue;
						if (Math.abs(untreatedList.get(p2).treatedPropensity - q1) < mindiff) {
							mindiff = Math.abs(untreatedList.get(p2).treatedPropensity - q1);
							minp2 = p2;
						}
					}
					if (minp2 == -1) break;
					untreatedMark[minp2] = true;
					treatedMark[p1] = true;
					SlashdotUserNewsPair pair1 = treatedList.get(p1);
					SlashdotUserNewsPair pair2 = untreatedList.get(minp2);
					if (pair1.adopted == 0) ++treatedCnt[0];
					else ++treatedCnt[1];
					if (pair2.adopted == 0) ++untreatedCnt[0];
					else ++untreatedCnt[1];
				}*/
				

				int p2 = 0;
				for (int p1 = 0; p1 < treatedList.size(); ++p1) {
					double q1 = treatedList.get(p1).treatedPropensity;
					double q2;
					/*
					double mindiff = 1.0;
					int minp2 = -1;
					
					for (int j =  0; j < untreatedList.size(); ++j) {
						if (untreatedMark[j] == true) continue;
						if (Math.abs(untreatedList.get(j).treatedPropensity - q1) < mindiff) {
							mindiff = Math.abs(untreatedList.get(j).treatedPropensity - q1);
							minp2 = j;
						}
					}
					*/
					while (p2 < untreatedList.size() && (untreatedList.get(p2).treatedPropensity < q1 || untreatedMark[p2] == true)) ++p2;
					int tmpp2 = p2;
					if (tmpp2 >= untreatedList.size()) tmpp2--;
					while (tmpp2 >= 0 && (untreatedList.get(tmpp2).treatedPropensity >= q1 || untreatedMark[tmpp2] ==  true)) --tmpp2;
					int fp2;
					if (p2 >= untreatedList.size()) {
						if (tmpp2 < 0) continue;
						else {
							fp2 = tmpp2;
						}
					}
					else {
						if (tmpp2 < 0) {
							fp2 = p2;
						}
						else {
							if (Math.abs(untreatedList.get(tmpp2).treatedPropensity - q1) < Math.abs(untreatedList.get(p2).treatedPropensity - q1)) {
								fp2 = tmpp2;
							}
							else {
								fp2 = p2;
							}
						}
					}
					q2 = untreatedList.get(fp2).treatedPropensity;
					//if (p2 < untreatedList.size()) q2 = untreatedList.get(p2).treatedPropensity;
					//else break;
					double diff = Math.abs(q2 - q1);
					/*
					if (p2 != 0 && untreatedMark[p2 - 1] == false) {
						if (diff > Math.abs(q1 - untreatedList.get(p2-1).treatedPropensity)) {
							--p2;
							diff = Math.abs(1 - untreatedList.get(p2).treatedPropensity);
						}
					}*/
					/*
					if (Math.abs(untreatedList.get(minp2).treatedPropensity - q1) != Math.abs(untreatedList.get(fp2).treatedPropensity - q1)) {
						System.out.println("ERROR p2! minp2=" + minp2 + ", p2=" + fp2);
						System.out.println(untreatedList.get(minp2).treatedPropensity + ", " + untreatedList.get(fp2).treatedPropensity);
						System.out.println(Math.abs(untreatedList.get(minp2).treatedPropensity - q1) + ", " + Math.abs(untreatedList.get(fp2).treatedPropensity - q1));
					}
					*/
					
					if (diff < propensitySigma * 2) {
						treatedMark[p1] = true;
						untreatedMark[fp2] = true;
						SlashdotUserNewsPair pair1 = untreatedList.get(fp2);
						SlashdotUserNewsPair pair2 = treatedList.get(p1);
						if (pair1.adopted == 0) ++untreatedCnt[0];
						else ++untreatedCnt[1];
						if (pair2.adopted == 0) ++treatedCnt[0];
						else ++treatedCnt[1];


						/*============DEBUG=========================
						System.out.println(pair1.treatedPropensity + "  vs  " + pair2.treatedPropensity);
						ArrayList<NetworkEdge> elist = nwdata.edgeFromList.get(pair1.uid);
						long tt1 = 0, tt2 = 0;
						int adfri = 0;
						for (NetworkEdge e: elist) {
							if (e.weight < 0) continue;
							if (e.v1id != pair1.uid) {
								System.out.println("ERROR1");
								System.exit(0);
							}
							if (cmdata.userComment.get(e.v1id).containsKey(pair1.nid)) tt1 = cmdata.userComment.get(e.v1id).get(pair1.nid).get(0).t.getTime();
							else tt1 = MAXTIME;
							if (cmdata.userComment.get(e.v2id).containsKey(pair1.nid)) tt2 = cmdata.userComment.get(e.v2id).get(pair1.nid).get(0).t.getTime();
							else tt2 = MAXTIME;
							//System.out.println(tt1 + "," + tt2);
							if (tt2 < tt1) {
								++adfri;
							}
						}
						if (cmdata.userComment.get(pair1.uid).containsKey(pair1.nid)) tt1 = cmdata.userComment.get(pair1.uid).get(pair1.nid).get(0).t.getTime();
						else tt1 = MAXTIME;
						System.out.println("Pair1, adopted friend: " + adfri);
						if (pair1.adopted == 1 && tt1 == MAXTIME) {
							System.out.println("ERROR2");
							System.exit(0);
						}
						else if (pair1.adopted == 0 && tt1 < MAXTIME) {
							System.out.println(tt1);
							System.out.println("ERROR3");
							System.exit(0);
						}
						
						elist = nwdata.edgeFromList.get(pair2.uid);
						adfri = 0;
						for (NetworkEdge e: elist) {
							if (e.weight < 0) continue;
							if (e.v1id != pair2.uid) System.out.println("ERROR1");
							if (cmdata.userComment.get(e.v1id).containsKey(pair2.nid)) tt1 = cmdata.userComment.get(e.v1id).get(pair2.nid).get(0).t.getTime();
							else tt1 = MAXTIME;
							if (cmdata.userComment.get(e.v2id).containsKey(pair2.nid)) tt2 = cmdata.userComment.get(e.v2id).get(pair2.nid).get(0).t.getTime();
							else tt2 = MAXTIME;
							//System.out.println(tt1 + "," + tt2);
							if (tt2 < tt1) {
								++adfri;
							}
						}
						if (cmdata.userComment.get(pair2.uid).containsKey(pair2.nid)) tt1 = cmdata.userComment.get(pair2.uid).get(pair2.nid).get(0).t.getTime();
						else tt1 = MAXTIME;
						System.out.println("Pair2, adopted friend: " + adfri);
						if (pair2.adopted == 1 && tt1 == MAXTIME) {
							System.out.println("ERROR2");
							System.exit(0);
						}
						else if (pair2.adopted == 0 && tt1 < MAXTIME) {
							System.out.println("ERROR3");
							System.exit(0);
						}
						
						System.out.println(pair1.ins);
						System.out.println(pair2.ins);
						System.out.println("===");
						=========================================*/
						
						//++p2;
					}
					//while (p2 != 0 && untreatedList.get(p2).treatedPropensity >= q1) --p2;
					p2 = tmpp2;
				}		
				long t5 = System.currentTimeMillis();
				bw.write(treatedCnt[1] + "\t" + treatedCnt[0] + "\t");
				bw.write(untreatedCnt[1] + "\t" + untreatedCnt[0] + "\t");
				bw.write("\n");
				bw.flush();
				long t6 = System.currentTimeMillis();
				timecounter[0] += t1 - t0;
				timecounter[1] += t2 - t1;
				timecounter[2] += t3 - t2;
				timecounter[3] += t4 - t3;
				timecounter[4] += t5 - t4;
				timecounter[5] += t6 - t5;
				long totaltime = 0;
				for (int i = 0; i < 6; ++i) {
					System.out.println("Timecost " + i + ":  " + timecounter[i] + " ms.");
					totaltime += timecounter[i];
				}
				System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
				System.out.println("Estimated time left: " + totaltime  * (nnum - nid - 1) / (nid + 1) / 1000 / 60 + " mins");
			}
			bw.close();

		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void StudyCorrByAssortativeMixingRandomSamplePos(String outputFileName, int threshold) {
		
		//Extract attributes and add to dataset.
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			System.out.println("Adding...");
			int nnum = cmdata.getNewsNum();
			long[] timecounter = new long[6];
			for (int nid = 0; nid < nnum; ++nid) {
				long t0 = System.currentTimeMillis();
				ArrayList<Attribute> attrv = new ArrayList<Attribute>();
				attrv.add(new Attribute("friendsnum"));
				attrv.add(new Attribute("foesnum"));
				attrv.add(new Attribute("fansnum"));
				attrv.add(new Attribute("freaksnum"));
				//attrv.add(new Attribute("adoptedfri"));
				ArrayList<String> classList = new ArrayList<String>();
				classList.add("YES");
				classList.add("NO");
				attrv.add(new Attribute("class", classList));
				Instances dataset = new Instances("trainset", attrv, nwdata.getNodeNum());
				dataset.setClassIndex(dataset.numAttributes() - 1);
				long t1 = System.currentTimeMillis();
				UpdateCheckAdoptedVector(nid);
				long t2 = System.currentTimeMillis();
				ArrayList<SlashdotUserNewsPair> userNewsPairList = new ArrayList<SlashdotUserNewsPair>();
				int posnum = 0, negnum = 0;
				for (int uid = 0; uid < nwdata.getNodeNum(); ++uid) {
					SlashdotUserNewsPair userNewsPair = new SlashdotUserNewsPair();
					userNewsPair.uid = uid;
					userNewsPair.nid = nid;
					long adptime = checkAdoptedVector[uid];
					if (adptime == MAXTIME) {
						userNewsPair.adopted = 0;
					}
					else {
						userNewsPair.adopted = 1;
					}
					
					SparseInstance ins = makePositiveInstance(uid, nid, dataset, threshold);
					if (ins == null) continue;
					//dataset.add(ins);
					userNewsPair.ins = ins;
					userNewsPair.treated = (int) (1 - ins.classValue());
					if (userNewsPair.treated == 1) {
						++posnum;
					}
					else {
						++negnum;
					}
					userNewsPairList.add(userNewsPair);
				}
				
				long t3 = System.currentTimeMillis();
				
				SlashdotUserNewsPair[] treatedArray = new SlashdotUserNewsPair[posnum];
				SlashdotUserNewsPair[] untreatedArray = new SlashdotUserNewsPair[negnum];
				posnum = 0;
				negnum = 0;
				for (int j = 0; j < userNewsPairList.size(); ++j) {
					SlashdotUserNewsPair pair = userNewsPairList.get(j);
					if (pair.treated == 1) treatedArray[posnum++] = pair;
					else untreatedArray[negnum++] = pair;
				}
				
				for (int j = 0; j < untreatedArray.length - 1; ++j) {
					int k = (int)(Math.floor(Math.random() * (untreatedArray.length - j - 1))) + j + 1;
					SlashdotUserNewsPair tmpPair = untreatedArray[k];
					untreatedArray[k] = untreatedArray[j];
					untreatedArray[j] = tmpPair;
				}

				long t4 = System.currentTimeMillis();
					
				int[] treatedCnt = new int[2];
				int[] untreatedCnt = new int[2];

				for (int j = 0; j < posnum; ++j) {
					SlashdotUserNewsPair pair1 = treatedArray[j];
					SlashdotUserNewsPair pair2 = untreatedArray[j];
					if (pair1.adopted == 1)  ++treatedCnt[1];
					else ++treatedCnt[0];
					if (pair2.adopted == 1)  ++untreatedCnt[1];
					else ++untreatedCnt[0];
				}				
				
				long t5 = System.currentTimeMillis();
				bw.write(treatedCnt[1] + "\t" + treatedCnt[0] + "\t");
				bw.write(untreatedCnt[1] + "\t" + untreatedCnt[0] + "\t");
				bw.write("\n");
				bw.flush();
				long t6 = System.currentTimeMillis();
				timecounter[0] += t1 - t0;
				timecounter[1] += t2 - t1;
				timecounter[2] += t3 - t2;
				timecounter[3] += t4 - t3;
				timecounter[4] += t5 - t4;
				timecounter[5] += t6 - t5;
				long totaltime = 0;
				for (int i = 0; i < 6; ++i) {
					System.out.println("Timecost " + i + ":  " + timecounter[i] + " ms.");
					totaltime += timecounter[i];
				}
				System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
				System.out.println("Estimated time left: " + totaltime  * (nnum - nid - 1) / (nid + 1) / 1000 / 60 + " mins");
			}
			bw.close();

		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	
	public void StudyPositiveInfluence(String outputFilename, int threshold) {
		int nnum = cmdata.getNewsNum();
		int unum = nwdata.getNodeNum();
		
		int[] frinum = new int[unum], foenum = new int[unum], fannum = new int[unum], frknum = new int[unum];
		int[] cmnnum = new int[unum], fricmnnum = new int[unum], foecmnnum = new int[unum], fancmnnum = new int[unum], frkcmnnum = new int[unum];
		for (int uid = 0; uid < unum; ++uid) {
			frinum[uid] = 0;
			foenum[uid] = 0;
			fannum[uid] = 0;
			frknum[uid] = 0;
			fricmnnum[uid] = 0;
			foecmnnum[uid] = 0;
			fancmnnum[uid] = 0;
			frkcmnnum[uid] = 0;
			cmnnum[uid] = cmdata.userComment.get(uid).size();
			ArrayList<NetworkEdge> elist;
			elist = nwdata.edgeFromList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++frinum[uid];
					fricmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
				else {
					++foenum[uid];
					foecmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
			}
			elist = nwdata.edgeToList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++fannum[uid];
					fancmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
				else {
					++frknum[uid];
					frkcmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
			}
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFilename));
			long ts0, ts1;
			long tstot = 0;
			for (int nid = 0; nid < nnum; ++nid) {
				ts0 = System.currentTimeMillis();
				MatchedSampleTester tester = new MatchedSampleTester();
				UpdateCheckAdoptedVector(nid);
				System.out.println("Adding...");
				for (int uid = 0; uid < unum; ++uid) {
					ArrayList<Double> featureVector = new ArrayList<Double>();
					boolean adopted, treated;
					featureVector.add((double)frinum[uid]);
					featureVector.add((double)foenum[uid]);
					featureVector.add((double)fannum[uid]);
					featureVector.add((double)frknum[uid]);	
					featureVector.add(frinum[uid]==0?0:(double)fricmnnum[uid]/frinum[uid]);
					featureVector.add(foenum[uid]==0?0:(double)foecmnnum[uid]/foenum[uid]);
					featureVector.add(fannum[uid]==0?0:(double)fancmnnum[uid]/fannum[uid]);
					featureVector.add(frknum[uid]==0?0:(double)frkcmnnum[uid]/frknum[uid]);	
					featureVector.add((double)cmnnum[uid]);
					ArrayList<NetworkEdge> elist;
					elist = nwdata.edgeFromList.get(uid);
					int adfri = 0, adfoe = 0;
					for (NetworkEdge e: elist) {
						if (e.weight < 0) {
							if (checkAdoptedVector[e.v2id] < MAXTIME) ++adfoe;
						}
						else {
							if (checkAdoptedVector[e.v2id] < MAXTIME) ++adfri;
						}
					}
					featureVector.add((double)adfoe);
					if (adfri >= threshold) {
						treated = true;
					}
					else if (adfri == 0){
						treated = false;
					}
					else continue;
					if (checkAdoptedVector[uid] == MAXTIME) adopted = false;
					else adopted = true;
					tester.AddSubject(featureVector, adopted, treated);
				}
				ArrayList<Integer> treatedRes = new ArrayList<Integer>();
				ArrayList<Integer> untreatedRes = new ArrayList<Integer>();
				System.out.println("Testing...");
				tester.MatchedSampleTest(treatedRes, untreatedRes);
				if (treatedRes.size() == 0) {
					bw.write("0\t0\t0\t0\t");
				}
				else {
					bw.write(treatedRes.get(1) + "\t" + treatedRes.get(0) + "\t");
					bw.write(untreatedRes.get(1) + "\t" + untreatedRes.get(0) + "\t");
				}
				bw.write("\n");
				bw.flush();
				ts1 = System.currentTimeMillis();
				System.out.println("News Timecost: " + (ts1-ts0) + "ms.");
				tstot += (ts1-ts0);
				System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
				System.out.println("Estimated time left: " + tstot  * (nnum - nid - 1) / (nid + 1) / 1000 / 60 + " mins");
			}
			bw.close();
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void StudyNegativeInfluence(String outputFilename, int threshold) {
		int nnum = cmdata.getNewsNum();
		int unum = nwdata.getNodeNum();
		
		int[] frinum = new int[unum], foenum = new int[unum], fannum = new int[unum], frknum = new int[unum];
		int[] cmnnum = new int[unum], fricmnnum = new int[unum], foecmnnum = new int[unum], fancmnnum = new int[unum], frkcmnnum = new int[unum];
		for (int uid = 0; uid < unum; ++uid) {
			frinum[uid] = 0;
			foenum[uid] = 0;
			fannum[uid] = 0;
			frknum[uid] = 0;
			fricmnnum[uid] = 0;
			foecmnnum[uid] = 0;
			fancmnnum[uid] = 0;
			frkcmnnum[uid] = 0;
			cmnnum[uid] = cmdata.userComment.get(uid).size();
			ArrayList<NetworkEdge> elist;
			elist = nwdata.edgeFromList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++frinum[uid];
					fricmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
				else {
					++foenum[uid];
					foecmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
			}
			elist = nwdata.edgeToList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++fannum[uid];
					fancmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
				else {
					++frknum[uid];
					frkcmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
			}
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFilename));
			long ts0, ts1;
			long tstot = 0;
			for (int nid = 0; nid < nnum; ++nid) {
				ts0 = System.currentTimeMillis();
				MatchedSampleTester tester = new MatchedSampleTester();
				UpdateCheckAdoptedVector(nid);
				System.out.println("Adding...");
				for (int uid = 0; uid < unum; ++uid) {
					ArrayList<Double> featureVector = new ArrayList<Double>();
					boolean adopted, treated;
					featureVector.add((double)frinum[uid]);
					featureVector.add((double)foenum[uid]);
					featureVector.add((double)fannum[uid]);
					featureVector.add((double)frknum[uid]);	
					featureVector.add(frinum[uid]==0?0:(double)fricmnnum[uid]/frinum[uid]);
					featureVector.add(foenum[uid]==0?0:(double)foecmnnum[uid]/foenum[uid]);
					featureVector.add(fannum[uid]==0?0:(double)fancmnnum[uid]/fannum[uid]);
					featureVector.add(frknum[uid]==0?0:(double)frkcmnnum[uid]/frknum[uid]);	
					featureVector.add((double)cmnnum[uid]);
					ArrayList<NetworkEdge> elist;
					elist = nwdata.edgeFromList.get(uid);
					int adfri = 0, adfoe = 0;
					for (NetworkEdge e: elist) {
						if (e.weight < 0) {
							if (checkAdoptedVector[e.v2id] < MAXTIME) ++adfoe;
						}
						else {
							if (checkAdoptedVector[e.v2id] < MAXTIME) ++adfri;
						}
					}
					featureVector.add((double)adfri);
					if (adfoe >= threshold) {
						treated = true;
					}
					else if (adfoe == 0){
						treated = false;
					}
					else continue;
					if (checkAdoptedVector[uid] == MAXTIME) adopted = false;
					else adopted = true;
					tester.AddSubject(featureVector, adopted, treated);
				}
				ArrayList<Integer> treatedRes = new ArrayList<Integer>();
				ArrayList<Integer> untreatedRes = new ArrayList<Integer>();
				System.out.println("Testing...");
				tester.MatchedSampleTest(treatedRes, untreatedRes);
				if (treatedRes.size() == 0) {
					bw.write("0\t0\t0\t0\t");
				}
				else {
					bw.write(treatedRes.get(1) + "\t" + treatedRes.get(0) + "\t");
					bw.write(untreatedRes.get(1) + "\t" + untreatedRes.get(0) + "\t");
				}
				bw.write("\n");
				bw.flush();
				ts1 = System.currentTimeMillis();
				System.out.println("News Timecost: " + (ts1-ts0) + "ms.");
				tstot += (ts1-ts0);
				System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
				System.out.println("Estimated time left: " + tstot  * (nnum - nid - 1) / (nid + 1) / 1000 / 60 + " mins");
			}
			bw.close();
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void StudyNegativeInfluenceWithSlices(String outputFilename, int threshold, long timeSliceLength, int timeSliceNum) {
		int nnum = cmdata.getNewsNum();
		int unum = nwdata.getNodeNum();
		
		int[] frinum = new int[unum], foenum = new int[unum], fannum = new int[unum], frknum = new int[unum];
		int[] cmnnum = new int[unum], fricmnnum = new int[unum], foecmnnum = new int[unum], fancmnnum = new int[unum], frkcmnnum = new int[unum];
		for (int uid = 0; uid < unum; ++uid) {
			frinum[uid] = 0;
			foenum[uid] = 0;
			fannum[uid] = 0;
			frknum[uid] = 0;
			fricmnnum[uid] = 0;
			foecmnnum[uid] = 0;
			fancmnnum[uid] = 0;
			frkcmnnum[uid] = 0;
			cmnnum[uid] = cmdata.userComment.get(uid).size();
			ArrayList<NetworkEdge> elist;
			elist = nwdata.edgeFromList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++frinum[uid];
					fricmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
				else {
					++foenum[uid];
					foecmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
			}
			elist = nwdata.edgeToList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++fannum[uid];
					fancmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
				else {
					++frknum[uid];
					frkcmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
			}
		}
		try {
			BufferedWriter[] bw = new BufferedWriter[3];
			for (int slice = 0; slice < timeSliceNum; ++slice) {
				bw[slice] = new BufferedWriter(new FileWriter(outputFilename + "_slices_" + slice + ".txt"));
			}
			long ts0, ts1;
			long tstot = 0;
			for (int nid = 0; nid < nnum; ++nid) {
				UpdateCheckAdoptedVector(nid);
				long earliestTime = MAXTIME;
				int commentUserNum = 0;
				for (int uid = 0; uid < unum; ++uid) {
					if (checkAdoptedVector[uid] < earliestTime) earliestTime = checkAdoptedVector[uid];
					if (checkAdoptedVector[uid] < MAXTIME) {
						++commentUserNum;
						//Date date = new Date(checkAdoptedVector[uid]); 
						//System.out.println(date.toGMTString());
					}
				}
				
				long lowerTime = earliestTime;
				long upperTime = earliestTime + timeSliceLength;
				for (int slice = 0; slice < timeSliceNum; ++slice) {
					ts0 = System.currentTimeMillis();
					MatchedSampleTester tester = new MatchedSampleTester();
					for (int uid = 0; uid < unum; ++uid) {
						if (checkAdoptedVector[uid] < lowerTime) continue;
						ArrayList<Double> featureVector = new ArrayList<Double>();
						boolean adopted, treated;
						featureVector.add((double)frinum[uid]);
						featureVector.add((double)foenum[uid]);
						featureVector.add((double)fannum[uid]);
						featureVector.add((double)frknum[uid]);	
						featureVector.add(frinum[uid]==0?0:(double)fricmnnum[uid]/frinum[uid]);
						featureVector.add(foenum[uid]==0?0:(double)foecmnnum[uid]/foenum[uid]);
						featureVector.add(fannum[uid]==0?0:(double)fancmnnum[uid]/fannum[uid]);
						featureVector.add(frknum[uid]==0?0:(double)frkcmnnum[uid]/frknum[uid]);	
						featureVector.add((double)cmnnum[uid]);
						ArrayList<NetworkEdge> elist;
						elist = nwdata.edgeFromList.get(uid);
						int adfri = 0, adfoe = 0;
						for (NetworkEdge e: elist) {
							if (e.weight < 0) {
								if (checkAdoptedVector[e.v2id] < upperTime && checkAdoptedVector[e.v2id] < checkAdoptedVector[uid]) ++adfoe;
							}
							else {
								if (checkAdoptedVector[e.v2id] < upperTime && checkAdoptedVector[e.v2id] < checkAdoptedVector[uid]) ++adfri;
							}
						}
						featureVector.add((double)adfri);
						if (adfoe >= threshold) {
							treated = true;
						}
						else if (adfoe == 0){
							treated = false;
						}
						else continue;
						if (checkAdoptedVector[uid] >= upperTime) adopted = false;
						else adopted = true;
						tester.AddSubject(featureVector, adopted, treated);
					}
					ArrayList<Integer> treatedRes = new ArrayList<Integer>();
					ArrayList<Integer> untreatedRes = new ArrayList<Integer>();
					System.out.println("Testing...");
					tester.MatchedSampleTest(treatedRes, untreatedRes);
					if (treatedRes.size() == 0) {
						bw[slice].write("0\t0\t0\t0\t");
					}
					else {
						bw[slice].write(treatedRes.get(1) + "\t" + treatedRes.get(0) + "\t");
						bw[slice].write(untreatedRes.get(1) + "\t" + untreatedRes.get(0) + "\t");
					}
					bw[slice].write("\n");
					bw[slice].flush();
					upperTime += timeSliceLength;
					lowerTime += timeSliceLength;
					ts1 = System.currentTimeMillis();
					System.out.println("News Timecost: " + (ts1-ts0) + "ms.");
					tstot += (ts1-ts0);
					System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
					System.out.println("Estimated time left: " + tstot  * (nnum - nid - 1) / (nid + 1) / 1000 / 60 + " mins");
				}
			}
			for (int slice = 0; slice < timeSliceNum; ++slice) {
				bw[slice].close();
			}
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void StudyPositiveInfluenceNewsBatch(String outputFilename, int threshold) {
		int nnum = cmdata.getNewsNum();
		int unum = nwdata.getNodeNum();
		
		int[] frinum = new int[unum], foenum = new int[unum], fannum = new int[unum], frknum = new int[unum];
		int[] cmnnum = new int[unum], fricmnnum = new int[unum], foecmnnum = new int[unum], fancmnnum = new int[unum], frkcmnnum = new int[unum];
		
		System.out.println("Initialization.");
		for (int uid = 0; uid < unum; ++uid) {
			frinum[uid] = 0;
			foenum[uid] = 0;
			fannum[uid] = 0;
			frknum[uid] = 0;
			fricmnnum[uid] = 0;
			foecmnnum[uid] = 0;
			fancmnnum[uid] = 0;
			frkcmnnum[uid] = 0;
			cmnnum[uid] = cmdata.userComment.get(uid).size();
			ArrayList<NetworkEdge> elist;
			elist = nwdata.edgeFromList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++frinum[uid];
					fricmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
				else {
					++foenum[uid];
					foecmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
			}
			elist = nwdata.edgeToList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++fannum[uid];
					fancmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
				else {
					++frknum[uid];
					frkcmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
			}
		}
		
		
		System.out.println("Testing...");
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFilename));
			long ts0, ts1;
			long tstot = 0;
			
			ts0 = System.currentTimeMillis();
			MatchedSampleTester tester = new MatchedSampleTester();
			final int BATCH = 1;
			
			for (int nid = 0; nid < nnum; ++nid) {
				if (nid % BATCH == 0) tester = new MatchedSampleTester();
				UpdateCheckAdoptedVector(nid);
				System.out.println("Adding...");
				for (int uid = 0; uid < unum; ++uid) {
					ArrayList<Double> featureVector = new ArrayList<Double>();
					boolean adopted, treated;
					featureVector.add((double)frinum[uid]);
					featureVector.add((double)foenum[uid]);
					featureVector.add((double)fannum[uid]);
					featureVector.add((double)frknum[uid]);	
					featureVector.add((double)fricmnnum[uid]/frinum[uid]);
					featureVector.add((double)foecmnnum[uid]/foenum[uid]);
					featureVector.add((double)fancmnnum[uid]/fannum[uid]);
					featureVector.add((double)frkcmnnum[uid]/frknum[uid]);	
					featureVector.add((double)cmnnum[uid]);
					ArrayList<NetworkEdge> elist;
					elist = nwdata.edgeFromList.get(uid);
					int adfri = 0, adfoe = 0;
					for (NetworkEdge e: elist) {
						if (e.weight < 0) {
							if (checkAdoptedVector[e.v2id] < checkAdoptedVector[uid]) ++adfoe;
						}
						else {
							if (checkAdoptedVector[e.v2id] < checkAdoptedVector[uid]) ++adfri;
						}
					}
					featureVector.add((double)adfoe);
					if (adfri == threshold) {
						treated = true;
					}
					else if (adfri == 0){
						treated = false;
					}
					else continue;
					if (checkAdoptedVector[uid] == MAXTIME) adopted = false;
					else adopted = true;
					tester.AddSubject(featureVector, adopted, treated);
				}
				ts1 = System.currentTimeMillis();
				System.out.println("News Timecost: " + (ts1-ts0) + "ms.");
				tstot += (ts1-ts0);
				System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
				System.out.println("Estimated time left: " + tstot  * (nnum - nid - 1) / (nid + 1) / 1000 / 60 + " mins");
				if (nid == nnum - 1 || (nid + 1) % BATCH == 0) {
					System.out.println("Testing...");
					ArrayList<Integer> treatedRes = new ArrayList<Integer>();
					ArrayList<Integer> untreatedRes = new ArrayList<Integer>();
					tester.MatchedSampleTest(treatedRes, untreatedRes);
					if (treatedRes.size() == 0) {
						bw.write("0\t0\t0\t0\t");
					}
					else {
						bw.write(treatedRes.get(1) + "\t" + treatedRes.get(0) + "\t");
						bw.write(untreatedRes.get(1) + "\t" + untreatedRes.get(0) + "\t");
					}
					bw.write("\n");
					bw.flush();
				}
			}

			bw.close();
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	public void StudyPositiveInfluenceWithSlices(String outputFilename, int threshold, long timeSliceLength, int timeSliceNum) {
		int nnum = cmdata.getNewsNum();
		int unum = nwdata.getNodeNum();
		
		int[] frinum = new int[unum], foenum = new int[unum], fannum = new int[unum], frknum = new int[unum];
		int[] cmnnum = new int[unum], fricmnnum = new int[unum], foecmnnum = new int[unum], fancmnnum = new int[unum], frkcmnnum = new int[unum];
		for (int uid = 0; uid < unum; ++uid) {
			frinum[uid] = 0;
			foenum[uid] = 0;
			fannum[uid] = 0;
			frknum[uid] = 0;
			fricmnnum[uid] = 0;
			foecmnnum[uid] = 0;
			fancmnnum[uid] = 0;
			frkcmnnum[uid] = 0;
			cmnnum[uid] = cmdata.userComment.get(uid).size();
			ArrayList<NetworkEdge> elist;
			elist = nwdata.edgeFromList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++frinum[uid];
					fricmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
				else {
					++foenum[uid];
					foecmnnum[uid] += cmdata.userComment.get(e.v2id).size();
				}
			}
			elist = nwdata.edgeToList.get(uid);
			for (NetworkEdge e: elist) {
				if (e.weight > 0) {
					++fannum[uid];
					fancmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
				else {
					++frknum[uid];
					frkcmnnum[uid] += cmdata.userComment.get(e.v1id).size();
				}
			}
		}
		try {
			BufferedWriter[] bw = new BufferedWriter[3];
			for (int slice = 0; slice < timeSliceNum; ++slice) {
				bw[slice] = new BufferedWriter(new FileWriter(outputFilename + "_slices_" + slice + ".txt"));
			}
			long ts0, ts1;
			long tstot = 0;
			for (int nid = 0; nid < nnum; ++nid) {
				UpdateCheckAdoptedVector(nid);
				long earliestTime = MAXTIME;
				int commentUserNum = 0;
				for (int uid = 0; uid < unum; ++uid) {
					if (checkAdoptedVector[uid] < earliestTime) earliestTime = checkAdoptedVector[uid];
					if (checkAdoptedVector[uid] < MAXTIME) {
						++commentUserNum;
						//Date date = new Date(checkAdoptedVector[uid]); 
						//System.out.println(date.toGMTString());
					}
				}
				
				long lowerTime = earliestTime;
				long upperTime = earliestTime + timeSliceLength;
				for (int slice = 0; slice < timeSliceNum; ++slice) {
					ts0 = System.currentTimeMillis();
					MatchedSampleTester tester = new MatchedSampleTester();
					for (int uid = 0; uid < unum; ++uid) {
						if (checkAdoptedVector[uid] < lowerTime) continue;
						ArrayList<Double> featureVector = new ArrayList<Double>();
						boolean adopted, treated;
						featureVector.add((double)frinum[uid]);
						featureVector.add((double)foenum[uid]);
						featureVector.add((double)fannum[uid]);
						featureVector.add((double)frknum[uid]);	
						featureVector.add(frinum[uid]==0?0:(double)fricmnnum[uid]/frinum[uid]);
						featureVector.add(foenum[uid]==0?0:(double)foecmnnum[uid]/foenum[uid]);
						featureVector.add(fannum[uid]==0?0:(double)fancmnnum[uid]/fannum[uid]);
						featureVector.add(frknum[uid]==0?0:(double)frkcmnnum[uid]/frknum[uid]);	
						featureVector.add((double)cmnnum[uid]);
						ArrayList<NetworkEdge> elist;
						elist = nwdata.edgeFromList.get(uid);
						int adfri = 0, adfoe = 0;
						for (NetworkEdge e: elist) {
							if (e.weight < 0) {
								if (checkAdoptedVector[e.v2id] < upperTime && checkAdoptedVector[e.v2id] < checkAdoptedVector[uid]) ++adfoe;
							}
							else {
								if (checkAdoptedVector[e.v2id] < upperTime && checkAdoptedVector[e.v2id] < checkAdoptedVector[uid]) ++adfri;
							}
						}
						featureVector.add((double)adfoe);
						if (adfri >= threshold) {
							treated = true;
						}
						else if (adfri == 0){
							treated = false;
						}
						else continue;
						if (checkAdoptedVector[uid] >= upperTime) adopted = false;
						else adopted = true;
						tester.AddSubject(featureVector, adopted, treated);
					}
					ArrayList<Integer> treatedRes = new ArrayList<Integer>();
					ArrayList<Integer> untreatedRes = new ArrayList<Integer>();
					System.out.println("Testing...");
					tester.MatchedSampleTest(treatedRes, untreatedRes);
					if (treatedRes.size() == 0) {
						bw[slice].write("0\t0\t0\t0\t");
					}
					else {
						bw[slice].write(treatedRes.get(1) + "\t" + treatedRes.get(0) + "\t");
						bw[slice].write(untreatedRes.get(1) + "\t" + untreatedRes.get(0) + "\t");
					}
					bw[slice].write("\n");
					bw[slice].flush();
					upperTime += timeSliceLength;
					lowerTime += timeSliceLength;
					ts1 = System.currentTimeMillis();
					System.out.println("News Timecost: " + (ts1-ts0) + "ms.");
					tstot += (ts1-ts0);
					System.out.println("Threshold=" + threshold + "   Checked news: " + (nid+1) + "/" + nnum + ":" + (double)(nid+1)/nnum*100 + "percents...");
					System.out.println("Estimated time left: " + tstot  * (nnum - nid - 1) / (nid + 1) / 1000 / 60 + " mins");
				}
			}
			for (int slice = 0; slice < timeSliceNum; ++slice) {
				bw[slice].close();
			}
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	
	static public void main(String[] args) {
		for (int sid = 0; sid < args.length; ++sid) {
			System.out.println(args[sid]);
		}
		//String nwdataFileName = "D:\\Users\\honglei\\negative\\data\\slashdot_network\\total.txt";
		//String cmdataFileName = "D:\\Users\\honglei\\negative\\data\\slashdot_content\\usercomment.txt";
		//String outFileName = "result2.txt";
		String nwdataFileName = "D:\\SRTWORK\\negative\\data\\slashdot_network.txt";
		String cmdataFileName = "D:\\SRTWORK\\negative\\data\\slashdot_comment_sampled.txt";
		
		SlashdotStudy slashdotStudy = new SlashdotStudy(nwdataFileName, cmdataFileName);
		//slashdotStudy.StudyNegativeCorrByAssortativeMixing(outFileName);
		//slashdotStudy.StudyCorrByAssortativeMixing(outFileName);
		//for (int threshold = 1; threshold <= 5; ++threshold) {
		String task = args[0];
		int threshold = Integer.parseInt(args[1]);
		if (task.equalsIgnoreCase("negative")) {
			String outFileName = "D:\\SRTWORK\\negative\\study\\influence_test\\matched_sample_neg_new_" + threshold + ".txt";
			//String outFileName = "matched_sample_neg_" + threshold + ".txt";
			//slashdotStudy.StudyNegativeInfluence(outFileName, threshold);
			//slashdotStudy.StudyNegativeInfluence(outFileName, threshold);
			slashdotStudy.StudyNegativeInfluenceWithSlices(outFileName, threshold, 1000*60*60*1, 3);
		}
		else if (task.equalsIgnoreCase("positive")){
			String outFileName = "D:\\SRTWORK\\negative\\study\\influence_test\\matched_sample_pos_new_" + threshold + ".txt";
			//String outFileName = "matched_sample_pos_" + threshold + ".txt";
			slashdotStudy.StudyPositiveInfluenceWithSlices(outFileName, threshold, 1000*60*60*1, 3);
		}
		//String outFileName = "matched_sample_pos_" + threshold + ".txt";
		//slashdotStudy.StudyCorrByAssortativeMixingMatchedSamplePos(outFileName, threshold);
		//slashdotStudy.StudyCorrByAssortativeMixingMatchedSamplePos(outFileName, threshold);
		
		//}
	}
	
}
