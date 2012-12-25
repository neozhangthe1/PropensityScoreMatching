package Slashdot;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;



public class SlashdotContent {
	
	private static HashMap<String, Integer> monthMap = new HashMap<String, Integer>();
	
	
	
	
	ArrayList<SlashdotComment> commentList = new ArrayList<SlashdotComment>();
	
	public SlashdotContent() {
		monthMap.put("January", 0);
		monthMap.put("February", 1);
		monthMap.put("March", 2);
		monthMap.put("April", 3);
		monthMap.put("May", 4);
		monthMap.put("June", 5);
		monthMap.put("July", 6);
		monthMap.put("August", 7);
		monthMap.put("September", 8);
		monthMap.put("October", 9);
		monthMap.put("November", 10);
		monthMap.put("December", 11);
	}
	
	
	@SuppressWarnings("deprecation")
	public void parseContent(String filename, String outfile) {
		int lastyear = 2012;
		try {
			BufferedReader bw = new BufferedReader(new FileReader(filename));
			BufferedWriter wr = new BufferedWriter(new FileWriter(outfile));
			String line;
			int linenum = 0;
			while ((line = bw.readLine()) != null){
				linenum++;
				if (linenum % 100000 == 0) {
					System.out.println("Processed: " + linenum + " lines...");
				}
				String[] strList = line.split("[|][|][|]");
				SlashdotComment slashdotComment = new SlashdotComment();
				slashdotComment.userId = Integer.parseInt(strList[0]);
				slashdotComment.newsId = Integer.parseInt(strList[1]);
				slashdotComment.cat1 = strList[2];
				slashdotComment.cat2 = strList[3];
				String[] dateStr = strList[4].split(" @");
				String[] dateTokens = dateStr[0].split(" |,");
				//System.out.println("uid:" + slashdotComment.userId);
				//System.out.println("newsid:" + slashdotComment.newsId);
				//System.out.println("cat1:" + slashdotComment.cat1);
				//System.out.println("cat2:" + slashdotComment.cat2);
				//System.out.println("datetokens:");
				Date date = new Date();
				
				/*if (dateTokens.length < 5) {
					for (int j = 0; j < dateTokens.length; ++j) {
						System.out.println(dateTokens[j]);
					}
				}*/
				date.setMonth(monthMap.get(dateTokens[2]));
				date.setDate(Integer.parseInt(dateTokens[3]));
				if (dateTokens.length < 5) {
					date.setYear(lastyear - 1900);
				}
				else {
					lastyear = Integer.parseInt(dateTokens[4]);
					date.setYear(Integer.parseInt(dateTokens[4]) - 1900);
				}
				

				//System.out.println("time:");
				String timeStr = dateStr[1];
				//System.out.println(timeStr);
				
				int h = Integer.parseInt(timeStr.substring(0, 2));
				int m = Integer.parseInt(timeStr.substring(3, 5));
				if (timeStr.substring(5, 7).equals("PM") && h != 12) {
					h += 12;
				}
				date.setHours(h);
				date.setMinutes(m);
				date.setSeconds(0);
				//System.out.println("parsed time:");
				//System.out.println(date.toString());
				slashdotComment.t = date;
				//System.out.println("=======");
				wr.write(slashdotComment.userId + "\t" + slashdotComment.newsId + "\t" + (slashdotComment.t.getTime()/1000) + "\n");
				wr.flush();
			}
			
			wr.close();
			bw.close();
			
		} catch(Exception e) {
			e.printStackTrace();
			return;
		}
		
	}
	
	public static void main(String[] args) {
		//String filename = "D:\\SRTWORK\\negative\\src\\tmp\\content_sample.txt";
		//String outfile = "D:\\SRTWORK\\negative\\src\\tmp\\content_sample_out.txt";
		String filename = "D:\\Users\\honglei\\influence\\exp\\slashdot\\slashdot\\content.txt";
		String outfile = "D:\\Users\\honglei\\influence\\exp\\slashdot\\slashdot\\usercomment.txt";
		SlashdotContent slashdotContent = new SlashdotContent();
		slashdotContent.parseContent(filename, outfile);
	}
}
