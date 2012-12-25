package Slashdot;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Map.Entry;

import Network.NetworkData;


public class SlashdotCommentData {
	HashMap<Integer, Integer> newsname2id;
	ArrayList<Integer> newsid2name;
	ArrayList<Long> newsFirstCommentTime;
	NetworkData nwdata;
	
	ArrayList<SlashdotComment> commentsList;
	
	/**
	 * userComment[user_id][news_id] = a list of user_id's comments on news_id.
	 */
	ArrayList<HashMap<Integer, ArrayList<SlashdotComment>>> userComment;
	
	
	
	
	
	public SlashdotCommentData(NetworkData nwdata) {
		this.nwdata = nwdata;
		newsid2name = new ArrayList<Integer>();
		newsname2id = new HashMap<Integer, Integer>();
		commentsList = new ArrayList<SlashdotComment>();
		userComment = new ArrayList<HashMap<Integer,ArrayList<SlashdotComment>>>();
		for (int i = 0; i < nwdata.getNodeNum(); ++i) userComment.add(new HashMap<Integer, ArrayList<SlashdotComment>>());
	}
	
	
	public int getNewsNum() {
		return newsid2name.size();
	}
	
	
	public int getIdWithAdding(int name) {
		if (newsname2id.containsKey(name)){
			return newsname2id.get(name);
		}
		else {
			newsid2name.add(name);
			newsname2id.put(name, newsid2name.size() - 1);
			return newsid2name.size() - 1;
		}
	}
	
	
	public void readSlashdotComments(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String str;
			while (null != (str = br.readLine())) {
				String[] tokens = str.split("\\t");
				if (tokens.length < 3) continue;
				int uname = Integer.parseInt(tokens[0]);
				int nname = Integer.parseInt(tokens[1]);
				if (!nwdata.name2id.containsKey(uname)) continue;
				int uid = nwdata.name2id.get(uname);
				int nid = getIdWithAdding(nname);
				long timenum = Integer.parseInt(tokens[2]);
				SlashdotComment comment = new SlashdotComment();
				comment.t = new Date(timenum * 1000);
				comment.userId = uid;
				comment.newsId = nid;
				commentsList.add(comment);
				if (!userComment.get(uid).containsKey(nid)) {
					userComment.get(uid).put(nid, new ArrayList<SlashdotComment>());
				}
				userComment.get(uid).get(nid).add(comment);
			}
			br.close();
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
		
		//Sort users' comments according to time
		for (int uid = 0; uid < nwdata.getNodeNum(); ++uid) {
			for (Entry<Integer, ArrayList<SlashdotComment>> e : userComment.get(uid).entrySet()) {
				ArrayList<SlashdotComment> clist = e.getValue();
				Collections.sort(clist, new Comparator<SlashdotComment>() {
					public int compare(SlashdotComment c1, SlashdotComment c2) {
						return c1.t.compareTo(c2.t);
					} 
				});
			}
		}
	}

}
