package Network;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class NetworkData {
	public HashMap<Integer, Integer> name2id;
	public ArrayList<Integer> id2name;
	public ArrayList<NetworkEdge> edgeList;
	public ArrayList<ArrayList<NetworkEdge>> edgeFromList;
	public ArrayList<ArrayList<NetworkEdge>> edgeToList;
	
	public NetworkData() {
		edgeList = new ArrayList<NetworkEdge>();
		edgeFromList = new ArrayList<ArrayList<NetworkEdge>>();
		edgeToList = new ArrayList<ArrayList<NetworkEdge>>();
		name2id = new HashMap<Integer, Integer>();
		id2name = new ArrayList<Integer>();
	}
	
	
	public int getNodeNum() {
		return id2name.size();
	}
	
	public int getEdgeNum() {
		return edgeList.size();
	}
	
	public int getIdWithAdding(int name) {
		if (name2id.containsKey(name)){
			return name2id.get(name);
		}
		else {
			id2name.add(name);
			edgeFromList.add(new ArrayList<NetworkEdge>());
			edgeToList.add(new ArrayList<NetworkEdge>());
			name2id.put(name, id2name.size() - 1);
			
			return id2name.size() - 1;
		}
	}
	
	public void ReadNetwork(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			String str;
			while (null != (str = br.readLine())) {
				String[] tokens = str.split("\\t");
				boolean flag = false;
				for (int j = 0; j < tokens.length; ++j) if (tokens[j].equals("")) flag = true;
				if (flag == true) continue;
				if (tokens.length < 3) continue;
				int v1name = Integer.parseInt(tokens[0]);
				int v2name = Integer.parseInt(tokens[1]);
				int v1id = getIdWithAdding(v1name);
				int v2id = getIdWithAdding(v2name);
				int sign = Integer.parseInt(tokens[2]);
				NetworkEdge edge = new NetworkEdge();
				edge.v1id = v1id;
				edge.v2id = v2id;
				edge.weight = sign;
				edge.eid = edgeList.size();
				edgeList.add(edge);
				edgeFromList.get(v1id).add(edge);
				edgeToList.get(v2id).add(edge);
			}
			br.close();
		}catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	
	static public void main(String[] args) {
		NetworkData nwdata = new NetworkData();
		nwdata.ReadNetwork("D:\\SRTWORK\\negative\\data\\slashdot_network.txt");
	}
}
