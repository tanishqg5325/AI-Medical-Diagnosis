#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <assert.h>
#define pb push_back
using namespace std;

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<string> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<float> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
	// Constructor- a node is initialised with its name and its categories
    Graph_Node(string name, int n, vector<string> vals)
	{
		Node_Name = name;
		nvalues = n;
		values = vals;
	}

	string get_name()
	{
		return Node_Name;
	}

	vector<int> get_children()
	{
		return Children;
	}

	vector<string> get_Parents()
	{
		return Parents;
	}

	vector<float> get_CPT()
	{
		return CPT;
	}

	int get_nvalues()
	{
		return nvalues;
	}

	vector<string> get_values()
	{
		return values;
	}

	void set_CPT(vector<float> new_CPT)
	{
		CPT.clear();
		CPT = new_CPT;
	}

    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents = Parent_Nodes;
    }

    // add another node in a graph as a child of this node
    bool add_child(int new_child_index)
    {
        for(int i=0;i<Children.size();i++)
            if(Children[i] == new_child_index)
                return 0;
        Children.pb(new_child_index);
        return 1;
    }
};


 // The whole network represted as a list of nodes
class network{

	list<Graph_Node> Pres_Graph;

public:

	void addNode(Graph_Node node)
	{
		Pres_Graph.pb(node);
	}
    
	int netSize()
	{
		return Pres_Graph.size();
	}

    // get the index of node with a given name
    int get_index(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name) == 0)
                return count;
            count++;
        }
        return -1;
    }
	
	// get the node at nth index
    list<Graph_Node>::iterator get_nth_node(int n)
    {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(count == n) return listIt;
            count++;
        }
        return listIt; 
    }
	
    //get the iterator of a node with a given name
    list<Graph_Node>::iterator search_node(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
            if(listIt->get_name().compare(val_name) == 0)
                return listIt;
        cout<<"node not found\n";
        return listIt;
    }
};

network read_network(string file_name)
{
	network Alarm;
	string line, temp, name;
  	ifstream myfile(file_name); 
  	vector<string> values;
    while(!myfile.eof())
    {
    	stringstream ss;
    	getline(myfile, line);
    	ss.str(line);
    	ss >> temp;
    	if(temp.compare("variable") == 0)
    	{
    		ss >> name;
    		getline (myfile,line);
    		stringstream ss2;
    		ss2.str(line);
    		for(int i=0;i<4;i++)
    			ss2 >> temp;
    		values.clear();
    		while(temp.compare("};") != 0)
    		{
    			values.pb(temp);
    			ss2 >> temp;
    		}
    		Graph_Node new_node(name, values.size(), values);
    		Alarm.addNode(new_node);
    	}
    	else if(temp.compare("probability") == 0)
    	{
    		ss >> temp >> temp;
            list<Graph_Node>::iterator listIt, listIt1;
    		listIt = Alarm.search_node(temp);
            int index = Alarm.get_index(temp);
            ss >> temp;
            values.clear();
    		while(temp.compare(")") != 0)
    		{
                listIt1 = Alarm.search_node(temp);
                listIt1->add_child(index);
    			values.pb(temp);
    			ss >> temp;
    		}
            listIt->set_Parents(values);
    		getline (myfile,line);
    		stringstream ss2;
            
    		ss2.str(line);
    		ss2 >> temp >> temp;
    		vector<float> curr_CPT;
            string::size_type sz;
    		while(temp.compare(";") != 0)
    		{
    			curr_CPT.pb(atof(temp.c_str()));
    			ss2 >> temp;
    		}
            listIt->set_CPT(curr_CPT);
    	}
    }
    myfile.close();
  	return Alarm;
}

vector<string> processLine(string line, int n)
{
    vector<string> ans; ans.reserve(n);
    stringstream ss(line); string word;
    while(ss >> word) ans.pb(word);
    return ans;
}

int main(int argc, char const *argv[])
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    assert(argc == 3);

    string bif_file_name = argv[1], data_file_name = argv[2]; 
    network Alarm = read_network(bif_file_name);
    // list<Graph_Node>::iterator it = Alarm.get_nth_node(0);
    // for(auto i : it->get_values()) cout << i << "\n";

    ifstream data_file(data_file_name);
    int n = Alarm.netSize();
    vector<vector<string>> data;
    string line;
    while(getline(data_file, line))
    {
        if(line == "") break;
        vector<string> row = processLine(line, n);
        data.pb(row);
    }
    data_file.close();
    // for(string str : data[0]) cout << str << " "; cout << "\n";
    // for(auto &v : data) assert(v.size() == n);
    
    
    return 0;
}
