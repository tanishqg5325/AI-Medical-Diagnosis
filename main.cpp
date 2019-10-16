#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <assert.h>
#include <random>
#define pb push_back
using namespace std;

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<int> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<double> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

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

	vector<int>& get_children()
	{
		return Children;
	}

	vector<int>& get_Parents()
	{
		return Parents;
	}

	vector<double>& get_CPT()
	{
		return CPT;
	}

	int get_nvalues()
	{
		return nvalues;
	}

	vector<string>& get_values()
	{
		return values;
	}

	void set_CPT(vector<double> new_CPT)
	{
		CPT.clear();
		CPT = new_CPT;
	}

    void set_Parents(vector<int> Parent_Nodes)
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

	int getIndexForValue(string val)
	{
		for(int i=0;i<nvalues;i++)
			if(values[i] == val)
				return i;
		return -1;
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
            values.clear(); vector<int> parents;
    		while(temp.compare(")") != 0)
    		{
                listIt1 = Alarm.search_node(temp);
                listIt1->add_child(index);
				parents.pb(Alarm.get_index(temp));
    			// values.pb(temp);
    			ss >> temp;
    		}
            listIt->set_Parents(parents);
    		getline (myfile,line);
    		stringstream ss2;
            
    		ss2.str(line);
    		ss2 >> temp >> temp;
    		vector<double> curr_CPT;
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

vector<list<Graph_Node>::iterator> nodes;

void fillAllNodes(network &Alarm)
{
	list<Graph_Node>::iterator it = Alarm.get_nth_node(0);
	int n = Alarm.netSize(); nodes.reserve(n);
	for(int i=0;i<n;i++) nodes.pb(it++);
}

default_random_engine generator;
uniform_real_distribution<double> distribution(0.0, 1.0);

// given a probability distribution in the form of a vector, returns the index of sampled value
int sample(vector<double> &v)
{
	int n = v.size();
	if(n == 1) return 0;
	double number = distribution(generator), sum = 0;
	for(int i=0;i<n;i++)
	{
		sum += v[i];
		if(number <= sum) return i;
	}
	assert(1 == 0);
}

double computeProbabilityGivenParents(network &alarm, vector<int> &row, int nodeIndex)
{
	list<Graph_Node>::iterator node = nodes[nodeIndex];
	vector<int> &parents = node->get_Parents();
	int numParents = parents.size(), cptIndex = 0, prev = 1;
	for(int l=numParents-1;l>=0;l--)
	{
		cptIndex += row[parents[l]] * prev;
		prev *= nodes[parents[l]]->get_nvalues();
	}
	cptIndex += row[nodeIndex] * prev;
	assert(cptIndex < node->get_CPT().size());
	return node->get_CPT()[cptIndex];
}

void computeProbabilityGivenAllOther(network &alarm, vector<int> &row, int missingIndex)
{
	int n = alarm.netSize();
	list<Graph_Node>::iterator missingNode = nodes[missingIndex];
	int k = missingNode->get_nvalues(); double normalizing_constant;
	vector<double> missingNodeDistribution(k);
	vector<int> &children = missingNode->get_children();

	for(int i=0;i<k;i++)
	{
		row[missingIndex] = i;
		missingNodeDistribution[i] = computeProbabilityGivenParents(alarm, row, missingIndex);
		for(int &j : children)
			missingNodeDistribution[i] *= computeProbabilityGivenParents(alarm, row, j);
		normalizing_constant += missingNodeDistribution[i];
	}

	for(int i=0;i<k;i++)
		missingNodeDistribution[i] /= normalizing_constant;
	row[missingIndex] = sample(missingNodeDistribution);
}

void fillMissingValues(network &alarm, vector<vector<int>> &data, vector<int> &missingIndexes)
{
	int rows = data.size();
	for(int i=0;i<rows;i++)
		computeProbabilityGivenAllOther(alarm, data[i], missingIndexes[i]);
}

int main(int argc, char const *argv[])
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    assert(argc == 3);

    string bif_file_name = argv[1], data_file_name = argv[2]; 
    network Alarm = read_network(bif_file_name);
	fillAllNodes(Alarm); int n = Alarm.netSize();

    // list<Graph_Node>::iterator it = nodes[1];
    // for(auto i : it->get_values()) cout << i << "\n";

    ifstream data_file(data_file_name);
    vector<vector<int>> data;
	vector<int> missingIndexes;
    string line, word; vector<int> row(n);
    while(getline(data_file, line))
    {
		stringstream ss(line); int i = 0;
		while(ss >> word)
		{
			int j = nodes[i]->getIndexForValue(word);
			if(j == -1) missingIndexes.pb(i);
			row[i++] = j;
		}
		// assert(i == n);
        data.pb(row);
    }
    data_file.close();
    // for(int i : data[0]) cout << i << " "; cout << "\n";
    // for(auto &v : data) assert(v.size() == n);
    // cout << data.size() << "\n";
    return 0;
}
