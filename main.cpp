#include <iostream>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <iomanip>
#define pb push_back
using namespace std;

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

class Graph_Node;
vector<list<Graph_Node>::iterator> nodes;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node{

private:
	string Node_Name;  // Variable name 
	vector<int> Children; // Children of a particular node - these are index of nodes in graph.
	vector<int> Parents; // Parents of a particular node- note these are names of parents
	int nvalues;  // Number of categories a variable represented by this node can take
	vector<string> values; // Categories of possible values
	vector<double> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning
	vector<double> count;

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

	vector<double>& get_count()
	{
		return count;
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

network alarm;

void read_network(string file_name)
{
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
    		getline(myfile, line);
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
			alarm.addNode(new_node);
    	}
    	else if(temp.compare("probability") == 0)
    	{
    		ss >> temp >> temp;
            list<Graph_Node>::iterator listIt, listIt1;
			listIt = alarm.search_node(temp);
            int index = alarm.get_index(temp);
            ss >> temp;
            values.clear(); vector<int> parents;
    		while(temp.compare(")") != 0)
    		{
                listIt1 = alarm.search_node(temp);
                listIt1->add_child(index);
				parents.pb(alarm.get_index(temp));
    			ss >> temp;
    		}
            listIt->set_Parents(parents);
    		getline(myfile, line);
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
}

void write_output(string input_file_name) {
	ifstream input_file(input_file_name);
	string line, temp, output;
	list<Graph_Node>::iterator it;
	while(!input_file.eof()) {
		stringstream ss;
		getline(input_file, line);
		ss.str(line); ss >> temp;
		output += line + "\n";
		if(temp.compare("probability") == 0) {
			ss >> temp >> temp;
			getline(input_file, line); getline(input_file, line);
			output += "\ttable ";
			it = alarm.search_node(temp);
			for(double &i : it->get_CPT()) {
				if(fabs(i) < 1e-4) i = 0.0001;
				stringstream ss2;
				ss2 << fixed << setprecision(4) << i;
				output += ss2.str() + " ";
			}
			output += ";\n}\n";
		}
	}
	input_file.close(); output.pop_back();
	ofstream output_file("solved_alarm.bif");
	output_file << output;
	output_file.close();
}

void fillAllNodes()
{
	list<Graph_Node>::iterator it = alarm.get_nth_node(0);
	int n = alarm.netSize(); nodes.reserve(n);
	for(int i=0;i<n;i++) nodes.pb(it++);
}

int getCPTIndex(vector<int> &row, int nodeIndex)
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
	return cptIndex;
}

vector<double> computeProbabilityGivenAllOther(vector<int> &row, int missingIndex)
{
	int n = alarm.netSize();
	list<Graph_Node>::iterator missingNode = nodes[missingIndex];
	int k = missingNode->get_nvalues(); double normalizing_constant;
	vector<double> missingNodeDistribution(k);
	vector<int> &children = missingNode->get_children();

	for(int i=0;i<k;i++)
	{
		row[missingIndex] = i;
		missingNodeDistribution[i] = missingNode->get_CPT()[getCPTIndex(row, missingIndex)];
		for(int &j : children)
			missingNodeDistribution[i] *= nodes[j]->get_CPT()[getCPTIndex(row, j)];
		normalizing_constant += missingNodeDistribution[i];
	}

	for(int i=0;i<k;i++)
		missingNodeDistribution[i] /= normalizing_constant;
	return missingNodeDistribution;
}

vector<vector<double>> fillMissingValues(vector<vector<int>> &data, vector<int> &missingIndexes)
{
	int rows = data.size();
	vector<vector<double>> miss_dist;
	miss_dist.reserve(rows);
	vector<double> tmp;
	for(int i=0;i<rows;i++)
	{
		if(missingIndexes[i] != -1) {
			miss_dist.pb(computeProbabilityGivenAllOther(data[i], missingIndexes[i]));
		}
		else {
			miss_dist.pb(tmp);
		}
	}
	return miss_dist;
}

void normalize()
{
	int n = alarm.netSize(), prod_parent_values;
	double sum;
	for(int i=0;i<n;i++)
	{
		prod_parent_values = nodes[i]->get_CPT().size() / nodes[i]->get_nvalues();
		for(int j=0;j<prod_parent_values;j++)
		{
			sum = 0;
			for(int k=j;k<nodes[i]->get_CPT().size();k+=prod_parent_values)
				sum += nodes[i]->get_count()[k];
			if(sum == 0) {
				nodes[i]->get_count()[j] = sum = 1;
			}
			for(int k=j;k<nodes[i]->get_CPT().size();k+=prod_parent_values)
				nodes[i]->get_CPT()[k] =  1.0 * nodes[i]->get_count()[k] / sum;
		}
	}
}

void update_CPT(vector<vector<int>> &data, vector<vector<double>> &prev_miss_dist, vector<vector<double>> &now_miss_dist, vector<int> missingIndexes)
{
	int rows = data.size(), new_value;
	for(int i=0;i<rows;i++)
	{
		if(prev_miss_dist[i].empty() || prev_miss_dist[i] == now_miss_dist[i]) continue;
		int k = nodes[missingIndexes[i]]->get_nvalues();
		for(int j=0;j<k;j++)
		{
			data[i][missingIndexes[i]] = j;
			nodes[missingIndexes[i]]->get_count()[getCPTIndex(data[i], missingIndexes[i])] += (now_miss_dist[i][j] - prev_miss_dist[i][j]);
			for(int &c : nodes[missingIndexes[i]]->get_children()) {
				nodes[c]->get_count()[getCPTIndex(data[i], c)] += (now_miss_dist[i][j] - prev_miss_dist[i][j]);
			}
		}
	}
	normalize();
}

void initialize_CPT(vector<vector<int>> &data, vector<vector<double>> &prev_miss_dist, vector<int> &missingIndexes, double smoothing_parameter)
{
	int rows = data.size(); int n = alarm.netSize();
	for(int i=0;i<n;i++)
	{
		nodes[i]->get_count().resize(nodes[i]->get_CPT().size());
		fill(nodes[i]->get_count().begin(), nodes[i]->get_count().end(), smoothing_parameter);
	}
	for(int i=0;i<rows;i++)
	{
		if(missingIndexes[i] != -1)
		{
			int k = nodes[missingIndexes[i]]->get_nvalues();
			for(int j=0;j<k;j++)
			{
				data[i][missingIndexes[i]] = j;
				for(int l=0;l<n;l++)
					nodes[l]->get_count()[getCPTIndex(data[i], l)] += prev_miss_dist[i][j];
			}
		}
		else
		{
			for(int l=0;l<n;l++)
				nodes[l]->get_count()[getCPTIndex(data[i], l)]++;
		}
	}
	normalize();
}

vector<vector<double>> computePrior(vector<vector<int>> &data)
{
	int n = alarm.netSize(), rows = data.size();
	vector<vector<double>> prior(n);
	for(int i=0;i<n;i++) prior[i].resize(nodes[i]->get_nvalues());
	for(int i=0;i<rows;i++)
		for(int j=0;j<n;j++)
			if(data[i][j] != -1) prior[j][data[i][j]]++;
	for(int i=0;i<n;i++)
	{
		double sum = 0; int k = nodes[i]->get_nvalues();
		for(int j=0;j<k;j++) sum += prior[i][j];
		for(int j=0;j<k;j++) prior[i][j] /= sum;
	}
	return prior;
}

void handleWithNoParents(double smoothing_parameter)
{
	int prod_parent_values, n = alarm.netSize(); double sum;
	for(int i=0;i<n;i++)
	{
		prod_parent_values = nodes[i]->get_CPT().size() / nodes[i]->get_nvalues();
		vector<bool> v(nodes[i]->get_CPT().size(), 0);
		for(int j=0;j<prod_parent_values;j++)
		{
			bool flag = 1;
			for(int k=j;k<nodes[i]->get_CPT().size();k+=prod_parent_values)
				if(nodes[i]->get_count()[k] != smoothing_parameter)
					{flag = 0; break;}

			if(flag == 1) {
				for(int k=j;k<nodes[i]->get_CPT().size();k+=prod_parent_values)
					v[k] = 1;
			}
		}
		for(int j=0;j<nodes[i]->get_nvalues();j++)
		{
			double sum = 0; int cnt = 0;
			for(int k=0;k<prod_parent_values;k++)
			{
				int ind = j * prod_parent_values + k;
				if(!v[ind])
				{
					cnt++;
					sum += nodes[i]->get_CPT()[ind];
				}
			}
			if(cnt == 0) continue;
			for(int k=0;k<prod_parent_values;k++)
			{
				int ind = j * prod_parent_values + k;
				if(v[ind])
					nodes[i]->get_CPT()[ind] = sum / cnt;
			}
		}
	}
}

int main(int argc, char const *argv[])
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    assert(argc == 3);

    string bif_file_name = argv[1], data_file_name = argv[2]; 
    read_network(bif_file_name);
	fillAllNodes(); int n = alarm.netSize();

    ifstream data_file(data_file_name);
    vector<vector<int>> data;
	vector<vector<double>> prev_miss_dist, now_miss_dist;
	vector<int> missingIndexes;
	data.reserve(11100); missingIndexes.reserve(11100); prev_miss_dist.reserve(11100);
    string line, word; vector<int> row(n);
	bool isMissingFound; int i, j;
    while(getline(data_file, line))
    {
		stringstream ss(line); i = 0;
		isMissingFound = false;
		while(ss >> word)
		{
			j = nodes[i]->getIndexForValue(word);
			if(j == -1) {
				missingIndexes.pb(i);
				isMissingFound = true;
			}
			row[i++] = j;
		}
		if(!isMissingFound) missingIndexes.pb(-1);
        data.pb(row);
    }
    data_file.close(); vector<double> temp; int rows = data.size();
	vector<vector<double>> prior = computePrior(data);
	for(int i=0;i<rows;i++)
	{
		if(missingIndexes[i] != -1) prev_miss_dist.pb(prior[missingIndexes[i]]);
		else prev_miss_dist.pb(temp);
	}

	double smoothing_parameter = 0.05;
	initialize_CPT(data, prev_miss_dist, missingIndexes, smoothing_parameter);
	double max_time = 115.0;
	while(((double)clock()/CLOCKS_PER_SEC) < max_time)
	// for(int i=0;i<500;i++)
	{
		now_miss_dist = fillMissingValues(data, missingIndexes);
		update_CPT(data, prev_miss_dist, now_miss_dist, missingIndexes);
		prev_miss_dist = now_miss_dist;
	}

	handleWithNoParents(smoothing_parameter);

	write_output(bif_file_name);
    return 0;
}
