#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>

class Node;
class Network;

class Link
{
	public:
		double weight;
		int dest;
};

class Node 
{
	public:
		double selection;
		vector<Link> Links;
		int clique;
		int id;
		
		double occupancy, new_occupancy;
				
		void generateLinks(int nlinks, double bias, Network N);
		void normalizeLinks();		
};

class Network
{
	public:
		vector<Node> Nodes;
		void generateTimeSequence(char *outName, int length);
		void saveTransitionMatrix(char *outName);
		void initialConditions();
};

void Node::generateLinks(int nlinks, double bias, Network N)
{
	int i,j;
	
	for (i=0;i<nlinks;i++)
	{
		j=rand()%N.Nodes.size();
		
		Link L;
		
		L.weight = 0.5 + 0.5*bias*( 2*(N.Nodes[j].clique == clique) - 1 );
		L.dest = j;
		
		Links.push_back(L);
	}
}

void Node::normalizeLinks()
{
	double total = 0;
	int i;
	
	for (i=0;i<Links.size();i++)
	{
		total += Links[i].weight;
	}
	
	for (i=0;i<Links.size();i++)
		Links[i].weight /= total;
}

Network GenerateNetwork(int nnodes, int ncliques, double clique_bias)
{
	Network N;
	int i;
	
	for (i=0;i<nnodes;i++)
	{
		Node thisNode;
		
		thisNode.clique = rand()%ncliques;
		thisNode.selection = 1.1;
		thisNode.id = i;
		
		N.Nodes.push_back(thisNode);
	}
	
	int links_per_node = pow(nnodes, 0.5);
	
	for (i=0;i<nnodes;i++)
	{
		N.Nodes[i].generateLinks(links_per_node, clique_bias, N);
		N.Nodes[i].normalizeLinks();
	}	
	
	return N;
}

void Network::initialConditions()
{
	int i;
	
	for (i=0;i<Nodes.size();i++)
	{
		Nodes[i].occupancy = (rand()%1000001)/1000000.0;
	}
}

void Network::saveTransitionMatrix(char *FName)
{
	FILE *f=fopen(FName,"wb");
	int i,j,k;
	double val;
	
	for (j=0;j<Nodes.size();j++)
	{
		for (i=0;i<Nodes.size();i++)
		{
			val = 0;
			
			for (k=0;k<Nodes[i].Links.size();k++)
			{
				if (Nodes[Nodes[i].Links[k].dest].id == j)
					val += Nodes[i].Links[k].weight * Nodes[i].selection;
			}
			
			fprintf(f,"%.6g ", val);
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
}

void Network::generateTimeSequence(char *FName, int length)
{
	int t,i,j;
	FILE *f=fopen(FName,"wb");
	
	fprintf(f,"# %d\n",Nodes.size());
	fflush(f);
	
	initialConditions();
	
	for (t=0;t<length;t++)
	{
		for (i=0;i<Nodes.size();i++)
			Nodes[i].new_occupancy = 0;
			
		for (i=0;i<Nodes.size();i++)
			for (j=0;j<Nodes[i].Links.size();j++)
			{
				Nodes[Nodes[i].Links[j].dest].new_occupancy += Nodes[i].occupancy * Nodes[i].selection * Nodes[i].Links[j].weight;
			}

		for (i=0;i<Nodes.size();i++)
			Nodes[i].occupancy = Nodes[i].new_occupancy;
			
		fprintf(f,"%d ",t);
		for (i=0;i<Nodes.size();i++)
			fprintf(f,"%.6g ",Nodes[i].occupancy);
		fprintf(f,"\n");
	}
	
	fclose(f);
}

int main(int argc, char **argv)
{
	int i;
	char Str[512];
	Network N = GenerateNetwork(1000,8,0.99);
	
	N.saveTransitionMatrix("transition.mat");
	
	for (i=0;i<1000;i++)
	{
		sprintf(Str,"timeseries/%.4d.txt", i);
		N.generateTimeSequence(Str, 5000);
        cout << i << endl;
	}
}
