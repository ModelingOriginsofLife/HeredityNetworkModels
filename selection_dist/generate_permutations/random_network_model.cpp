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
		double selection, perturbation;
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
		
		void iterate();
		void generatePerturbation(double amp);
		void generateTimeSequence(char *outName, int length);
		void saveTransitionMatrix(char *outName);
		void initialConditions();
		void normalize();
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
		thisNode.selection = 1.1+0.03*(rand()%200001-100000.0)/100000.0;
		thisNode.id = i;
		
		N.Nodes.push_back(thisNode);
	}
	
	int links_per_node = 2000;
	
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
	
//	Nodes[rand()%Nodes.size()].occupancy = 1;
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
					val += Nodes[i].Links[k].weight * ( Nodes[i].selection + Nodes[i].perturbation );
			}
			
			fprintf(f,"%.6g ", val);
		}
		fprintf(f,"\n");
	}
	
	fclose(f);
}

void Network::generatePerturbation(double amp)
{
	int i;
	
	for (i=0;i<Nodes.size();i++)
	{
		if (rand()%2==0)
			Nodes[i].perturbation = -Nodes[i].selection + 0.015;
		else Nodes[i].perturbation = 0;
	}
}

void Network::normalize()
{
	double total=0;
	int i;
	
	for (i=0;i<Nodes.size();i++)
		total+=Nodes[i].occupancy;
	
	for (i=0;i<Nodes.size();i++)
		Nodes[i].occupancy/=total;
}

void Network::iterate()
{
	int i,j;
	
	for (i=0;i<Nodes.size();i++)
		Nodes[i].new_occupancy = 0;
		
	for (i=0;i<Nodes.size();i++)
		for (j=0;j<Nodes[i].Links.size();j++)
		{
			Nodes[Nodes[i].Links[j].dest].new_occupancy += Nodes[i].occupancy * (Nodes[i].selection + 0.01 * ( rand()%2000001-1000000.0)/1000000.0  + Nodes[i].perturbation ) * Nodes[i].Links[j].weight;
		}

	for (i=0;i<Nodes.size();i++)
		Nodes[i].occupancy = Nodes[i].new_occupancy;
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
		iterate();
			
		fprintf(f,"%d ",t);
		for (i=0;i<Nodes.size();i++)
			fprintf(f,"%.6g ",Nodes[i].occupancy);
		fprintf(f,"\n");
	}
	
	fclose(f);
}

void loadEigenvectors(vector<double> *first, vector<double> *second, vector<double> *firstI, vector<double> *secondI)
{
	first->clear(); second->clear();	
	firstI->clear(); secondI->clear();
		
	FILE *f=fopen("first.txt","rb");
	double fbuf;
	
	while (fscanf(f,"%lf\n",&fbuf)!=EOF)
	{
		first->push_back(fbuf);
	}
	
	fclose(f);

	f=fopen("second.txt","rb");
	
	while (fscanf(f,"%lf\n",&fbuf)!=EOF)
	{
		second->push_back(fbuf);
	}
	
	fclose(f);

	f=fopen("firstI.txt","rb");
	
	while (fscanf(f,"%lf\n",&fbuf)!=EOF)
	{
		firstI->push_back(fbuf);
	}
	
	fclose(f);

	f=fopen("secondI.txt","rb");
	
	while (fscanf(f,"%lf\n",&fbuf)!=EOF)
	{
		secondI->push_back(fbuf);
	}
	
	fclose(f);
}

int main(int argc, char **argv)
{
	int i,j;
	char Str[512];
	Network Base = GenerateNetwork(10000,300,0.999);
	Network Perturbed;
	
	srand(time(NULL));
	
	for (j=0;j<1000;j++)
	{
		Base.initialConditions();
		Base.generatePerturbation(1.5);
		
		for (i=0;i<250;i++)
		{
			Base.iterate();
			Base.normalize();
		}
		Base.normalize();

		FILE *f=fopen("vectors.txt","a");
		for (i=0;i<Base.Nodes.size();i++)
			fprintf(f,"%.6g ",Base.Nodes[i].occupancy);
		fprintf(f,"\n");
		fclose(f);
	}
	
/*	for (j=0;j<1000;j++)
	{
		Perturbed.Nodes = Base.Nodes;
		
		for (i=0;i<Perturbed.Nodes.size();i++)
			Perturbed.Nodes[i].occupancy += 0.5*(rand()%1000001)/1000000.0;
			
		Perturbed.normalize();
		
		for (i=0;i<200;i++)
			Perturbed.iterate();
		
		Perturbed.normalize();
		
		FILE *f=fopen("vectors.txt","a");
		for (i=0;i<Perturbed.Nodes.size();i++)
			fprintf(f,"%.6g ",Perturbed.Nodes[i].occupancy);
		fprintf(f,"\n");
		fclose(f);
	}*/
}
