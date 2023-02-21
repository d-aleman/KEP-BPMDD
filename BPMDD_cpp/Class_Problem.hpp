//
//  Class_Problem.hpp
//  CycleFormulation
//
//  Created by Carolina Riascos Alvarez on 2019-12-02.
//  Copyright © 2019 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef Class_Problem_hpp
#define Class_Problem_hpp

#include <time.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <map>
#include <deque>
#include <ilcplex/ilocplex.h>


//#include <lemon/list_graph.h> //Clase para crear la lista, hay varias pero esta es rápida
//#include <lemon/dijkstra.h> //Clase para que usemos dijkstra
//#include <lemon/connectivity.h>
ILOSTLBEGIN


struct state{
    int id;     //fake
    int actual; // arco real
    map<state, vector<int>> mActual;   // donde voy dependiendo de donde vengo map<state, int> actual ARCO!;
    double MaxWeight;   //pues no sé jeje es para el longest
    vector<int>::iterator Itestoy;  //vecino donde estoy
    vector<int>::iterator ItEnd;  //vecino donde estoy
    vector<int>::iterator ItmActualCualestoy; //
    vector<int>::iterator ItmActualCualEnd;
    map<int,bool> vVisited;   // arcos visitados
    vector<int> vSuce;     // la posición de los estados suce
    vector<int> vPrede;    // la posición de los estados prede
    map<int,bool>NoMergeables;
    bool firstTime;
    bool pathLlego;
    int predLongest; //position of state giving the longest path in a previous layer
    int mActualCual; // para ver qué mactual es del -8, maximum actual(node) given the longest path for -8, el -8 no lo tiene, sólo los estados que se conectan con -8
    vector<int> VpredLongest;
    vector<int> VmActualCual;
    bool cadeRepe;
    bool meConectoDos;
    int totalVeci;
    int vamosVeci;
    vector<tuple<int,int,double>> vTupleVarios;
    //////Functions//////
    state (int _id, int _actu, map<int,vector<int>>& Comp, state& vengo);
    state (int _id, int _origen, map<int,vector<int>>& Comp);
    void limpiar();
    void sizeBytes();
    void eliminar88(vector<state>& vLosNew);
    bool operator<(const state & Esta) const{ return id < Esta.id;}
};

struct cycle{
    vector<int> acycle;
    float fraction;
    double TotalWeight;
    int indexAllC;
    cycle(){}
    cycle(vector<int> _cycle, float _fraction, int _indexAllC){acycle = _cycle, fraction = _fraction, indexAllC = _indexAllC;}
    cycle(double _weight, vector<int> _cycle, int _indexAllC){acycle = _cycle, TotalWeight = _weight, fraction = 1, indexAllC = _indexAllC;}
};

struct arc{
    //pair<int,int> tailhead;
    vector<int> WhereinMDDs; // in which copies of the MDD the arc is repeated
    //bool BranchOnIt; //if the tail is repeated
    float fraction;
    
    //Functions//
    arc(){fraction = 0;}
    arc(float _fraction, vector<int> _WhereinMDDs){fraction = _fraction, WhereinMDDs = _WhereinMDDs;}
};
        
struct node{
    int id;
    float weight;
    //bool BranchOnIt; //if the tail is repeated
    vector<int> WhereinMDDs; // in which copies of the MDD the arc is repeated
    
    //Functions//
    node(int _id, float _weight, vector<int> _WhereinMDDs){id = _id, weight = _weight, WhereinMDDs = _WhereinMDDs;}
};

struct CopySelection{
    bool selected = false;
    double max_score = 0;
    int sel_freq = 0;
    float sel_prob = 0;
};

struct CycleSelection{
    bool selected = false;
    double max_score = 0;
    float sel_prob = 0.5;
};
        
struct BBnode{
    int id;
    int level;
    int nextNeighbor;
    double UB; //Given by the Lagrangian relaxation
    double LB; //Given by a feasible solution of the Lagrangian Relaxation
    pair<int,int>parentNode; // position in Tree of the parent node, i->row j->column
    vector<pair<int,int>>ArcGroupZero;
    vector<pair<int,int>>ArcGroupOne;
    map<int,vector<int>>NodeGroupZero; // copy MDD and node removed from that copy
    vector<int>NodeGroupOne;
    vector<int>NeighNextLevel; //Position of neighbors in the next level of Tree
    map<int,cycle>CyclesInSolution;
    vector<vector<int>>ExtraCycles;
    vector<double>ExtraWeight;
    
    //Constructor//
    BBnode(int _id, int _level, int _nextNeighbor, double _UB, vector<pair<int,int>>_ArcGroupZero, vector<pair<int,int>>_ArcGroupOne, pair<int,int>_parentNode);
    BBnode(int _id, int _level, int _nextNeighbor, double _UB, map<int,vector<int>>_NodeGroupZero, vector<int>_NodeGroupOne, pair<int,int>_parentNode);
};
        
class Problem{
public:
    IloEnv env;
    IloModel MP;
    IloInt Pairs;
    IloInt NDDs;//altruists
    IloInt Nodes; //altruists + pairs
    IloInt NumArcs;
    IloInt CycleLength;//K
    IloInt ChainLength;//F
    IloInt FeedNode;
    IloInt biggestDegree = 0;
    IloInt TotalCycles = 0;
    IloNumArray2 WeightMatrix;
    IloNumArray2 AdjacencyList;//Successors
    string FileName;
    string FolderName;
    fstream file;
    map<int,vector<int>>CycleNode;
    vector<vector<int>>PredList;
    map<pair<int,int>,double>Weights;
    vector<map<int,vector<int>>>Components; //This array is removed in LongestPath
    map<int,vector<int>>Comp;//Successors
    map<int,vector<int>>CompPred;//predecessor
    //vector<vector<vector<Network>>>Net;//índice del grafo, layer, nodes in layer
    map<int,bool>NonPresentNodes; //node,first layer encountered. This array is removed in Lagrangian
    vector<vector<map<int,int>>>Posid;//
    vector<map<int,vector<vector<int>>>>ListOfEquals;
    
    
    ///////////////////Lagrangian Relaxation////////////////
    int exploredBBnodes = 0;
    int CadaCuantoMerge = 0;
    double LongestTime = 0;
    double MPTime = 0;
    double CFTime = 0;
    double SPTime = 0;
    double MDDTime = 0;
    double TBBtime = 0;
    double UpperBound;
    long MergedNodes = 0;
    long nstatesMDD = 0;
    double graphDensity;
    double CFObj = 0;
    double diffheu = 0;
    bool FirstTimeChains = true;
    bool FirstTimeCycles = true;
    bool FirstConversion = true;
    vector<int>nNodes;
    vector<int>nArcs;
    vector<vector<int>>vFeasSol;
    vector<vector<int>>vBestSol;
    map<int,bool>NeverConstrained;
    IloArray<IloNumVarArray> x;
    IloArray<IloNumVarArray> y;
    IloModel ch;
    IloModel cy;
    IloObjective ObjCH;
    IloObjective ObjCY;
    string DegreeType;
    IloNumVarArray pi;
    IloNumArray solpi;
    IloNumVarArray theta;
    IloNumArray sol_theta;
    IloBool FirstTimeLPath = true;
    IloBool FirstTimeMIPCY = true;
    IloBool UpdateMaxW = true;
    IloBool FastExit = false;
    IloNum TimeLimit;
    IloCplex cex;
    IloNumVarArray IncumbentVars;
    IloNumArray IncumbentVals;
    IloRangeArray cons;
    IloNumArray2 x_sol;
    int countSame = 0;
    vector<pair<int,double>>Order;
    vector<vector<pair<int,double>>> LongestPath;
    vector<map<int,vector<int>>>JustVisitedByState; // level<actual,vector<pos>>
    vector<pair<vector<state>, vector<vector<int>>>> vMDD;   //vector de copias, first vector<state> and second vector<vector<int>> vLvEstados
    vector<pair<int,int>> ordenFor;
    tuple<int,int,double,int,double> tCicloFin;
    vector<tuple<int,int,double,int,double>> ttCicloFin;
    void CreateIncumbent(IloNumVarArray& z, vector<int> thisSol);
    vector<vector<int>>LeftCycles();
    
    ////////Functions/////
    vector<pair<int,int>> ordenDegre();
    void SetName(IloNumVar& var, const char* prefix, IloInt i);
    pair<vector<state>, vector<vector<int>>> BuildCycleMDD(int iori, bool& hay, map<int,bool>& mOrden, int voyPeque, int cuantosMDDsPeque);   //dani
    pair<vector<state>, vector<vector<int>>> BuildChainMDD(int iori, bool& hay, map<int,bool>& mOrden);
    pair<int, double> predMax(vector<state>& vLosNew, int indexestoy, bool cadena, int CuantosPahtsQueres);
    int DijktrasAlgorithm(map<int,vector<int>>& Comp2, int source, int sink);
    //void mergeNodes (vector<state>& vLosNew, int lv, int estoy);
    //void mergeNodes();
    void BuildMDDs();
    void llenarComp();
    void Lagrange();
    void CountStatesinMDDs();
    void HeadingLagRelaxation();
    void PrintLag(string status);
    void InterchangeableNodes(map<int,vector<int>>Comp, int origen);
    double Heuristic();
    
    
    ////Functions////
    Problem (string _FolderName, string _FileName, string _DegreeType, IloInt _cycleLength, IloInt _chainLength, IloNum _TimeLimit, IloNumArray2 WeightMatrix, IloNumArray2 AdjacencyList, IloInt Pairs, IloInt NDDs);
    void FindCycles();
    IloBool IsxinStack (IloInt test, vector<int>& xinTrial);
    void inVLcompSP (IloEnv env, IloNumArray2 AdjaList, IloInt origin);
    void PairWiseShortestPath(int sp);
    IloNum PathWeight (vector<int>& Stack);
    IloInt degree(IloInt Node, string DegreeType);
    vector<map<int,int>> FindEquals(map<int,vector<int>>& Comp2, int oriMDD);
    int ReadData();
    void MainCycleFinder();
    vector<pair<int,int>> OptimalVFS();
    vector<vector<int>> SubCycleFinder (IloEnv env, IloNumArray2 AdjaList, IloInt origin);
    
    
    bool IsxinStackTuple(IloInt node, vector<tuple<int,int,double>> v);
    void FindLongestPath();
    void NewFindLongestPath();
    tuple<int,int,double,int,double> maxMenosOcho(vector<state>& vLosNew, state& estoy, int i, bool cadena, int CuantosPahtsQueres);
    vector<tuple<int,int,double,int,double>> NewmaxMenosOcho(vector<state>& vLosNew, state& estoy, int i, bool cadena);
    vector<pair<int, double>> NewpredMax(vector<state>& vLosNew, int indexestoy, bool cadena);
    void SetName2(IloNumVar& var, const char* prefix, IloInt i, IloInt j);
    vector<vector<int>> Chains();
    vector<vector<int>> FindChains(IloNumArray2 x_sol);
    vector<vector<int>> FindFeasibleCycles(IloNumArray2& x_sol);
    vector<int> FindSubtour(IloNumArray2 x_sol);
    bool isInOrderm(vector<pair<int,double>>::iterator& itVeci, int veci);

    /////////////////////////////BRANCH&BOUND//////////////////////////
    vector<vector<BBnode>>Tree;
    vector<pair<pair<int,int>,arc>>LagArcs;
    vector<node>LagNodes;
    vector<pair<int,int>>ArcBranches;
    vector<int>NodeBranches;
    vector<pair<int,vector<cycle>>> ListCycles; //copy, cycles in copy
    float GlobalLB = -10000000;
    float GlobalUB = -10000000;
    clock_t StartSolTime;
    IloNum TrunTime = 0;
    bool FirstBBnode = true;
    bool insideBB = false;
    bool TimedOut = false;
    bool ThirdPhase = false;
    int BranchingMethod = 2; //CHANGE ACCORDINGLY: 1->By Node, 2->By Arc
    pair<int,int>VoyBB; //row, column in Tree
    int posMDD = -1;
    IloBool FirstIterLag = true;
    map<int,cycle>CyclesInSol;
    double AccumulatedWeight = 0;
    map<int, vector<int>>NodeinCycle;
    vector<float>NodeinCycleWeight;
    double nCycleBBnode = 0;
    IloNumArray z_sol;
    vector<vector<int>>AllCycles;
    vector<vector<int>>AllCyclesVFS;
    vector<double>AllCyclesWeight;
    IloInt LastinAllCycles;
    map<int,bool> mOrden;
    
    //Printing variables
    int nColCycles = 0;
    int nColChains = 0;
    int itPH2LP = 0;
    int itPH2MIP = 0;
    int itPH2Sep = 0;
    int nColChainPH2LP = 0;
    int nColChainPH2MIP = 0;
    int nColChainPH2Sep = 0;
    int nColCyclePH2 = 0;
    int nColCyclePH3Sch = 0;
    int nColCyclePH3MIP = 0;
    int itPH3Search = 0;
    int itPH3MIP = 0;
    double UB_PH1 = 0;
//    int nChainsSol = 0;
//    int nCyclesSol = 0;
    double TimePH2 = 0;
    double TimePH3 = 0;
    
    /////////////////////////////Functions//////////////////////////
    void BBTree();
    void getChildren();
    void CycleFormulation();
    void PrintMDDs();
        
    //void buscarLongestPath(vector<state>& vLosNew, vector<vector<int>>& vLvEstados, map<pair<int,int>,double>& Weights, bool cadena);
private:
};

bool sortNodes(node& node1, node& node2);
bool sortNodesRep(node& node1, node& node2);
bool sortMap (pair<int,vector<cycle>>&m1, pair<int,vector<cycle>>&m2);
bool sortCycles(cycle& c1, cycle& c2);
bool sortArcs(pair<pair<int,int>,arc>& a1, pair<pair<int,int>,arc>& a2);
bool sortArcsFrac(pair<pair<int,int>,arc>& a1, pair<pair<int,int>,arc>& a2);
bool sortBBnode(BBnode& n1, BBnode& n2);
bool sortLowHigh (int i,int j);
bool sortOrder(pair<int,double>& a1, pair<int,double>& a2);
bool sortsCopies(CopySelection& c1, CopySelection& c2);

#endif /* Class_Problem_hpp */
