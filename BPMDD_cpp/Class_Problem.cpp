//
//  Class_Problem.cpp
//
//  Created by Carolina Riascos Alvarez on 2019-12-02.
//  Copyright © 2019 Carolina Riascos Alvarez. All rights reserved.
//

#include "Class_Problem_VFS.hpp"

Problem::Problem(string _FilePath, string _OutputPath, string _DegreeType, IloInt _cycleLength, IloInt _chainLength, IloNum _TimeLimit, IloNumArray2 _WeightMatrix, IloNumArray2 _AdjacencyList, IloInt _Pairs, IloInt _NDDs, string _Preference){
    FilePath = _FilePath;
    OutputPath = _OutputPath;
    DegreeType = _DegreeType;
    CycleLength = _cycleLength;
    ChainLength = _chainLength;
    TimeLimit = _TimeLimit;
    WeightMatrix = _WeightMatrix;
    AdjacencyList = _AdjacencyList;
    Pairs = _Pairs;
    NDDs = _NDDs;
    Preference = _Preference;
    tCicloFin = make_tuple(-1, -1, -10000, -1, -10000);
}
void Problem::SetName(IloNumVar& var, const char* prefix, IloInt i){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i);
    const char* varName = name.c_str();
    var.setName(varName);
}
IloInt Problem::degree(IloInt Node, string DegreeType){
    IloInt size;
    if (DegreeType == "Indegree"){
        size = PredList[Node].size();
        if (size > biggestDegree){
            FeedNode = Node;
            biggestDegree = PredList[Node].size();
        }
    }
    else if (DegreeType == "Degree"){
        size = AdjacencyList[Node].getSize();
        if (size > biggestDegree){
            FeedNode = Node;
            biggestDegree = PredList[Node].size();
        }
    }
    else{//Total degree
        size = PredList[Node].size() + AdjacencyList[Node].getSize();
        if (size > biggestDegree){
            FeedNode = Node;
            biggestDegree = PredList[Node].size() + AdjacencyList[Node].getSize();
        }
    }

    return FeedNode;
}
bool Problem::IsxinStackTuple(IloInt id, vector<tuple<int,int,double>> v){
    for (int h = 0; h < v.size(); h++){
        if (get<0>(v[h]) == id){
            return true;
        }
    }
    return false;
}
void Problem::PrintLag(string status){
    file.open(OutputPath, fstream::app);
    file << "Instance: " << FileName << endl;
    file << "Max Cyle Length: " << CycleLength << endl;
    file << "Max Chain Length: " << ChainLength << endl;
    file << "Degree Type: " << DegreeType << endl;
    file << "Preference: " << Preference << endl;
    file << "Time Taken: " << MDDTime + TBBtime << endl;
    file << "Time Limit: " << TimeLimit << endl;
    file << "Status: " << status << endl;
    file << "Objective Value: " << GlobalLB << endl;
    
    
    double density = ((Pairs)*(Pairs - 1) + (Pairs*NDDs));
    density = NumArcs/density*100;

    int c_chain = 0;
    int c_cycle = 0;
    int maxcycle = 0;
    int mincycle = 1000;
    int maxchain = 0;
    int minchain = 1000;
    int nTransChains = 0;
    int nTransCycles = 0;
    double aveLenChain = 0;
    double aveLenCycle = 0;
    int TotalHPatients = 0;
    for (int i = 0; i < PairsType.getSize(); i++) TotalHPatients++;
    
    //Print out
    for (int m = 0; m  < 2; m++){
        for (int i = 0; i < vBestSol.size(); i++){
            if (m == 1){
                if (vBestSol[i][0] > Pairs - 1){
                    c_chain++;
                    aveLenChain+= int(vBestSol[i].size() - 1);
                    file << "Chain " << c_chain << ":" << '\t';
                    for (int j = 0; j < vBestSol[i].size(); j++){
                        file << vBestSol[i][j] << '\t';
                    }
                    file << endl;
                    if (vBestSol[i].size() < minchain) minchain = int(vBestSol[i].size() - 1);
                    if (vBestSol[i].size() > maxchain) maxchain = int(vBestSol[i].size() - 1);
                    if (PairsType.getSize() > 0) StatisticsHighlySPatiens(vBestSol[i]);
                }
            }
            else{
                if (vBestSol[i][0] < Pairs){
                    c_cycle++;
                    aveLenCycle+= int(vBestSol[i].size());
                    file << "Cycle " << c_cycle << ":" << '\t';
                    for (int j = 0; j < vBestSol[i].size(); j++){
                        file << vBestSol[i][j] << '\t';
                    }
                    file << endl;
                    if (vBestSol[i].size() < mincycle) mincycle = int(vBestSol[i].size());
                    if (vBestSol[i].size() > maxcycle) maxcycle = int(vBestSol[i].size());
                    if (PairsType.getSize() > 0) StatisticsHighlySPatiens(vBestSol[i]);
                }
            }
        }
    }
    nTransChains = aveLenChain;
    nTransCycles = aveLenCycle;
    aveLenChain = aveLenChain/c_chain;
    aveLenCycle = aveLenCycle/c_cycle;
    
    //Finalize statistics on Highly Sensitized patients
    if (TotalHPatients > 0){
        perHinSol = perHinSol/TotalHPatients;
        perHinChains = perHinChains/TotalHPatients;
        perHin4PlusChains = perHin4PlusChains/TotalHPatients;
        perHinCycles = perHinCycles/TotalHPatients;
        perH4PlusCycles = perH4PlusCycles/TotalHPatients;
    }
    
    file << endl;
    file << "Average Chain Length: " << aveLenChain << endl;
    file << "Average Cycle Length: " << aveLenCycle << endl;
    file << "Percentage of highly-sensitized patients in solution (%): " << perHinSol*100 << endl;
    
    file.close();
    
}
void Problem::StatisticsHighlySPatiens(vector<int>match){
    if (match[0] > Pairs - 1){
        for (int i = 1; i < match.size(); i++){
            if (PairsType[match[i]] == 1){
                perHinSol++;
                perHinChains++;
                if (match.size() >= 5) perHin4PlusChains++; //Chains of size 4 or longer
            }
        }
    }
    else{
        for (int i = 0; i < match.size(); i++){
            if (PairsType[match[i]] == 1){
                perHinSol++;
                perHinCycles++;
                if (match.size() >= 4) perH4PlusCycles++; //Cycles of size 4 or longer
            }
        }
    }
}
state::state (int _id, int _actu, map<int,vector<int>>& Comp, state& vengo){
    id = _id; MaxWeight = 0; firstTime = true; pathLlego = false; predLongest = 0;
    mActual[vengo].push_back(_actu); actual = -88; meConectoDos = false;
    ItEnd = Comp[mActual[vengo].back()].end();
    Itestoy = Comp[mActual[vengo].back()].begin();
    mActualCual = -11;
}
state::state (int _id, int _actu, map<int,vector<int>>& Comp){
    id = _id; MaxWeight = 0; firstTime = true;  pathLlego = false; predLongest = 0; meConectoDos = false;
    if (id != -88){
        ItEnd = Comp[_actu].end();
        Itestoy = Comp[_actu].begin();
        totalVeci = int(Comp[_actu].size());
        vamosVeci = 0;
    }
    actual = _actu;
    mActualCual = -11;
}
BBnode::BBnode(int _id, int _level, int _nextNeighbor, double _UB, vector<pair<int,int>>_ArcGroupZero, vector<pair<int,int>>_ArcGroupOne, pair<int,int>_parentNode){
    id = _id;
    level = _level;
    nextNeighbor = _nextNeighbor;
    UB = _UB; //Given by the Lagrangian relaxation
    LB = -1000; //Given by a feasible solution of the Lagrangian Relaxation
    ArcGroupZero = _ArcGroupZero; // node, copy where a node is removed
    ArcGroupOne = _ArcGroupOne;   // arc, copy where an arc is removed
    parentNode = _parentNode;
    NeighNextLevel = vector<int>();
}
BBnode::BBnode(int _id, int _level, int _nextNeighbor, double _UB, map<int,vector<int>> _NodeGroupZero, vector<int> _NodeGroupOne, pair<int,int>_parentNode){
    id = _id;
    level = _level;
    nextNeighbor = _nextNeighbor;
    UB = _UB; //Given by the Lagrangian relaxation
    LB = -1000; //Given by a feasible solution of the Lagrangian Relaxation
    NodeGroupZero = _NodeGroupZero;   // arc, copy where an arc is removed
    NodeGroupOne = _NodeGroupOne;
    parentNode = _parentNode;
    NeighNextLevel = vector<int>();
}
void state::limpiar(){
    //clock_t start = clock();
    id = -88; actual = -88; Itestoy = ItEnd;
    MaxWeight = -88; vVisited.clear(); vSuce.clear(); vPrede.clear(); mActual.clear();
    //Tlimpiar += clock() - start;
}
void state::sizeBytes(){
    int total = 0; int s = 0;
    s = sizeof(id); total += s; cout << "id =\t" << s << endl;
    s = sizeof(actual); total += s; cout << "actual =\t" << s << endl;
    s = sizeof(mActual); total += s; cout << "mActual =\t" << s << endl;
    s = sizeof(MaxWeight); total += s; cout << "MaxWeight =\t" << s << endl;
    s = sizeof(Itestoy); total += s; cout << "Itestoy =\t" << s << endl;
    s = sizeof(ItEnd); total += s; cout << "ItEnd =\t" << s << endl;
    s = sizeof(vVisited); total += s; cout << "vVisited =\t" << s << endl;
    s = sizeof(vSuce); total += s; cout << "vSuce =\t" << s << endl;
    s = sizeof(vPrede); total += s; cout << "vPrede =\t" << s << endl;
    s = sizeof(firstTime); total += s; cout << "firstTime =\t" << s << endl;
    s = sizeof(pathLlego); total += s; cout << "pathLlego =\t" << s << endl;
    cout << "suma: " << total << " bytes" << endl;
    int real = sizeof(*this);
    cout << "real: " << real << " bytes = " << real/double(1000000) << " megabytes" <<  endl;
}
void state::eliminar88(vector<state>& vLosNew){   // de los prede y suce
    if (id != -88){
        for (auto i = 0; i < vPrede.size(); i++){
            if (vLosNew[vPrede[i]].id == -88) {
                vPrede.erase(vPrede.begin() + i);
            }
        }
        for (auto i = 0; i < vSuce.size(); i++){
            if (vLosNew[vSuce[i]].id == -88) {
                vSuce.erase(vSuce.begin() + i);
            }
        }
    }
}
//
void Problem::CountStatesinMDDs(){
    for (int i = 0; i < vMDD.size(); i++){
        nstatesMDD += vMDD[i].first.size();
    }
}
bool sortNodes(node& node1, node& node2){
    return (node1.weight > node2.weight);
}
bool sortNodesRep(node& node1, node& node2){
    return (node1.WhereinMDDs.size() > node2.WhereinMDDs.size());
}
bool sortMap (pair<int,vector<cycle>>&m1, pair<int,vector<cycle>>&m2){
    return (m1.second.size() > m2.second.size());
}
bool sortCycles(cycle& c1, cycle& c2){
    return (c1.fraction > c2.fraction);
}
bool sortArcs(pair<pair<int,int>,arc>& a1, pair<pair<int,int>,arc>& a2){
    if (a1.second.fraction < 1){
        return (a1.second.WhereinMDDs.size() > a2.second.WhereinMDDs.size());
    }
    else{
        return false;
    }
}
bool sortArcsFrac(pair<pair<int,int>,arc>& a1, pair<pair<int,int>,arc>& a2){
    float v1 = abs(a1.second.fraction - 0.5)/0.5;
    float v2 = abs(a2.second.fraction - 0.5)/0.5;
    return (v1 < v2);
}
bool sortBBnode(BBnode& n1, BBnode& n2){
    if (n1.UB > n2.UB){
        return true;
    }
    else if (n1.UB == n2.UB){
        return (n1.LB > n2.LB);
    }
    else{
        return false;
    }
}
bool sortLowHigh (int i,int j){
     return (i<j); 
}
bool sortOrder(pair<int,double>& a1, pair<int,double>& a2){
    return a1.second < a2.second;
}
bool sortsCopies(CopySelection& c1, CopySelection& c2){
    return (c1.sel_prob < c2.sel_prob);
}
void Problem::Verify(){
    for (int q = 0; q < LongestPath.size(); q++)
        for (int i = 0; i < LongestPath[q].size(); i++){
            if (LongestPath[q][0].first > Pairs - 1){
                bool bk = false;
                for (int j = i + 1; j < LongestPath[q].size(); j++){
                    if (LongestPath[q][i].first == LongestPath[q][j].first){
                        LongestPath.erase(LongestPath.begin() + q);
                        q--;
                        bk = true;
                        break;
                    }
                }
                if (bk == true) break;
            }
        }
}
