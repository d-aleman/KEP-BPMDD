//
//  HeaderCycle.hpp
//  CycleFormulation
//
//  Created by Carolina Riascos Alvarez on 2019-11-27.
//  Copyright © 2019 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef mainCycle_hpp
#define mainCycle_hpp

#include <iostream>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <string>
#include "ReadData.hpp"
#include "FindCycles.hpp"
#include "Class_Problem.hpp"

ILOSTLBEGIN

IloEnv env;
IloModel model;

//Structures for reading
IloInt Nodes, NDDs, Pairs, numMatch, NumArcs;
IloNumArray2 WeightMatrix(env), AdjacencyList(env);
IloNumArray3 mCopies(env);
vector<vector<vector<int>>> digraphD;
IloNumArray vikingos(env);

void Problem::llenarComp(){ //dani  map<int,vector<int>>Comp;//Successors  map<int,vector<int>>CompPred;//predecessor
    CompPred.clear();
    Comp.clear();
    for (int i = 0; i < AdjacencyList.getSize(); i++) {
        int k = int(AdjacencyList[i].getSize());
        if (k == 0) {
            Comp[i].clear();
            CompPred[i].clear();
        }
        for (int j =0; j < AdjacencyList[i].getSize(); j++) {
            auto aux = AdjacencyList[i][j] - 1;
            Comp[i].push_back(aux);
            CompPred[i];
            CompPred[aux].push_back(i);
            Comp[aux];
        }
    }
    for (auto i = 0; i < Nodes; i++) {
        if (Comp[i].size() == 0 && CompPred[i].size() == 0){
            auto it1 = Comp.find(i);
            Comp.erase(it1);
            auto it2 = CompPred.find(i);
            CompPred.erase(it2);
        }
    }
}
void cuantoTarda(clock_t tStart){
    double d = (clock() - tStart)/double(CLOCKS_PER_SEC);
    if (d <= 6000) {
        cout << d << " s." << endl;
    }else{
        cout << d/double(60) << " min." << endl;
    }
}
bool compaDegree(pair<int, int>& p1, pair<int, int>& p2){ return p1.second > p2.second;}
vector<pair<int,int>> Problem:: ordenDegre(){
    vector<pair<int,int>> vp;
    if (DegreeType == "Indegree"){
        for (auto it: CompPred){
            vp.push_back(make_pair(it.first, it.second.size()));
        }
    }else if (DegreeType == "Outdegree") {
        for (auto it: Comp){
            vp.push_back(make_pair(it.first, it.second.size()));
        }
    }else if (DegreeType == "Totaldegree") {
        auto it2 = CompPred.begin();
        for (auto it = Comp.begin(); it != Comp.end(); it++){
            vp.push_back(make_pair(it->first , it->second.size() + it2->second.size()));
            it2++;
        }
    } else if (DegreeType == "Increasing"){
        for (auto i = 0; i < Nodes; i++) {
            if (Comp.count(i) == 1){
                vp.push_back(make_pair(i, i));
            }
        }
        return vp;
    }
    else if (DegreeType == "OptimalVFS"){
        MainCycleFinder();
        vp = OptimalVFS();
        return vp;
    }
    sort(vp.begin(), vp.end(), compaDegree);
    return vp;
}
#endif /* HeaderCycle_hpp */
