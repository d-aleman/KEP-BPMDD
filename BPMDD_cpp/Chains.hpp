//
//  Chains.hpp
//  LagrangeBranch&Bound
//
//  Created by Carolina Riascos Alvarez on 2020-06-11.
//  Copyright © 2020 Carolina Riascos Alvarez. All rights reserved.
//

#ifndef Chains_hpp
#define Chains_hpp

#include <stdio.h>
#include "Class_Problem_VFS.hpp"
bool isInside(vector<int>chain, int node);
bool v2AlreadyinChain(vector<vChain> v1, int v2);
bool PtivePriceCycle(IloNumArray& solpi, map<pair<int,int>,double>&Weights, vector<int>& cycle);
int FindNewRow(int& posCol, int NewCol, IloNumArray2& p_sol, IloNumArray2& AdjacencyList);
void MarkArc(vector<CheckedArc>&ArcsinSol, vector<int>&cycle, int&posCol, int&start);
void FindNewNeighbor(vector<Chain>& PPChains);

#endif /* Chains_hpp */

