//
//  BBTree.cpp
//  LagrangeBranch&Bound
//
//  Created by Carolina Riascos Alvarez on 2020-04-03.
//  Copyright © 2020 Carolina Riascos Alvarez. All rights reserved.
//

#include "BBTree.hpp"
void Problem::getChildren(){
    
    //Verify if solution is feasible
    LagArcs.clear();
    LagNodes.clear();
    
    //Store arc information
    vector<pair<pair<int,int>,arc>>arcos; ///<i,j>, vector<copy1,...,copyn>
    map<pair<int,int>, arc>ListArcs;
    pair<int,vector<pair<int,int>>>RemovedArcs;
    map<int,vector<int>>nodes; //node, list of copies
    map<int, float> Fraction;
    vector<pair<int,map<int,int>>>NodeCountInCopy;//copy, nodes in copy
    bool IsFractionalSol = false;
    
    for (int i = 0; i < ListCycles.size(); i++){
        NodeCountInCopy.push_back(make_pair(i, map<int,int>()));
        for (int j = 0; j < ListCycles[i].second.size(); j++){
            for (int h = 0; h < ListCycles[i].second[j].acycle.size(); h++){
                int node = ListCycles[i].second[j].acycle[h];
                NodeCountInCopy.back().second[node]++;
                float div = ListCycles[i].second[j].fraction;
                Fraction[node]+=div;
                pair<int,int>p;
                if (h < ListCycles[i].second[j].acycle.size() - 1){
                    p = make_pair(ListCycles[i].second[j].acycle[h], ListCycles[i].second[j].acycle[h+1]);
                }else{
                    if (ListCycles[i].second[j].acycle[0] < Pairs){
                        p = make_pair(ListCycles[i].second[j].acycle[h], ListCycles[i].second[j].acycle[0]);
                    }
                }
                ListArcs[p].fraction+= ListCycles[i].second[j].fraction;
                ListArcs[p].WhereinMDDs.push_back(ListCycles[i].first); //last position is the copy where the cycle is present
            }
        }
    }
    ListCycles.clear();
    
    //Clasify arcs
    LagArcs.clear();
    for (auto it = ListArcs.begin(); it != ListArcs.end(); it++){
        LagArcs.push_back(make_pair(it->first, it->second));
    }
    sort(LagArcs.begin(), LagArcs.end(), sortArcsFrac);//Sort arcs closer to 0.5
    ListArcs.clear();
    ListArcs.end();
    
    //Classify nodes
    for (int i = 0; i < NodeCountInCopy.size(); i++){
        for (auto it = NodeCountInCopy[i].second.begin(); it != NodeCountInCopy[i].second.end(); it++){
            for (int j = 0; j < it->second; j++){
                nodes[it->first].push_back(i);
            }
        }
    }
    for (auto it = nodes.begin(); it != nodes.end(); it++){
        LagNodes.push_back(node(it->first, Fraction[it->first], it->second));
    }
    sort(LagNodes.begin(), LagNodes.end(), sortNodesRep);
    
    int PosNodeToBranch = -1;
    //Check whether a node is repeated among copies
    if (LagNodes[0].WhereinMDDs.size() > 1){
        IsFractionalSol = true;
        if (BranchingMethod == 1){//By Node
            //Find node to branch on
            while (PosNodeToBranch == -1){
                for (int i = 0; i < LagNodes.size(); i++){
                    int reference = LagNodes[0].WhereinMDDs[0];
                    if (LagNodes[i].WhereinMDDs.size() > 1){
                        bool unique = true;
                        for (int j = 1; j < LagNodes[i].WhereinMDDs.size(); j++){//From 1 to compare the other copies
                            if (LagNodes[i].WhereinMDDs[j] != reference){
                                unique = false;
                                break;
                            }
                        }
                        if (unique == false){
                            PosNodeToBranch = i;
                            break;
                        }
                    }else{
                        PosNodeToBranch = 0;
                        break;
                    }
                }
            }
        }
    }
    
    nodes.end();
    ListCycles.end();
    
    if (IsFractionalSol == true){
        if (UpperBound > GlobalLB){//BRANCH ONLY IF NEED BE
            if (FirstBBnode == true){
                VoyBB = make_pair(0,0);
                Tree.push_back(vector<BBnode>());
                Tree.back().push_back(BBnode(0, 0, -1, UpperBound, vector<pair<int,int>>(), vector<pair<int,int>>(), make_pair(-1,-1)));
                FirstBBnode = false;
            }
            //Create new level if need be
            if (VoyBB.first + 1 > Tree.size() - 1) Tree.push_back(vector<BBnode>());
            //pair<int,int>parent = Tree[VoyBB.first][VoyBB.first].parentNode;
            if (BranchingMethod == 1){//By nodes
                for (int i = 0; i < LagNodes[PosNodeToBranch].WhereinMDDs.size(); i++){
                    map<int,vector<int>>SetNodesZero;
                    //include parent's arcs map ZERO
                    if (VoyBB.first != 0 || VoyBB.second != 0){
                        SetNodesZero = Tree[VoyBB.first][VoyBB.second].NodeGroupZero;
                    }
                    vector<int>SetNodesOne;
                    //include parent's arcs map ONE
                    if (VoyBB.first != 0 || VoyBB.second != 0){
                        SetNodesOne = Tree[VoyBB.first][VoyBB.second].NodeGroupOne;
                    }
                    ////////////////////Create child/////////////////
                    SetNodesOne.push_back(LagNodes[PosNodeToBranch].id);
                    for (int j = 0; j < vMDD.size(); j++){
                        if (LagNodes[PosNodeToBranch].WhereinMDDs[i] != j){
                            SetNodesZero[j].push_back(LagNodes[PosNodeToBranch].id);
                        }
                    }
                    Tree.back().push_back(BBnode(int(Tree.back().size() - 1), int(Tree.size() - 1), -1, -1000, SetNodesZero, SetNodesOne, VoyBB));
                    Tree[VoyBB.first][VoyBB.second].NeighNextLevel.push_back(int(Tree.back().size() - 1));
                }
                //Last B&B node
                map<int,vector<int>>SetNodesZero;
                if (VoyBB.first != 0 || VoyBB.second != 0){
                    SetNodesZero = Tree[VoyBB.first][VoyBB.second].NodeGroupZero;
                }
                for (int i = 0; i < LagNodes[PosNodeToBranch].WhereinMDDs.size(); i++){
                    for (int j = 0; j < vMDD.size(); j++){
                        if (LagNodes[PosNodeToBranch].WhereinMDDs[i] == j){
                            SetNodesZero[j].push_back(LagNodes[PosNodeToBranch].id);
                        }
                    }
                }
                Tree.back().push_back(BBnode(int(Tree.back().size() - 1), int(Tree.size() - 1), -1, -1000, SetNodesZero, Tree[VoyBB.first][VoyBB.second].NodeGroupOne, VoyBB));
                Tree[VoyBB.first][VoyBB.second].NeighNextLevel.push_back(int(Tree.back().size() - 1));
            }
            else{//By arcs
                ArcBranches.push_back(LagArcs[0].first);
                vector<pair<int,int>>SetArcsZero;
                
                //include parent's arcs map ZERO
                if (VoyBB.first != 0 || VoyBB.second != 0){
                    SetArcsZero = Tree[VoyBB.first][VoyBB.second].ArcGroupZero;
                }
                vector<pair<int,int>>SetArcsOne;
                //include parent's arcs map ONE
                if (VoyBB.first != 0 || VoyBB.second != 0){
                    SetArcsOne = Tree[VoyBB.first][VoyBB.second].ArcGroupOne;
                }
                for (int i = 0; i < 2; i ++){
                    //Create child
                    if (i == 0){
                        SetArcsZero.push_back(LagArcs[0].first);
                        Tree.back().push_back(BBnode(int(Tree.back().size() - 1), int(Tree.size() - 1), -1, -1000, SetArcsZero, Tree[VoyBB.first][VoyBB.second].ArcGroupOne, VoyBB));
                    }
                    else{
                        SetArcsOne.push_back(LagArcs[0].first);
                        Tree.back().push_back(BBnode(int(Tree.back().size() - 1), int(Tree.size() - 1), -1, -1000, Tree[VoyBB.first][VoyBB.second].ArcGroupZero, SetArcsOne, VoyBB));
                    }
                    Tree[VoyBB.first][VoyBB.second].NeighNextLevel.push_back(int(Tree.back().size() - 1));
               }
            }
        }
    }
    else{/////////////////////INTEGER SOLUTION//////////////////
        Tree[VoyBB.first][VoyBB.second].LB = UpperBound;
    }
};

void Problem::BBTree(){
    
    //Modify Adjancency List: Super Source node
    AdjacencyList.add(IloNumArray(env,0));
    for (int i = Pairs; i < Nodes; i++){
        AdjacencyList[AdjacencyList.getSize() - 1].add(i + 1);
        Weights[make_pair(AdjacencyList.getSize() - 1, i)] = 0;
        PredList[i].push_back(AdjacencyList.getSize() - 1);
    }
    
    //Sink node
    AdjacencyList.add(IloNumArray(env,0));
    PredList.push_back(vector<int>());//source
    PredList.push_back(vector<int>());//sink
    for (int i = 0; i < Pairs; i++){
        AdjacencyList[i].add(AdjacencyList.getSize());
        Weights[make_pair(i, AdjacencyList.getSize() - 1)] = 0;
        //cout << AdjacencyList.getSize() - 1 << " " << PredList.size() << endl;
        PredList[AdjacencyList.getSize() - 1].push_back(i);
    }
    
    
    StartSolTime = clock();
    //Lagranagian relaxation
    VoyBB = make_pair(0, 0);
    Lagrange();
    
    if (TimedOut == false){
        GlobalUB = UpperBound;
        GlobalLB = CFObj;
        if (CFObj < UpperBound - 0.00001){
            //Save Solution
            vBestSol = vFeasSol;
            //Get children
            getChildren();
            Tree[VoyBB.first][VoyBB.second].LB = GlobalLB;
            insideBB = true;
            while(true){
                //Next neighbor
                if (Tree[VoyBB.first][VoyBB.second].nextNeighbor + 1 < Tree[VoyBB.first][VoyBB.second].NeighNextLevel.size()){//Still nodes to explore in the current B&B Node
                    Tree[VoyBB.first][VoyBB.second].nextNeighbor++;
                    VoyBB = make_pair(Tree[VoyBB.first][VoyBB.second].level + 1, Tree[VoyBB.first][VoyBB.second].NeighNextLevel[Tree[VoyBB.first][VoyBB.second].nextNeighbor]);
                    //Lagranagian relaxation
                    Lagrange();
                    //Update UB
                    Tree[VoyBB.first][VoyBB.second].UB = UpperBound;
                    Tree[VoyBB.first][VoyBB.second].LB = CFObj;
                    if (Tree[VoyBB.first][VoyBB.second].LB > GlobalLB){
                        GlobalLB = Tree[VoyBB.first][VoyBB.second].LB;
                        //Save solution
                        vBestSol = vFeasSol;
                        //cout << "We have found the best solution so far";
                        if (GlobalLB == GlobalUB) break;
                    }
                    if (TimedOut == true) {
                        break;
                    }
                    getChildren();
                    //Update node explored
                    VoyBB = Tree[VoyBB.first][VoyBB.second].parentNode;
                }
                else{
                    //Go to the next node in the same level
                    if (VoyBB.second + 1 < Tree[VoyBB.first].size()){
                        //Is there a node with a better upper bound than the best lower bound found so far?
                        bool better = false;
                        for (int h = VoyBB.second + 1; h < Tree[VoyBB.first].size(); h++){
                            if(Tree[VoyBB.first][h].UB > GlobalLB) {
                                better = true;
                                VoyBB = make_pair(VoyBB.first, h);
                                break;
                            }
                        }
                        if (better == false){
                            //Go down one level
                            if (VoyBB.first + 1 < Tree.size()){
                                //Set new level
                                VoyBB = make_pair(VoyBB.first + 1, 0);
                                //Update GlobalUB
                                double NewGlobalUB = 0;
                                for (int h = 0; h < Tree[VoyBB.first].size(); h++){
                                    if (Tree[VoyBB.first][h].UB > NewGlobalUB) NewGlobalUB = Tree[VoyBB.first][h].UB;
                                }
                                GlobalUB = NewGlobalUB;
                            }
                            else{//Nothing else to look for!
                                break;
                            }
                        }
                    }
                    else{//Go down one level
                        if (VoyBB.first + 1 < Tree.size()){
                            //Set new level
                            VoyBB = make_pair(VoyBB.first + 1, 0);
                            //Update GlobalUB
                            double NewGlobalUB = 0;
                            for (int h = 0; h < Tree[VoyBB.first].size(); h++){
                                if (Tree[VoyBB.first][h].UB > NewGlobalUB) NewGlobalUB = Tree[VoyBB.first][h].UB;
                            }
                            GlobalUB = NewGlobalUB;
                        }
                        else{
                            break;
                        }
                    }
                }
            }
             //Print best solution
            if (TimedOut == true){
                //vBestSol = vFeasSol;
                TBBtime = (clock() - StartSolTime)/double(CLOCKS_PER_SEC);
                cout << endl << "Feasible Solution " << endl << "Objective Value: " << GlobalLB << endl;
                PrintLag("Feasible");
            }
            else{
                //vBestSol = vFeasSol;
                TBBtime = (clock() - StartSolTime)/double(CLOCKS_PER_SEC);
                cout << endl << "Optimal Solution " << endl << "Objective Value: " << GlobalLB << endl;
                PrintLag("Optimal");
            }
            
        }
        else{
            //Save solution
            vBestSol = vFeasSol;
            TBBtime = (clock() - StartSolTime)/double(CLOCKS_PER_SEC);
            if (DegreeType == "Best_K-VFS" && TooManyVFSCycles == true){
                cout << endl << "BP_MDD unable to start...More than 4M cycles, thus no K-VFS was obtained." << endl;
                PrintLag("Unstarted");
            }else{
                //Print optimal solution
                cout << endl << "Optimal Solution " << endl << "Objective Value: " << GlobalLB << endl;
                PrintLag("Optimal");
            }
        }
    }
    else{
        //Save solution
        GlobalUB = UpperBound;
        GlobalLB = CFObj;
        vBestSol = vFeasSol;
        TBBtime = (clock() - StartSolTime)/double(CLOCKS_PER_SEC);
        //Print best solution
        cout << endl << "Feasible Solution " << endl << "Objective Value: " << GlobalLB << endl;
        if (CGNotOptimal == false){
            if (GlobalUB == GlobalLB){
                PrintLag("Optimal");
            }
            else{
                PrintLag("Feasible");
            }
        }
        else{
            PrintLag("Unknown_UB");
        }
    }
}
