//
//  CycleFormulation.cpp
//
//  Created by Carolina Riascos Alvarez on 2019-12-03.
//  Copyright © 2019 Carolina Riascos Alvarez. All rights reserved.
//

#include "FindCycles.hpp"


void SetName2(IloNumVar& var, const char* prefix, IloInt i, IloInt j){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i) + "." + to_string(j);
    const char* varName = name.c_str();
    var.setName(varName);
}
vector<pair<int,int>> Problem::OptimalVFS(){
    vector<pair<int,int>>OptVFS;
    IloModel vfs(env);
    IloCplex cp_vfs(vfs);
    cp_vfs.setParam(IloCplex::TiLim, 1800);
    IloNum tStartVFS = cp_vfs.getTime();
    
    //Create variables
    IloNumVarArray d(env, Pairs, 0, 1, ILOINT);
    for (int i = 0; i < Pairs; i++){
        SetName(d[i], "d", i + 1);
        //cout << d[i].getName() << endl;
    }
    
    IloRangeArray Covering(env, AllCyclesVFS.size());
    for (int i = 0; i < AllCyclesVFS.size(); i++){
        IloExpr vertices (env, 0);
        for (int j = 0; j < AllCyclesVFS[i].size(); j++){
            vertices+= d[AllCyclesVFS[i][j]];
        }
        Covering[i] = IloRange(env, 1, vertices, IloInfinity);
        vertices.end();
    }
    vfs.add(Covering);
    
    IloObjective Obj (env, IloSum(d), IloObjective::Minimize);
    vfs.add(Obj);
    
    
    //Solve
    cp_vfs.setOut(env.getNullStream());
    cp_vfs.solve();
    Status_VFS_IP = cp_vfs.getStatus();
    cout << "StatusVFS: " << Status_VFS_IP << endl;
    Size_K_VFS = cp_vfs.getObjValue();
    IloNumArray Sol_vfs(env);
    cp_vfs.getValues(Sol_vfs, d);
    
    for (int i = 0; i < Sol_vfs.getSize(); i++){
        if (Sol_vfs[i] > 0.9) OptVFS.push_back(make_pair(i,i));
    }
    
    Time_VFS_IP = cp_vfs.getTime() - tStartVFS;
    cout << "VFS-IP time: " << Time_VFS_IP << endl;
    vfs.end();
    cp_vfs.end();
    d.end();
    Obj.end();
    Covering.end();
    
    return OptVFS;
}
void Problem::MainCycleFinder(){
    cout << "Finding K-Cycles..." << endl;
    //Pick a source node
    int origin;
    
    //Copy Adjacency List
    IloNumArray2 AdjaList(env, AdjacencyList.getSize());
    for (int l = 0; l < AdjacencyList.getSize(); l++){
        AdjaList[l] = IloNumArray(env);
        for (int i = 0; i < AdjacencyList[l].getSize(); i++){
            AdjaList[l].add(AdjacencyList[l][i]);
        }
    }
    AllCycles.clear();
    bool terminate = false;
    origin = -1;
    while (terminate == false){
        //Set origin
        int u = origin + 1;
        int prev_origin = origin;
        for (; u < Pairs; u++){
            if (AdjaList[u].getSize() > 0){
                origin = u;//subtracted 1, so 100 is 99
                break;
            }
        }
        if (origin == prev_origin){
            break;
        }
        
        //Find Cycles
        int ThisMany = int(ListCycles.size());
        vector<vector<int>> ListC;
        ListC = SubCycleFinder(env, AdjaList, origin);
        for (int i = 0; i < ListC.size(); i++){
            AllCyclesVFS.push_back(ListC[i]);
            if (AllCyclesVFS.size() > 4000000) TooManyVFSCycles = true;
            for (int y = 0; y < ListC[i].size(); y++){
                CycleNode[ListC[i][y]].push_back(ListCycles.size() - 1);
            }
        }
        if (TooManyVFSCycles == true) terminate = true;
        
        //Remove node
        for (int l = 0; l < AdjaList.getSize(); l++){
            if (l == origin){
                AdjaList[l] = IloNumArray(env,0);
            }
            else{
                for (int i = 0; i < AdjaList[l].getSize(); i++){
                    if (AdjaList[l][i] == origin + 1){
                        AdjaList[l].remove(i);
                    }
                }
            }
        }
        
    }
}
vector<vector<int>> Problem::SubCycleFinder (IloEnv env, IloNumArray2 AdjaList, IloInt origin){
    vector<vector<int>> ListC;
    vector<int> nodesBestCycle;
    IloNumArray stepinSU (env, AdjaList.getSize());
    for (int c = 0; c < stepinSU.getSize(); c++) stepinSU[c] = 0;
    vector<int> Stack;
    Stack.push_back(origin);
    IloInt whoi = origin;
    stepinSU[origin] = 0;
    IloBool VLpresent = true; //because we start with the origin
    IloBool cont = true;
    IloInt prevWhoi = 1000000;
    IloNum BigW = 0;
    IloNum weight = 0;
    IloBool NewCycle = false;
    if (AdjaList[whoi].getSize() == 0){
        cont = false;
    }

    while (cont == true){
        if (IsxinStack(AdjaList[whoi][stepinSU[whoi]] - 1, Stack) == true || Stack.size() >= CycleLength){
            if (AdjaList[whoi][stepinSU[whoi]] == (origin + 1) && VLpresent == true && NewCycle == false){
                // they form a cycle together!
                weight = PathWeight(Stack);
                ListC.push_back(Stack);
                if (weight > BigW){
                    BigW = weight;
                    nodesBestCycle = Stack;
                }
                //Backtrack
                NewCycle = true;
            }
            else{
                NewCycle = false;
                while(true){
                    //Backtrack
                    stepinSU[whoi]++;
                    if (stepinSU[whoi] >= AdjaList[whoi].getSize()){
                        if (Stack.back() == origin){ VLpresent = false;}
                        stepinSU[whoi] = 0;
                        Stack.pop_back();
                        if (Stack.size() == 0){
                            cont  = false;
                            break;
                        }else{
                            whoi = Stack.back();
                        }
                    }
                    else{
                        break;
                    }
                }
            }
        }
        else{
            if (AdjaList[whoi][stepinSU[whoi]] == origin + 1) VLpresent = true;
            Stack.push_back(AdjaList[whoi][stepinSU[whoi]] - 1);
            prevWhoi = whoi;
            whoi = AdjaList[whoi][stepinSU[whoi]] - 1;
            if (AdjaList[whoi].getSize() == 0){
                //bactrack
                if (Stack.back() == origin){ VLpresent = false;}
                stepinSU[whoi] = 0;
                Stack.pop_back();
                whoi = Stack.back();
                stepinSU[whoi]++;
                while(true){
                    if (stepinSU[whoi] >= AdjaList[whoi].getSize()){
                        if (Stack.back() == origin){VLpresent = false;}
                        stepinSU[whoi] = 0;
                        Stack.pop_back();
                        if (Stack.size() == 0){
                            cont  = false;
                        }else{
                            whoi = Stack.back();
                            stepinSU[whoi]++;
                        }
                    }
                    else{
                        break;
                    }
                }
            }
        }
    }
    return ListC;
}
IloNum Problem::PathWeight (vector<int>& Stack){
    IloNum weight = 0;
    //arc weights
    for (int l = 0; l < Stack.size(); l++){
        if (l >= 1){
            for (int i = 0; i < PredList[Stack[l]].size(); i++){
                if (PredList[Stack[l]][i] == Stack[l - 1]){
                    weight += Weights[make_pair(Stack[l - 1], Stack[l])];
                    break;
                }
            }
        }
        else{
            for (int i = 0; i < PredList[Stack[l]].size(); i++){
                if (PredList[Stack[l]][i] == Stack[Stack.size() - 1]){
                    weight += Weights[make_pair(Stack[Stack.size() - 1], Stack[l])];
                    //weight += PredList[Stack[l]][i].second;
                    break;
                }
            }
        }
    }
    return weight;
}
IloBool Problem::IsxinStack (IloInt test, vector<int>& xinTrial){
    for (int i = 0; i < xinTrial.size(); i++){
        if (test == xinTrial[i])    return true;
    }
    return false;
}
void Problem::CycleFormulation(){
    //Find all arcs selected in the solution
    map<pair<int,int>,bool> Allarcs;
    for (int i = 0; i  < AllCycles.size(); i++){
        if (vMDD[AllCycles[i].back()].first[1].actual < Pairs){
            Allarcs[make_pair(AllCycles[i][AllCycles[i].size() - 2], AllCycles[i][0])] = true;
        }
        for (int j = 0; j  < AllCycles[i].size() - 2; j++){
            Allarcs[make_pair(AllCycles[i][j], AllCycles[i][j + 1])] = true;
        }
    }
    
    typedef IloArray<IloNumVarArray> NumVarArray2;
    
    //Create model
    IloModel cf(env);
    IloCplex cplex(cf);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Param::TimeLimit, 1800);
    cplex.setParam(IloCplex::Param::Threads, 1);
    //IloNum time1 = cplex.getCplexTime();
    
    //Create variables
    IloNumVarArray z(env, AllCyclesWeight.size(), 0, 1, ILOFLOAT);

    for (int i = 0; i  < AllCyclesWeight.size(); i++){
        SetName(z[i], "z", i);
    }
    
    
    int pimany = 0;
    if (ChainLength == 0){
        pimany = int(Pairs);
    }
    else{
        pimany = int(Nodes);
    }
    
    //Create variables
    NumVarArray2 x(env, AdjacencyList.getSize());
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        x[i] = IloNumVarArray(this->env, AdjacencyList[i].getSize(), 0, 1, ILOFLOAT);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            SetName2(x[i][j], "x", i + 1, AdjacencyList[i][j]);
            //cout << r[i][j].getName() << endl;
            if (Allarcs.count(make_pair(i,AdjacencyList[i][j] -1)) == true){
                //x[i][j].setUB(0);
            }
        }
    }
    
    //Constraints
    IloRangeArray assignment(env);
    for (int i = 0; i < pimany; i++){
        string name = "a." + to_string(i);
        const char* cName = name.c_str();
        IloExpr exp(env,0);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            if (Allarcs.count(make_pair(i,AdjacencyList[i][j] -1)) == false){
                exp+= x[i][j];
            }
        }
        assignment.add(IloRange(env, 0, exp,1, cName));
        exp.end();
    }
    cf.add(assignment);
    
    IloRangeArray flow(env);
    for (int i = 0; i < pimany; i++){
        IloExpr units(env,0);
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            units+= x[i][j];
        }
        for (int j = 0; j < PredList[i].size(); j++){
            int h = 0;
            for (h = 0; h < AdjacencyList[PredList[i][j]].getSize(); h++){
                if (AdjacencyList[PredList[i][j]][h] == i + 1) {
                   break;
                }
            }
            units-= x[PredList[i][j]][h];
        }
        cout << units << endl;
        flow.add(IloRange(env,-IloInfinity, units, 0));
        units.end();
    }
    
    //Objective value
    IloExpr Obj(env, 0);
    
    for (int i = 0; i < AdjacencyList.getSize(); i++){
        for (int j = 0; j < AdjacencyList[i].getSize(); j++){
            if (Allarcs.count(make_pair(i,AdjacencyList[i][j] -1)) == false){
                Obj += Weights[make_pair(i,AdjacencyList[i][j] -1)]*x[i][j];
            }
        }
    }
    //cout << Obj << endl;
    cf.add(IloMaximize(env,Obj));

    cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        env.out() << "No solution" << endl;
    }
    else {
        double Obj_B = cplex.getObjValue();
        
        env.out() << "Objective: " << Obj_B << endl;
        
        z_sol = IloNumArray (env, ListCycles.size());
        cplex.getValues(z_sol,z);
    }
    x.end();
    assignment.end();
    
}
void Problem::InterchangeableNodes(map<int,vector<int>>Comp, int origin){
    for (int i = 0; i < Pairs; i++){
        sort(Comp[i].begin(), Comp[i].end(), sortLowHigh);
    }
    int Threshold = 0;
    map<int,bool>selected;
    if (origin > Pairs - 1){//chain
        Threshold = int(ChainLength) - 1;
    }
    else{
        Threshold = int(CycleLength) - 2;
    }
    
    ListOfEquals.push_back(map<int,vector<vector<int>>>());
    //Check whether origin's neighbors have the same successors
    
    int round = 1;
    while(round <= Threshold){
        map<int,bool>NewToCheck;
        ListOfEquals.back()[round] = vector<vector<int>>();
        for (int i = 0; i < Comp[origin].size(); i++) {
            ListOfEquals.back()[round].push_back(vector<int>());
            ListOfEquals.back()[round].back().push_back(Comp[origin][i]);
            for (int j = i + 1; j < Comp[origin].size(); j++) {
                bool sameNeighbors = true;
                if (Comp[origin][i] != Comp[origin][j] && Comp[Comp[origin][i]].size() == Comp[Comp[origin][j]].size() && selected.count(Comp[origin][j]) == false){
                    for (int h = 0; h < Comp[Comp[origin][i]].size(); h++) {
                        if (Comp[Comp[origin][i]][h] != origin){
                            NewToCheck[Comp[Comp[origin][i]][h]] = true;
                        }
                        if (Comp[Comp[origin][i]][h] != Comp[Comp[origin][j]][h]){
                            sameNeighbors = false;
                            break;
                        }
                    }
                    if (sameNeighbors == true){
                        ListOfEquals.back()[round].back().push_back(Comp[origin][j]);
                        selected[Comp[origin][j]] = true;
                    }
                }
            }
            if (ListOfEquals.back()[round].back().size() == 1) ListOfEquals.back()[round].erase(ListOfEquals.back()[round].begin() + ListOfEquals.back()[round].size() - 1);
        }
        round++;
        Comp[origin].clear();
        for (auto it = NewToCheck.begin(); it != NewToCheck.end(); it++){
            Comp[origin].push_back(it->first);
        }
    }
    
}



