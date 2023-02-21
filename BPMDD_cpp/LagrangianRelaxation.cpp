//
//  LagrangianRelaxation.cpp
//
//  Created by Carolina Riascos Alvarez on 2019-12-18.
//  Copyright © 2020 Carolina Riascos Alvarez. All rights reserved.
//

#include "LagrangianRelaxation.hpp"
#include <algorithm>

void Problem::Lagrange(){
    //Create model
    IloEnv env;
    IloModel Lag(env);
    IloCplex cplex(Lag);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::RootAlg, 1);
    cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
    //cplex.setParam(IloCplex::PreInd, 0);
    cplex.setParam(IloCplex::Param::TimeLimit, 1800);
    cplex.setParam(IloCplex::Param::Threads, 1);
    clock_t tStartMP = clock();
    tStartMP = clock();
    bool FirstTime = true;
    
    //Create variables
    IloNumVar y(env, 0, 1, ILOFLOAT);
    IloNumVarArray z(env, 0, 0, 1, ILOFLOAT);//AllCycles.size()
    
    //Create variables
    int pimany = 0;
    if (ChainLength == 0){
        pimany = int(Pairs);
    }
    else{
        pimany = int(Nodes);
    }
    
    //Arrays for constraints and objective
    IloRangeArray onecycle(env);
    IloObjective Obj(env);
    //Objective
    Obj = IloAdd(Lag, IloMaximize(env,-10000*y));
    for (int i = 0; i < pimany; i++){
       string name = "c." + to_string(i);
       const char* cName = name.c_str();
       onecycle.add(IloRange(env, -IloInfinity, y,1, cName));
    }
    Lag.add(onecycle);
    
    //Cutting plane
    IloBool NewCycleAdded = true;
    IloBool WarmUp = true;
    IloRangeArray cycles(env);
    solpi = IloNumArray(this->env,pimany);
    MPTime += (clock() - tStartMP)/double(CLOCKS_PER_SEC);
    int MIPiter = 0;
    vector<int>ColumnsAdded;
    vector<vector<int>>LastCol;
    ThirdPhase = false;
    
    while(NewCycleAdded == true){

        cplex.setOut(env.getNullStream());
        clock_t tStartMP2 = clock();
        tStartMP2 = clock();
        cplex.solve();
        MIPiter++;
        MPTime += (clock() - tStartMP2)/double(CLOCKS_PER_SEC);
        
        if (cplex.getStatus() == IloAlgorithm::Infeasible) {
            cout << "Infeasible";
        }
        else{
            //Check running time
            IloNum trun = (clock() - StartSolTime)/double(CLOCKS_PER_SEC);
            if (trun >= TimeLimit){
                CGNotOptimal = true;
                break;
            }
            //Get objective value
            double pastUpperBound = UpperBound;
            UpperBound = cplex.getObjValue();
            cout << "Upper Bound:" << UpperBound << endl;
            
            //Get dual variables
            cplex.getDuals(solpi,onecycle);
            Order.clear();
            for (int i = 0; i < vMDD.size(); i++){
                Order.push_back(make_pair(i, solpi[vMDD[i].first[1].actual]));
            }
            //Sort Order
            sort(Order.begin(), Order.end(), sortOrder);

            NewCycleAdded = false;
            bool next = false;
            int counter = 0;
            countSame = 0;
            
            int nCHCols = 0;
            if (WarmUp == true){
                for (int i = 0; i < Order.size(); i++){
                    posMDD = Order[i].first;
                    next = false;
                    if (next == false){
                        //Solve subproblems
                        clock_t tStartSP = clock();
                        LongestPath.clear();
                        FindLongestPath();
                        SPTime += (clock() - tStartSP)/double(CLOCKS_PER_SEC);
                        IloExpr NewCut(env,0);
                        for (int t = 0; t < LongestPath.size(); t++){
                            if (LongestPath[t].back().second > 0.005){//positive-price column
                                next = false;
                                if (pastUpperBound == UpperBound){
                                    countSame++;
                                    //check whether it is the same col
                                    next = false;
                                    for (int b = 0; b < LastCol.size(); b++){
                                        next = true;
                                        int sizec = 0;
                                        if (vMDD[posMDD].first[1].actual > Pairs - 1){
                                            sizec = int(LastCol[b].size() - 1);
                                        }
                                        else{
                                            sizec = int(LastCol[b].size());
                                        }
                                        if (LongestPath[t].size() == sizec && posMDD == LastCol[b].back()){
                                            for (int b2 = 0; b2 < LastCol[b].size() - 1; b2++){
                                                if (LastCol[b][b2] != LongestPath[t][b2].first){
                                                    next = false;
                                                    break;
                                                }
                                            }
                                            if (next == true){
                                                break;
                                            }
                                        }
                                        else{
                                            next = false;
                                        }
                                    }
                                }
                                else{
                                    LastCol.clear();
                                    countSame = 0;
                                }
                                if (vMDD[posMDD].first[1].actual > Pairs - 1 && Preference == "CY" && nCHCols > 0) next = true;
                                if (next == false){
                                    if (vMDD[posMDD].first[1].actual > Pairs - 1){
                                      nColChains++;
                                    }
                                    else{
                                      nColCycles++;
                                    }
                                    AllCycles.push_back(vector<int>());
                                    float weight = 0;
                                    //Create new column
                                    IloNumColumn col(env);
                                    for (int j = 0; j < LongestPath[t].size() - 1; j++){
                                        int nodei = LongestPath[t][j].first;
                                        int nodej = LongestPath[t][j +  1].first;
                                        //Update cycle constraints
                                        col+= onecycle[nodei](1);
                                        //Store node in new cycle
                                        AllCycles.back().push_back(nodei);
                                        NodeinCycle[nodei].push_back(int(AllCycles.size() - 1));
                                        if (vMDD[posMDD].first[1].actual > Pairs - 1 && j == LongestPath[t].size() - 2){
                                            col+= onecycle[nodej](1);
                                            AllCycles.back().push_back(nodej);
                                            NodeinCycle[nodej].push_back(int(AllCycles.size() - 1));
                                        }
                                        weight += Weights[make_pair(nodei,nodej)];
                                    }
                                    AllCyclesWeight.push_back(weight);
                                    AllCycles.back().push_back(posMDD); //Add copy the cycle belongs to
                                    if (pastUpperBound == UpperBound) LastCol.push_back(AllCycles.back());
                                    NewCycleAdded = true;
                                    //Add weight of new column
                                    col += Obj(weight);
                                    //Create new variable
                                    z.add(IloNumVar(col));
                                    string name = "z." + to_string(AllCycles.size());
                                    const char* varName = name.c_str();
                                    z[z.getSize() - 1].setName(varName);
                                    z[z.getSize() - 1].setBounds(0, 1);
                                    col.end();
                                    counter++;
                                    if (AllCycles.size() > 6000 && ChainLength > 0) {
                                        break;
                                    }
                                    if (vMDD[posMDD].first[1].actual > Pairs - 1 && Preference == "CY"){
                                        nCHCols++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            if (NewCycleAdded == true){
                FirstTime = false;
            }
            else{
                UB_PH1 = UpperBound;
                //Chains//
                if (ChainLength == 0 || NDDs == 0) ThirdPhase = true;
                WarmUp = false;
                vector<vector<int>>vChains;
                while(true){
                    if (ThirdPhase == false){
                        clock_t tStartSP = clock();
                        vChains = Chains();
                        TimePH2 += (clock() - tStartSP)/double(CLOCKS_PER_SEC);
                        if (Preference == "CY"){
                            nCHCols = 0;
                            for (int i = 0; i < vChains.size(); i++){
                                if (vChains[i][0] > Pairs - 1) nCHCols++;
                                if (nCHCols > 1 && vChains[i][0] > Pairs - 1){
                                    vChains.erase(vChains.begin() + i);
                                    i--;
                                }
                            }
                        }
                    }
                    else{
                        if (CycleLength > 3 || ChainLength > 0){
                            clock_t tStartSP = clock();
                            vChains = LeftCycles();
                            TimePH3 += (clock() - tStartSP)/double(CLOCKS_PER_SEC);
                        }
                    }
                    float weight = 0;
                    countSame = 0;
                    for (int i = 0; i < vChains.size(); i++){
                        next = false;
                        if (pastUpperBound == UpperBound){
                            countSame++;
                            //check whether it is the same col
                            next = false;
                            for (int b = 0; b < LastCol.size(); b++){
                                next = true;
                                int sizec = 0;
                                if (vChains[i][0] > Pairs - 1){
                                    sizec = int(LastCol[b].size() - 1);
                                }
                                else{
                                    sizec = int(LastCol[b].size());
                                }
                                if (vChains[i].size() == sizec){
                                    for (int b2 = 0; b2 < LastCol[b].size() - 1; b2++){
                                        if (LastCol[b][b2] != vChains[i][b2]){
                                            next = false;
                                            break;
                                        }
                                    }
                                    if (next == true){
                                        break;
                                    }
                                }
                                else{
                                    next = false;
                                }
                            }
                        }
                        else{
                            LastCol.clear();
                            countSame  = 0;
                        }
                        if (next == false){
                            if (vChains[i][0] > Pairs - 1){
                               nColChains++;
                             }
                             else{
                               nColCycles++;
                             }
                            AllCycles.push_back(vector<int>());
                            weight = 0;
                            //Create new column
                            IloNumColumn col(env);
                            for (int j = 0; j < vChains[i].size() - 1; j++){
                                int nodei = vChains[i][j];
                                int nodej = vChains[i][j +  1];
                                //Update cycle constraints
                                col+= onecycle[nodei](1);
                                //Store node in new cycle
                                AllCycles.back().push_back(nodei);
                                NodeinCycle[nodei].push_back(int(AllCycles.size() - 1));
                                if (vChains[i][0] > Pairs - 1 && j == vChains[i].size() - 2){
                                    col+= onecycle[nodej](1);
                                    AllCycles.back().push_back(nodej);
                                    NodeinCycle[nodej].push_back(int(AllCycles.size() - 1));
                                }
                                weight += Weights[make_pair(nodei,nodej)];
                            }
                            AllCyclesWeight.push_back(weight);
                            AllCycles.back().push_back(0); //Add copy the cycle belongs to
                            if (pastUpperBound == UpperBound) LastCol.push_back(AllCycles.back());
                            NewCycleAdded = true;
                            //Add weight of new column
                            col += Obj(weight);
                            //Create new variable
                            z.add(IloNumVar(col));
                            string name = "z." + to_string(AllCycles.size());
                            const char* varName = name.c_str();
                            z[z.getSize() - 1].setName(varName);
                            z[z.getSize() - 1].setBounds(0, 1);
                            //col.end();
                        }
                    }
                    if (NewCycleAdded == false){
                        if (ThirdPhase == true){
                            break;
                        }
                        ThirdPhase = true;
                    }
                    else{
                        if (ChainLength > 0) ThirdPhase = false;
                        break;
                    }
                }
            }
        }
    }
    //Round Upper Bound down to closest integer if all arc weights are integer
    if (AllArcWeightsInt == true) UpperBound = floor (UpperBound + 0.1);
    
    
    //Get feasible solution
    clock_t tStartCF = clock();
    exploredBBnodes++;

    cplex.getValues(z_sol,z);
    ListCycles.clear();
    ListCycles = vector<pair<int,vector<cycle>>>(vMDD.size() + NDDs);
    int count = 0;
    for (int i = 0; i < z_sol.getSize(); i++){
        //cout << endl;
        if (z_sol[i] > 0.001){
            //cout << "Cycle " << i << ": " << z_sol[i] << endl;
            if (z_sol[i] == 1) count++;
            vector<int>vc;
            //cout << endl;
            for (int j = 0; j < AllCycles[i].size() - 1; j++){
                vc.push_back(AllCycles[i][j]);
                //cout << AllCycles[i][j] << "\t";
            }
            ListCycles[AllCycles[i].back()].first = AllCycles[i].back();
            ListCycles[AllCycles[i].back()].second.push_back(cycle(vc,z_sol[i],i));
            ListCycles[AllCycles[i].back()].second.back().TotalWeight = AllCyclesWeight[i];
        }
    }
    
    for (auto it = ListCycles.begin(); it != ListCycles.end();){
        if (it->second.size() == 0){
            it = ListCycles.erase(it);
        }
        else{
            it++;
        }
    }
    
    string name = "LogCF_K" + to_string(CycleLength) + "_L" + to_string(ChainLength) + "_" + Preference + "_" + FileName;
    const char* cName = name.c_str();
    ofstream logfile(cName, ios_base::app);
    cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 1);
    cplex.setParam(IloCplex::RootAlg, 0);
    double subTime = TimeLimit - SPTime - MPTime;

    cplex.setParam(IloCplex::Param::TimeLimit, 400);
    cplex.setOut(logfile);
    //Transform variables into integers
    Lag.add(IloConversion(env, z, ILOBOOL));
    //Set start MIP solution
    cplex.solve();
    CFTime += (clock() - tStartCF)/double(CLOCKS_PER_SEC);
    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        env.out() << "No solution" << endl;
    }
    else {
        CFObj = cplex.getObjValue();
        IloNumArray z_sol(this->env);
        cplex.getValues(z_sol,z);
        
        if (CFObj > GlobalLB){
            vFeasSol.clear(); vFeasSol.end();
            for (int i = 0; i < AllCycles.size(); i++){
                if (z_sol[i] > 0.9){
                    vFeasSol.push_back(vector<int>());
                    for (int j = 0; j < AllCycles[i].size() - 1; j++){
                        vFeasSol.back().push_back(AllCycles[i][j]);
                    }
                }
            }
        }
    }
    
    AllCycles.clear();
    AllCyclesWeight.clear();
    Lag.end();
    cplex.end();
    env.end();
    
    //Update time
    IloNum trun = (clock() - StartSolTime)/double(CLOCKS_PER_SEC);
    if (trun >= TimeLimit){
        TimedOut = true;
    }
}

