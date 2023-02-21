//
//  CycleFormulation.cpp
//  CycleFormulation
//
//  Created by Carolina Riascos Alvarez on 2019-12-03.
//  Copyright © 2019 Carolina Riascos Alvarez. All rights reserved.
//

#include "FindCycles.hpp"
#include <random>
random_device rd;
mt19937 generator(rd());

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
    cp_vfs.setParam(IloCplex::TiLim, 3600);
    
    
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
    cout << "StatusVFS: " << cp_vfs.getStatus() << endl;
    IloNumArray Sol_vfs(env);
    cp_vfs.getValues(Sol_vfs, d);
    
    for (int i = 0; i < Sol_vfs.getSize(); i++){
        if (Sol_vfs[i] > 0.9) OptVFS.push_back(make_pair(i,i));
    }
    
    vfs.end();
    cp_vfs.end();
    d.end();
    Obj.end();
    Covering.end();
    return OptVFS;
}
void Problem::MainCycleFinder(){
    //Pick a source node
    int origin;
    
    //for (int u = 0; u < comp.size(); u++){
        //Copy Adjacency List
    IloNumArray2 AdjaList(env, AdjacencyList.getSize());
    for (int l = 0; l < AdjacencyList.getSize(); l++){
        AdjaList[l] = IloNumArray(env);
        for (int i = 0; i < AdjacencyList[l].getSize(); i++){
            AdjaList[l].add(AdjacencyList[l][i]);
        }
        //cout << AdjaList[l].getSize() << " " << AdjacencyList[l].getSize() << endl;
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
            for (int y = 0; y < ListC[i].size(); y++){
                CycleNode[ListC[i][y]].push_back(ListCycles.size() - 1);
            }
        }
        //if (ListCycles.size() == ThisMany) terminate = true;
        
        //Remove node
        for (int l = 0; l < AdjaList.getSize(); l++){
            if (l == origin){
                //cout << AdjaList[l].getSize() << endl;
                AdjaList[l] = IloNumArray(env,0);
                //cout << AdjaList[l].getSize() << endl;
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
    //cout << endl << dicX[sp][this->vL[sp] - 1] << " " << this->vL[sp];
    while (cont == true){
        //cout << AdjaList[whoi][stepinSU[whoi]] << endl;
        if (IsxinStack(AdjaList[whoi][stepinSU[whoi]] - 1, Stack) == true || Stack.size() >= CycleLength){
            if (AdjaList[whoi][stepinSU[whoi]] == (origin + 1) && VLpresent == true && NewCycle == false){
                // they form a cycle together!
                weight = PathWeight(Stack);
                ListC.push_back(Stack);
//                for (int y = 0; y < Stack.size(); y++){
//                    CycleNode[Stack[y]].push_back(ListCycles.size() - 1);
//                }
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
                            //return BigW;
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
        //cout << z[i].getName() << endl;
    }
    
    //Add constraints
//    for (auto it = NodeinCycle.begin(); it != NodeinCycle.end(); it++){
//        IloExpr cycle (env,0);
//        auto pos = NodeinCycle.find(it->first);
//        if (pos != NodeinCycle.end()){
//            for (int j = 0; j < NodeinCycle[it->first].size(); j++){
//                cycle+= z[NodeinCycle[it->first][j]];
//            }
//        }
//        //Add constraint: node only in one cycle
//        cf.add(cycle <= 1);
//        cycle.end();
//    }
    
    
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
        if (NodeinCycle.count(i) == true){
            for (int j = 0; j < NodeinCycle[i].size(); j++){
                //exp+= z[NodeinCycle[i][j]];
            }
        }
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
            //cout << AdjacencyList[PredList[i][j]].getSize() << endl;
            units-= x[PredList[i][j]][h];
        }
        cout << units << endl;
        flow.add(IloRange(env,-IloInfinity, units, 0));
        units.end();
    }
    
    //Objective value
    IloExpr Obj(env, 0);
    for (int j = 0; j < AllCyclesWeight.size(); j++){
        //Obj += AllCyclesWeight[j]*z[j];
    }
    
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
        //DualBin = cplex.getDual(bin);
    }
    x.end();
    assignment.end();
    //bin.end();
    
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
int Problem::DijktrasAlgorithm(map<int,vector<int>>& Comp2, int source, int target) {
    map<int,double>minDistFromSourceToNode;
    for (auto it = Comp2.begin(); it != Comp2.end(); it++){
        minDistFromSourceToNode[it->first] = 1000000;
    }
    minDistFromSourceToNode[source] = 0; // l = source
    set <pair<int, int>> active_vertices;
    active_vertices.insert({ 0, source });
    

    // Shortest Path from source node to every other node i
    while (!active_vertices.empty()) {
        int where = active_vertices.begin()->second;
        if (where == target) break; // Cuando se ha llegado al nodo objetivo
        active_vertices.erase(active_vertices.begin());
        for (int i = 0; i < Comp2[where].size(); i++) {
            int node = Comp2[where][i];
            if (minDistFromSourceToNode[node] > minDistFromSourceToNode[where] + 1) { // El 1 es porque es un arco m·s wij = 1 para este Shortest Path
                active_vertices.erase({minDistFromSourceToNode[node], node});
                minDistFromSourceToNode[node] = minDistFromSourceToNode[where] + 1;
                active_vertices.insert({minDistFromSourceToNode[node], node});
            }
        }
    }
    return minDistFromSourceToNode[target];
}
vector<map<int,int>> Problem::FindEquals(map<int,vector<int>>& Comp2, int oriMDD){
    vector<int> equals;
    vector<map<int,int>>GroupofEquals;//node, group
    int minDist = 0, Threshold = 0;
    int source1, source2, target1, target2, count = 0, nElements;
    
    if (oriMDD > Pairs - 1){//chain
        Threshold = int(ChainLength);
    }
    else{
        Threshold = int(CycleLength);
    }
    
    for (auto it = ListOfEquals.back().begin(); it != ListOfEquals.back().end(); it++){
        GroupofEquals.push_back(map<int,int>());
        for (int g = 0; g < it->second.size(); g++){
            for (int s = 0; s < it->second[g].size(); s++){
                if (GroupofEquals.back().count(it->second[g][s]) == false){
                    count++;
                    equals.clear();
                    equals.push_back(it->second[g][s]);
                    GroupofEquals.back()[it->second[g][s]] = count;
                    nElements = int(GroupofEquals.back().size());
                    for (int j = s + 1; j < it->second[g].size(); j++){
                        if (GroupofEquals.back().count(it->second[g][j]) == false){
                            bool isEqual = true;
                            for (int h = 0; h < equals.size(); h++){
                                source1 = equals[h];
                                target1 = it->second[g][j];
                                minDist = DijktrasAlgorithm(Comp2, source1, target1);
                                minDist += it->first; // + because of the position from the MDD source node to its succesor (source)
                                if (minDist == Threshold) minDist++;
                                if (oriMDD < Pairs){
                                    if (minDist < Threshold){
                                        source2 = it->second[g][j];
                                        target2 = oriMDD; // MDD source node
                                        minDist+= DijktrasAlgorithm(Comp2, source2, target2);
                                    }
                                }
                                if (minDist > Threshold){//Check min distance starting with target1
                                    source1 = target1;
                                    target1 = equals[h];
                                    minDist = DijktrasAlgorithm(Comp2, source1, target1);
                                    minDist += it->first; // +1 because of the arc from the MDD source node to its succesor (source)
                                    if (minDist == Threshold) minDist++;
                                    if (oriMDD < Pairs){
                                        if (minDist < Threshold){
                                            source2 = it->second[g][j];
                                            target2 = oriMDD; // MDD source node
                                            minDist+= DijktrasAlgorithm(Comp2, source2, target2);
                                        }
                                    }
                                    if (minDist <= Threshold){//Not Equal
                                        isEqual = false;
                                        break;
                                    }
                                }
                                else{//Not equal
                                    isEqual = false;
                                    break;
                                }
                            }
                            if (isEqual  == true){
                                equals.push_back(it->second[g][j]);
                                GroupofEquals.back()[it->second[g][j]] = count; //Assign to a team
                            }
                        }
                    }
                    if (GroupofEquals.back().size() == nElements){// No one new in the current team
                        auto pos = GroupofEquals.back().find(it->second[g][s]); //it->second[s] still available to be matched
                        GroupofEquals.back().erase(pos);
                    }
                }
            }
        }
    }
    return GroupofEquals;
}
double Problem::Heuristic(){
    vector<CopySelection> sCopies(ListCycles.size());
    float sum = 0;
    int ThisMany = 0, times = 0;
    double n_rand, BestCyObj = 0;
    map<int,bool>sel_cycles, sel_cyclesB; //copy,cycle
    map<int,bool>sel_node;
    vector<int> estos;
    
    for (int i = 0; i < ListCycles.size(); i++){
        sort(ListCycles[i].second.begin(), ListCycles[i].second.end(), sortCycles);
        for (int j = 0; j < ListCycles[i].second.size(); j++){
//            if (ListCycles[i].second[j].fraction == 0){
//                sCopies[i].sel_prob += 0.15;
//            }
//            else{
                sCopies[i].sel_prob += ListCycles[i].second[j].fraction;
//            }
        }
        sum+= sCopies[i].sel_prob;
        if (sCopies[i].sel_prob > 0.5) ThisMany++;
    }
    //Sort List of Cycles
    sort(sCopies.begin(), sCopies.end(), sortsCopies);
    sCopies[0].sel_prob = sCopies[0].sel_prob/sum;
    
    for (int i = 1; i < sCopies.size(); i++){
        sCopies[i].sel_prob = (sCopies[i].sel_prob)/sum + sCopies[i - 1].sel_prob;
    }
    
    //cout << n_rand << endl;
    
    ThisMany = ThisMany - 6;
    int Total = max(0,ThisMany);
    uniform_real_distribution<double> dUniform(0.0,1.0);
    
    while(Total < sCopies.size() - 1){
        Total++;
        times = 0;
        while (times < 500){
            times++;
            int sel = 0;
            vector<int>OrderList;
            while(sel < Total){
                n_rand = dUniform(generator);
                if (n_rand > 0.5){
                    for (int i = sCopies.size() - 1; i >= 0; i--){
                        if (sCopies[i].sel_prob < n_rand && sCopies[i].selected == false){
                            sCopies[i].selected = true;
                            OrderList.push_back(i);
                            sel++;
                            break;
                        }
                    }
                }else{
                    for (int i = 0; i < sCopies.size(); i++){
                        if (sCopies[i].sel_prob < n_rand && sCopies[i].selected == false){
                            sCopies[i].selected = true;
                            OrderList.push_back(i);
                            sel++;
                            break;
                        }
                    }
                }
            }
            sel_node.clear();
            estos.clear();
            double obj = 0;
            sel_cycles.clear();
            
            //Pick Cycles
            for (int h = 0; h < OrderList.size(); h++){
                if (sCopies[h].selected == true){
                    uniform_int_distribution<int> dUniformInt(0,ListCycles[OrderList[h]].second.size() - 1);
                    vector<int>order;
                    while(order.size() < ListCycles[OrderList[h]].second.size()){
                        int nn = dUniformInt(generator);
                        order.push_back(nn);
                    }
                    for (int i = 0; i < order.size(); i++){
                        bool isSafe = true;
                        estos.clear();
                        for (int l = 0; l < ListCycles[OrderList[h]].second[order[i]].acycle.size(); l++){
                            if (sel_node.count(ListCycles[OrderList[h]].second[order[i]].acycle[l]) == false){
                                estos.push_back(ListCycles[OrderList[h]].second[order[i]].acycle[l]);
                            }
                            else{
                                isSafe = false;
                                break;
                            }
                        }
                        if (isSafe == true){
                            obj+= ListCycles[OrderList[h]].second[order[i]].TotalWeight;
                            sel_cycles[ListCycles[OrderList[h]].second[order[i]].indexAllC] = true;
                            for (int l = 0; l < estos.size(); l++){
                                sel_node[estos[l]] = true;
                            }
                        }
                    }
                }
            }
            if (BestCyObj < obj){
                BestCyObj = obj;
                sel_cyclesB = sel_cycles;
                //cout << "BestObj: " << BestCyObj << endl;
            }
            for (int l = 0; l < sCopies.size(); l++){
                sCopies[l].selected = false;
            }
        }
    }
    
    //IncumbentVars = IloNumVarArray(env,0);
    IncumbentVals = IloNumArray(env,0);

    /////////// x.l.i ///////////////
    for (int l = 0; l < AllCycles.size(); l++) {
        //IncumbentVars.add(z[l]);
        if (sel_cyclesB.count(l) == true){
            IncumbentVals.add(1);
        }
        else{
            IncumbentVals.add(0);
        }
    }
    
    return BestCyObj;
}

