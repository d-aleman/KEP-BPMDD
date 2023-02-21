//
//  Chains.cpp
//  LagrangeBranch&Bound
//
//  Created by Carolina Riascos Alvarez on 2020-06-11.
//  Copyright © 2020 Carolina Riascos Alvarez. All rights reserved.
//

#include "Chains.hpp"



vector<vector<int>> Problem::Chains(){
    
    vector<vector<int>>Bestchains;
    if (FirstTimeChains == true){
        FirstTimeChains = false;
        //Create model
        
        ch = IloModel(env);
        
        //Create variables
        x = IloArray<IloNumVarArray>(env, AdjacencyList.getSize());
        for (int i = 0; i < AdjacencyList.getSize(); i++){
           x[i] = IloNumVarArray(this->env, AdjacencyList[i].getSize(), 0, 1, ILOFLOAT);
           for (int j = 0; j < AdjacencyList[i].getSize(); j++){
              SetName2(x[i][j], "x", i + 1, AdjacencyList[i][j]);
               //cout << x[i][j].getName() << endl;
           }
        }
        
        cons = IloRangeArray(env);
        for (int i = 0; i < Nodes + 2; i++){
            IloExpr units(env,0);
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                units+= x[i][j];
            }
            ch.add(units <= 1);
            for (int j = 0; j < PredList[i].size(); j++){
                int h = 0;
                for (h = 0; h < AdjacencyList[PredList[i][j]].getSize(); h++){
                    if (AdjacencyList[PredList[i][j]][h] == i + 1) {
                       break;
                    }
                }
                units-= x[PredList[i][j]][h];
                //cout << units << endl;
            }
            //cout << units << endl;
            if (i < Nodes){
                cons.add(IloRange(env,0, units, 0));
            }
            else if (i == Nodes){
                cons.add(IloRange(env, 1, units , 1));
                //cout << units << endl;
            }
            else{
                cons.add(IloRange(env, -1, units , -1));
                //cout << units << endl;
            }
            units.end();
        }
        ch.add(cons);
        
        //Objective value
        ObjCH = IloObjective(env);
        IloExpr suma(env, 0);
        for (int i = 0; i < AdjacencyList.getSize() - 2; i++){
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                suma += (Weights[make_pair(i,AdjacencyList[i][j] - 1)] - solpi[i])*x[i][j];
            }
        }
        //cout << Obj << endl;
        //Objective
        ObjCH = IloAdd(ch, IloMaximize(env,suma));
        suma.end();
        
    }
    
    if (FirstTimeChains == false){
        IloNumArray vals(env);
        for (int i = 0; i < AdjacencyList.getSize() - 2; i++){
            vals = IloNumArray(env);
            if (AdjacencyList[i].getSize() > 0){
                for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                    vals.add(Weights[make_pair(i,AdjacencyList[i][j] - 1)] - solpi[i]);
                    //cout << i << " " << j << endl;
                }
                ObjCH.setLinearCoefs(x[i], vals);
            }
        }
        ch.add(ObjCH);
        vals.end();
    }
    
    if (insideBB == true){
        //Reset bounds
        for (int i = 0; i < AdjacencyList.getSize(); i++){
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                x[i][j].setUB(1);
            }
        }
        
        //Apply branching
        int nodei = -1;
        int nodej = -1;
        
        
        if (BranchingMethod == 2){
            for (int i = 0; i < Tree[VoyBB.first][VoyBB.second].ArcGroupZero.size(); i++){
                nodei = Tree[VoyBB.first][VoyBB.second].ArcGroupZero[i].first;
                nodej = Tree[VoyBB.first][VoyBB.second].ArcGroupZero[i].second;
                int j = 0;
                for (j; j < AdjacencyList[nodei].getSize(); j++){
                    if (AdjacencyList[nodei][j] - 1 == nodej){
                        break;
                    }
                }
                //cout << "Zero: " << x[nodei][j].getName() << endl;
                x[nodei][j].setUB(0);
            }
            for (int i = 0; i < Tree[VoyBB.first][VoyBB.second].ArcGroupOne.size(); i++){
                nodei = Tree[VoyBB.first][VoyBB.second].ArcGroupOne[i].first;
                nodej = Tree[VoyBB.first][VoyBB.second].ArcGroupOne[i].second;
                for (int j = 0; j < AdjacencyList[nodei].getSize(); j++){
                    if (AdjacencyList[nodei][j]  - 1 != nodej){
                        x[nodei][j].setUB(0);
                    }
                    else{
                      //x[nodei][j].setLB(1);
                      //cout << "One: " <<  x[nodei][j].getName() << endl;
                    }
                }
            }
        }
    }
    //cout << "Start Phase 2" << endl;
    //Solve LP
    IloCplex cex(ch);
    cex.setOut(env.getNullStream());
    cex.setParam(IloCplex::Param::TimeLimit, 1800);
    cex.setParam(IloCplex::Param::Threads, 1);
    cex.solve();
    if (cex.getStatus() == IloAlgorithm::Infeasible) {
        env.out() << "No chains found" << endl;
    }
    else {
        double Obj_B = cex.getObjValue();
        //cout << "Obj chains" << Obj_B << endl;
        itPH2LP++;
        
        x_sol = IloNumArray2(env, AdjacencyList.getSize());
        for (int i = 0; i < AdjacencyList.getSize(); i++){
            x_sol[i] = IloNumArray(env, 0);
            cex.getValues(x[i], x_sol[i]);
        }
        
        if (Obj_B > 0.0001){

            //Check solution
            InitializeVertexinSolChain(x_sol, VertexinSolChain);
            vector<vector<int>>Bestchains;
            Bestchains = FindChains(x_sol);
            
            if (Bestchains.size() == 0){
                //Check whether there are feasible cycles
                Bestchains = FindCyclesFromIntSol(x_sol);
                //Check whether cyccles are positive priced
                vector<vector<int>>::iterator iv;
                for (iv = Bestchains.begin(); iv != Bestchains.end();)
                    if (PtivePriceCycle(solpi, Weights, *iv) == false){
                        Bestchains.erase(iv);
                    }else{
                        iv++;
                    }
                
                if (Bestchains.size() > 0){
                    //Return BestChains
                    nColCyclePH2+=Bestchains.size();
                    return Bestchains; //COUNT nCylePh2
                }
                else{
                    //Impose integrality and arc constraints on chains and resolve
                    double RestrictedObj;
                    IloExpr sum(env,0);
                    
                    for (int j = 0; j < AdjacencyList.getSize() - 1; j++){
                        sum+= IloSum(x[j]);
                    }
                    IloRange arcs(IloRange(env, -IloInfinity, sum, ChainLength + 2));
                    ch.add(arcs);
                    
                    if (FirstConversion == true){
                        FirstConversion = false;
                        for (int j = 0; j < AdjacencyList.getSize() - 1; j++){
                            ch.add(IloConversion(env, x[j], ILOBOOL));
                        }
                    }
                    
                    //Resolve model
                    cex.solve();
                    
                    if (cex.getStatus() == IloAlgorithm::Infeasible) {
                        env.out() << "No solution" << endl;
                    }
                    else {
                        RestrictedObj = cex.getObjValue();
                        itPH2MIP++;
                        //cout <<  RestrictedObj << endl;
                        if (RestrictedObj < 0.0001){
                            //Return empty array. No positive price chains found
                            //Remove additional constraints
                            ch.remove(arcs);
                            //ch.remove(flow);
                            return Bestchains;
                        }
                        else{
                            //Find the positive price chain
                            for (int i = 0; i < AdjacencyList.getSize() - 2; i++){
                                x_sol[i] = IloNumArray(env, 0);
                                cex.getValues(x[i], x_sol[i]);
                            }
                            InitializeVertexinSolChain(x_sol, VertexinSolChain);
                            Bestchains = FindChains(x_sol);
                            
                            if (Bestchains.size() > 0){//If we can find a positive price chain
                                //Return chains
                                //Remove additional constraints
                                ch.remove(arcs);
                                //ch.remove(flow);
                                nColChainPH2MIP+= Bestchains.size();
                                //cout << " End Phase 2" << endl;
                                return Bestchains; //nChainPh2
                            }
                            else{
                                //There must exist subtours. Check solution again
                                Bestchains = FindCyclesFromIntSol(x_sol);
                                //Check whether cyccles are positive priced
                                vector<vector<int>>::iterator iv;
                                for (iv = Bestchains.begin(); iv != Bestchains.end();)
                                    if (PtivePriceCycle(solpi, Weights, *iv) == false){
                                        Bestchains.erase(iv);
                                    }else{
                                        iv++;
                                    }
                                
                                if (Bestchains.size() > 0){// If we find a feasible positive priced cycle
                                    //Remove additional constraints
                                    ch.remove(arcs);
                                    
                                    //Return cycle
                                    return Bestchains;//nCyclePh2
                                }
                                else{
                                    //Introduce a cut removing the infeasible cycles and resolve
                                    // The cutting plane should return either a positive chain or cycle, or a certificate that no one of such exists
                                    itPH2Sep++;
                                    vector<vector<int>>subtours;
                                    int rhs = 0;
                                    while(true){
                                        //Find subtout
                                        subtours = FindSubtour(x_sol);
                                        if (subtours.size() == 0) break;
                                        //cout << "Size subtour: " << subtour.size() << endl;
                                        for (int h = 0; h < subtours.size(); h++){
                                            IloExpr expr(env,0);
                                            for (int j = 0; j < subtours[h].size() - 1; j++){
                                                int s = 0;
                                                for (s; s < AdjacencyList[subtours[h][j]].getSize(); s++){
                                                    if (AdjacencyList[subtours[h][j]][s] - 1 == subtours[h][j + 1]){
                                                        break;
                                                    }
                                                }
                                                expr+= x[subtours[h][j]][s];
                                            }
                                            rhs = int(subtours[h].size() - 2);
                                            //cout << expr << endl;
                                            ch.add(expr <= rhs);
                                            //expr.end();
                                        }
                                        cex.solve();
                                        if (cex.getStatus() == IloAlgorithm::Infeasible) {
                                            env.out() << "No solution" << endl;
                                        }
                                        else {
                                            RestrictedObj = cex.getObjValue();
                                            if (RestrictedObj > 0.0001){
                                                for (int i = 0; i < AdjacencyList.getSize() - 2; i++){
                                                    x_sol[i] = IloNumArray(env, 0);
                                                    cex.getValues(x[i], x_sol[i]);
                                                }
                                                //Check for a positive priced chain
                                                InitializeVertexinSolChain(x_sol, VertexinSolChain);
                                                Bestchains = FindChains(x_sol);
                                                if (Bestchains.size() > 0){
                                                    //Remove additional constraints
                                                    ch.remove(arcs);
                                                    //ch.remove(flow);
                                                    nColChainPH2Sep+= Bestchains.size();
                                                    return Bestchains;//nChainPh2
                                                }
                                            }
                                            else{
                                                //Remove additional constraints
                                                ch.remove(arcs);
                                                //ch.remove(flow);
                                                return Bestchains;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }else{
                cex.end();
                nColChainPH2LP+= Bestchains.size();
                return Bestchains;//COUNT nChainPh2
            }
        }else{
            cex.end();
            //cout << "Zero";
        }
    }

    return Bestchains;

}
void Problem::InitializeVertexinSolChain(IloNumArray2& x_sol, vector<vChain>& VertexinSolChain){
    for (int i = AdjacencyList.getSize() - 3; i >= 0 ; i--){//2 dummy nodes not counted
        if (VertexinSolChain[Nodes - 1 - i].e_sol.size() > 0) VertexinSolChain[Nodes - 1 - i].e_sol = vector<int>();
        for (int j = 0; j < x_sol[i].getSize(); j++){
            if (x_sol[i][j] > 0.7 && AdjacencyList[i][j] - 1 < Nodes){
                VertexinSolChain[Nodes - 1 - i].e_sol.push_back(AdjacencyList[i][j] - 1);
            }
        }
        VertexinSolChain[Nodes - 1 - i].it = VertexinSolChain[Nodes - 1 - i].e_sol.begin();
        VertexinSolChain[Nodes - 1 - i].itEnd = VertexinSolChain[Nodes - 1 - i].e_sol.end();
    }
}
vector<vector<int>> Problem::FindChains(IloNumArray2 x_sol){
    vector<vector<int>>Bestchains;
    int u,v;
    double w = 0, price = 0;
    vector<vChain>::iterator itVinChain = VertexinSolChain.begin();
    vector<vChain>::iterator itVinChainAux;
    vector<Chain>PPChains;
    
    while (itVinChain->vertex > Pairs - 1){//By construction altruistic donors come first
        if (itVinChain->e_sol.size() > 0){
            PPChains.push_back(Chain(*itVinChain));
            PPChains.back().AccumWeight = -1*solpi[itVinChain->vertex];
            while(PPChains.back().Vnodes.size() >  0){//next neighbor

                u = PPChains.back().Vnodes.back().vertex;
                v = *(PPChains.back().Vnodes.back().it);
                w = Weights[make_pair(u, v)];
                price = PPChains.back().AccumWeight + w - solpi[v];
                //Find pointer to v in VertexinSolChain
                itVinChainAux = VertexinSolChain.begin();
                for (itVinChainAux; itVinChainAux != VertexinSolChain.end(); itVinChainAux++)
                    if (itVinChainAux->vertex == v){
                        break;
                    }
                if (PPChains.back().AccumWeight < price){
                    //If vertex not already in chain
                    if (v2AlreadyinChain(PPChains.back().Vnodes, v) == false){
                        //Add vertex and update AccumWeight
                        PPChains.back().Vnodes.push_back(*itVinChainAux);
//                        for (auto v: PPChains.back().Vnodes)
//                            v.it = v.e_sol.begin();
                        PPChains.back().AccumWeight = price;
                    }
                    //Else, do nothing. It will visit the next neighbor
                }
                else{
                    //If vertex not already in chain
                    if (v2AlreadyinChain(PPChains.back().Vnodes, v) == false){
                        //If PPChains.back().AccumWeight > 0, then store current chain and duplicate it to start a new one
                        if (PPChains.back().AccumWeight > 0.0001 && PPChains.back().Vnodes.size() >= 2){
                            Chain aux = PPChains.back();
                            PPChains.push_back(aux);
                            //Add vertex to new chain
                            PPChains.back().Vnodes.push_back(*itVinChainAux);
                            PPChains.back().AccumWeight = price;
                        }
                        else{//Else, add vertex to current chain
                            PPChains.back().Vnodes.push_back(*itVinChainAux);
                            PPChains.back().AccumWeight = price;
                        }
                    }
                    //Else, do nothing. It will visit the next neighbor
                };
                //Increase iterator
                if (PPChains.back().Vnodes.back().FirstPass == true){
                    PPChains.back().Vnodes.back().FirstPass = false;
                }else{PPChains.back().Vnodes.back().it++;}
                //If chain is about to become too long or we have reached the end of some vertex's neighbors
                if (PPChains.back().Vnodes.size() > ChainLength || PPChains.back().Vnodes.back().it == PPChains.back().Vnodes.back().e_sol.end()){
                    //If current chain  positive, store, duplicate and pop_back new one
                    if (PPChains.back().AccumWeight > 0.0001){
                        Chain aux = PPChains.back();
                        PPChains.push_back(aux);
                        PPChains.back().Vnodes.pop_back();
                        //Find new neighbor
                        FindNewNeighbor(PPChains);
                        
                    }else{//If current chain negative, pop_back
                        PPChains.back().Vnodes.pop_back();
                        //Find new neighbor
                        FindNewNeighbor(PPChains);
                    }
                }
            }
            PPChains.erase(PPChains.end() - 1);
        }
        itVinChain++;
    }
    //Transfer chains to Bestchains
    for (int i = 0; i < PPChains.size();i++){
        Bestchains.push_back(vector<int>());
        for (int j = 0; j < PPChains[i].Vnodes.size(); j++){
            Bestchains.back().push_back(PPChains[i].Vnodes[j].vertex);
        }
    }
    
    if (Bestchains.size() > 0){
        if (Bestchains[0].size() == 1){
            cout << "S.O.S" << endl;
        }
    }
    return Bestchains;
    
}
void FindNewNeighbor(vector<Chain>& PPChains){
    while (PPChains.back().Vnodes.back().it != PPChains.back().Vnodes.back().e_sol.end()){
        PPChains.back().Vnodes.back().it++;
        if (PPChains.back().Vnodes.back().it == PPChains.back().Vnodes.back().itEnd){
            PPChains.back().Vnodes.pop_back();
        }
        else{
            break;
        }
        if (PPChains.back().Vnodes.size() == 0) break;
    }
}
bool v2AlreadyinChain(vector<vChain> v1, int v2){
    for (int j = 0; j < v1.size(); j++){
        if (v1[j].vertex == v2){
            return true;
        }
    }
    return false;
}
bool isInside(vector<int>chain, int node){
    for (int j = 0; j < chain.size(); j++){
        if (node == chain[j]){
            return true;
        }
    }
    return false;
}
bool isInside2(vector<int>chain, int node){
    for (int j = 1; j < chain.size(); j++){
        if (node == chain[j]){
            return true;
        }
    }
    return false;
}
vector<vector<int>>Problem::LeftCycles(){
    FastExit = false;
    vector<vector<int>>Bestchains;
    if (FirstTimeCycles == true){
        FirstTimeCycles = false;
        x_sol = IloNumArray2(env, AdjacencyList.getSize());
        for (int i = 0; i < AdjacencyList.getSize(); i++){
            x_sol[i] = IloNumArray(env, AdjacencyList[i].getSize());
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                if (AdjacencyList[i][j] - 1 < Pairs){
                    x_sol[i][j] = 1;
                }
                else{
                    x_sol[i][j] = 0;
                }
            }
        }
    }

    if (insideBB == true){
        //Reset bounds
        for (int i = 0; i < Pairs; i++){
            for (int j = 0; j < AdjacencyList[i].getSize(); j++){
                x_sol[i][j] = 1;
            }
        }
        //Apply branching
        int nodei = -1;
        int nodej = -1;
        if (BranchingMethod == 2){
            for (int i = 0; i < Tree[VoyBB.first][VoyBB.second].ArcGroupZero.size(); i++){
                nodei = Tree[VoyBB.first][VoyBB.second].ArcGroupZero[i].first;
                nodej = Tree[VoyBB.first][VoyBB.second].ArcGroupZero[i].second;
                if (nodei < Pairs){
                    int j = 0;
                    for (j; j < AdjacencyList[nodei].getSize(); j++){
                        if (AdjacencyList[nodei][j] - 1 == nodej){
                            break;
                        }
                    }
                    //cout << "Zero: " << x[nodei][j].getName() << endl;
                    x_sol[nodei][j] = 0;
                }
            }
            for (int i = 0; i < Tree[VoyBB.first][VoyBB.second].ArcGroupOne.size(); i++){
                nodei = Tree[VoyBB.first][VoyBB.second].ArcGroupOne[i].first;
                nodej = Tree[VoyBB.first][VoyBB.second].ArcGroupOne[i].second;
                if (nodei < Pairs){
                    for (int j = 0; j < AdjacencyList[nodei].getSize() - 1; j++){
                        if (AdjacencyList[nodei][j]  - 1 != nodej){
                            x_sol[nodei][j] = 0;
                        }
                        else{
                        //cout << "One: " <<  x[nodei][j].getName() << endl;
                        }
                    }
                }
            }
        }
    }
    Bestchains = FindFeasibleCycles(x_sol);
    itPH3Search++;
    
    if (Bestchains.size() > 0){
        nColCyclePH3Sch+= Bestchains.size();
        return Bestchains;
    }
    else{
        
        //cout << " Start Phase 3" << endl;
        if (FirstTimeMIPCY == true){
             //Create model
             cy = IloModel(env);

             //Create variables
             y = IloArray<IloNumVarArray>(env, Pairs);
             for (int i = 0; i < Pairs; i++){
                y[i] = IloNumVarArray(this->env, AdjacencyList[i].getSize() - 1, 0, 1, ILOINT);
                for (int j = 0; j < AdjacencyList[i].getSize() - 1; j++){
                   SetName2(y[i][j], "y", i + 1, AdjacencyList[i][j]);
                    //cout << y[i][j].getName() << endl;
                }
             }

             //Constraints
             IloRangeArray cons(env);
             for (int i = 0; i < Pairs; i++){
                 IloExpr units(env,0);
                 for (int j = 0; j < AdjacencyList[i].getSize() - 1; j++){
                     units+= y[i][j];
                 }
                 cy.add(units <= 1);
                 for (int j = 0; j < PredList[i].size(); j++){
                     int h = 0;
                     if (PredList[i][j] < Pairs){
                         for (h = 0; h < AdjacencyList[PredList[i][j]].getSize(); h++){
                             if (AdjacencyList[PredList[i][j]][h] == i + 1) {
                                break;
                             }
                         }
                         units-= y[PredList[i][j]][h];
                     }
                 }
                 cons.add(IloRange(env,0, units, 0));
                 units.end();
             }
             cy.add(cons);
             
            //Limit number of arcs to K
             IloExpr sum(env,0);
             for (int j = 0; j < Pairs; j++){
                 sum+= IloSum(y[j]);
             }
             IloRange arcs(IloRange(env, -IloInfinity, sum, CycleLength));
             cy.add(arcs);
             sum.end();

             //Objective value
             ObjCY = IloObjective(env);
             IloExpr suma(env, 0);
             for (int i = 0; i < Pairs; i++){
                 for (int j = 0; j < AdjacencyList[i].getSize() - 1; j++){
                     suma += (Weights[make_pair(i,AdjacencyList[i][j] - 1)] - solpi[i])*y[i][j];
                 }
             }

             //Objective
             ObjCY = IloAdd(cy, IloMaximize(env,suma));
             suma.end();
        }
        if (FirstTimeMIPCY == false){
            IloNumArray vals(env);
            for (int i = 0; i < Pairs; i++){
                vals = IloNumArray(env);
                for (int j = 0; j < AdjacencyList[i].getSize() - 1; j++){
                    vals.add(Weights[make_pair(i,AdjacencyList[i][j] - 1)] - solpi[i]);
                    //cout << i << " " << j << endl;
                }
                if (AdjacencyList[i].getSize() - 1 > 0){
                    ObjCY.setLinearCoefs(y[i], vals);
                }
            }
            cy.add(ObjCY);
            vals.end();
        }
        FirstTimeMIPCY = false;
        
        //IloNumArray2 y_sol(env, Pairs);
        if (insideBB == true){
            //Reset bounds
            for (int i = 0; i < Pairs; i++){
                for (int j = 0; j < AdjacencyList[i].getSize() - 1; j++){
                    y[i][j].setUB(1);
                }
            }
            //Apply branching
            int nodei = -1;
            int nodej = -1;
            if (BranchingMethod == 2){
                for (int i = 0; i < Tree[VoyBB.first][VoyBB.second].ArcGroupZero.size(); i++){
                    nodei = Tree[VoyBB.first][VoyBB.second].ArcGroupZero[i].first;
                    nodej = Tree[VoyBB.first][VoyBB.second].ArcGroupZero[i].second;
                    if (nodei < Pairs){
                        int j = 0;
                        for (j; j < AdjacencyList[nodei].getSize(); j++){
                            if (AdjacencyList[nodei][j] - 1 == nodej){
                                break;
                            }
                        }
                        y[nodei][j].setUB(0);
                    }
                }
                for (int i = 0; i < Tree[VoyBB.first][VoyBB.second].ArcGroupOne.size(); i++){
                    nodei = Tree[VoyBB.first][VoyBB.second].ArcGroupOne[i].first;
                    nodej = Tree[VoyBB.first][VoyBB.second].ArcGroupOne[i].second;
                    if (nodei < Pairs){
                        for (int j = 0; j < AdjacencyList[nodei].getSize() - 1; j++){
                            if (AdjacencyList[nodei][j]  - 1 != nodej){
                                y[nodei][j].setUB(0);
                            }
                            else{
                              //cout << "One: " <<  x[nodei][j].getName() << endl;
                            }
                        }
                    }
                }
            }
        }
    
        //Solve LP
        IloCplex cex(cy);
        cex.setOut(env.getNullStream());
        cex.setParam(IloCplex::Param::TimeLimit, 1800);
        cex.setParam(IloCplex::Param::Threads, 1);
        cex.solve();
    
        //Solve
        cex.solve();
        itPH3MIP++;
        double Obj_B = 0;
        if (cex.getStatus() == IloAlgorithm::Infeasible) {
            env.out() << "It should NOT be here" << endl;
        }else{
            Obj_B = cex.getObjValue();
            //cout << "Phase 3" << endl;
            if (Obj_B > 0.0001){
                //Get solution values
                IloNumArray2 y_sol(env, Pairs);
                for (int i = 0; i < Pairs; i++){
                    y_sol[i] = IloNumArray(env, 0);
                    cex.getValues(y[i], y_sol[i]);
                }
                //Find a positive price cycle
                Bestchains = FindCyclesFromIntSol(y_sol);
                if (Bestchains.size() == 0){
                    cout << "This should NOT happen";
                }
                cex.end();
                nColCyclePH3MIP+= Bestchains.size();
                return Bestchains;
            }
            else{
                return Bestchains;
            }
        }
    }
    return Bestchains;
}
vector<vector<int>> Problem::FindCyclesFromIntSol(IloNumArray2& p_sol){
    vector<vector<int>>cycles;
    
    //Update values
    for (int i = 0; i < ArcsinSol.size(); i++){
        ArcsinSol[i].val = p_sol[ArcsinSol[i].row][ArcsinSol[i].column]; //It will get updated as the algorithm runs
    }
    
    //Find feasible cycles
    bool backToStart = false;
    int posCol = 0;
    for (int i = 0; i < ArcsinSol.size(); i++){
        if (ArcsinSol[i].val > 0.7 && ArcsinSol[i].selected == false){
            //New cycle
            cycles.push_back(vector<int>());
            cycles.back().push_back(ArcsinSol[i].row);
            backToStart = false;
            int NewRow = AdjacencyList[ArcsinSol[i].row][ArcsinSol[i].column] - 1;
            posCol = ArcsinSol[i].column;
            while (backToStart == false){
                //Add new vertex
                cycles.back().push_back(NewRow);
                //Mark arc as selected
                MarkArc(ArcsinSol,cycles.back(), posCol, i);
                //Find new row
                NewRow = FindNewRow(posCol,NewRow, p_sol, AdjacencyList);
                //Drop no feasible cycle
                if (cycles.back().size() > CycleLength || NewRow == -1){
                    cycles.erase(cycles.begin() + cycles.size() - 1);
                    backToStart = true;
                }
                if (cycles.size() > 0){
                    if (NewRow == cycles.back()[0]){
                        cycles.back().push_back(NewRow);
                        MarkArc(ArcsinSol,cycles.back(), posCol, i);
                        backToStart = true;
                    }
                }
            }
        }
    }
    
    //Restore original values
    for (int i = 0; i < ArcsinSol.size(); i++){
        ArcsinSol[i].val = 0;
        ArcsinSol[i].selected = false;
    }
    
    return cycles;
}
int FindNewRow(int& posCol, int NewCol, IloNumArray2& p_sol, IloNumArray2& AdjacencyList){
    for (int l = 0; l < p_sol[NewCol].getSize(); l++){
        if (p_sol[NewCol][l] > 0.9){
            posCol = l;
            return AdjacencyList[NewCol][l] - 1;
        }
    }
    return -1;
}
void MarkArc(vector<CheckedArc>&ArcsinSol, vector<int>&cycle, int&posCol, int& start){
    for (int l = start; l < ArcsinSol.size(); l++){
        if (ArcsinSol[l].row == cycle[cycle.size() - 2] && ArcsinSol[l].column == posCol){
            ArcsinSol[l].selected = true;
            break;
        }
    }
}
vector<vector<int>> Problem::FindFeasibleCycles(IloNumArray2& p_sol){
    
    //cout << mOrden.size() << endl;
    clock_t tStart = clock();
    double ThisTime = 0;
    vector<int>chain;
    double Totalchains = 0;
    vector<vector<int>>Bestchains;
    
    int counter = 0;
    for (auto it = Order.begin(); it != Order.end(); it++){
        counter++;
        //cout << counter << endl;
        int pos = it->first;
        if (vMDD[pos].first[1].actual < Pairs){
            //cout << i << endl;
            int donor = vMDD[pos].first[1].actual;
            double sum = 0;
            double maxsum = 0;
            chain.clear();
            vector<int>veci(AdjacencyList.getSize(), -1);
            veci[donor]++;
            //cout << veci[donor] << endl;
            for (int patient = 0; patient < AdjacencyList[donor].getSize() - 1; patient++){// -1 to not visit the dummy node
                if (chain.size() == 0){
                    chain.push_back(donor);
                }
                if (p_sol[donor][patient] >= 0.98 && chain.size() <= CycleLength + 1){
                    if (isInside2(chain, AdjacencyList[donor][veci[donor]] - 1) == true){
                        cout << "REPEATED NODES S.O.S";
                    }
                    chain.push_back(AdjacencyList[donor][patient] - 1);
                    sum=0;
                    for (int j = 0; j < chain.size() - 1; j++){
                        sum+= Weights[make_pair(chain[j], chain[j + 1])] - solpi[chain[j]];
                    }
                    if (sum > maxsum && sum > 0.005 && AdjacencyList[donor][patient] - 1 == vMDD[pos].first[1].actual && chain.size() <= CycleLength + 1){ //if it is a positive price cycle
                        maxsum = sum;
                        Bestchains.clear();
                        Bestchains.push_back(chain);
                        if (Bestchains.size() > 0) return Bestchains;
                    }
                    else if (sum == maxsum && sum > 0.005 && AdjacencyList[donor][patient] - 1 == vMDD[pos].first[1].actual && chain.size() <= CycleLength + 1){//if it is a positive price cycle
                        Bestchains.push_back(chain);
                    }
                    if (chain.size() == CycleLength + 1){
                        if (AdjacencyList[donor][patient] - 1 != vMDD[pos].first[1].actual){
                            veci[chain.back()] = -1;
                        }
                        chain.pop_back();
                        donor = chain.back();
                    }
                    else{
                        if (chain.back() == vMDD[pos].first[1].actual){
                            chain.pop_back();
                            veci[chain.back()] = -1;
                            chain.pop_back();
                            donor = chain.back();
                        }
                        else{
                            donor = AdjacencyList[donor][patient] - 1;
                        }
                    }
                }
                veci[donor]++;
                patient = veci[donor] - 1;
                if (veci[donor] < AdjacencyList[donor].getSize() - 1){
                    while (isInside2(chain, AdjacencyList[donor][veci[donor]] - 1) == true || isInOrderm(it, AdjacencyList[donor][veci[donor]] - 1) == true){
                        veci[donor]++;
                        patient = veci[donor] - 1;
                    }
                }
                ThisTime = (clock() - tStart)/double(CLOCKS_PER_SEC);
                if (ThisTime > 20){
                    break;
                }
                if (veci[donor] >= AdjacencyList[donor].getSize() - 1 && chain.size() > 1){
                    bool safe = false;
                    while (chain.size() >= 2 && safe == false){
                        chain.pop_back();
                        donor = chain.back();
                        while (true){
                            veci[donor]++;
                            patient = veci[donor] - 1;
                            if (veci[donor] < AdjacencyList[donor].getSize() - 1){
                                //cout << AdjacencyList[donor][veci[donor]] - 1 << endl;
                                if (isInside2(chain, AdjacencyList[donor][veci[donor]] - 1) == false && isInOrderm(it, AdjacencyList[donor][veci[donor]] - 1) == false){
                                    sum=0;
                                    for (int j = 0; j < chain.size() - 1; j++){
                                        sum+= Weights[make_pair(chain[j], chain[j + 1])] - solpi[chain[j]];
                                    }
                                    safe = true;
                                    break;
                                }
                            }
                            else{
                                break;
                            }
                        }
                    }
                }
            }
            if (maxsum > 0) Totalchains+=maxsum;
        }
        ThisTime = (clock() - tStart)/double(CLOCKS_PER_SEC);
        if (ThisTime > 10){
            FastExit = true;
            break;
        }
    }
    
    return Bestchains;
}
vector<vector<int>> Problem::FindSubtour(IloNumArray2 x_sol){
    
    vector<vector<int>>cycles;
    
    //Update values
    for (int i = 0; i < ArcsinSol.size(); i++){
            ArcsinSol[i].val = x_sol[ArcsinSol[i].row][ArcsinSol[i].column]; //It will get updated as the algorithm runs
    }
    
    //Find feasible cycles
    bool backToStart = false;
    int posCol = 0;
    for (int i = 0; i < ArcsinSol.size(); i++){
        if (ArcsinSol[i].val > 0.7 && ArcsinSol[i].selected == false){
            //New cycle
            cycles.push_back(vector<int>());
            cycles.back().push_back(ArcsinSol[i].row);
            backToStart = false;
            int NewRow = AdjacencyList[ArcsinSol[i].row][ArcsinSol[i].column] - 1;
            posCol = ArcsinSol[i].column;
            while (backToStart == false){
                //Add new vertex
                cycles.back().push_back(NewRow);
                //Mark arc as selected
                MarkArc(ArcsinSol,cycles.back(), posCol, i);
                //Find new row
                NewRow = FindNewRow(posCol,NewRow, x_sol, AdjacencyList);
                if (NewRow == cycles.back()[0]){
                    cycles.back().push_back(NewRow);
                    MarkArc(ArcsinSol,cycles.back(), posCol, i);
                    backToStart = true;
                    if (cycles.back().size() <= CycleLength + 1){
                        cycles.erase(cycles.begin() + cycles.size() - 1);
                    }
                }
                else if (NewRow == -1){//Drop no feasible cycle
                     cycles.erase(cycles.begin() + cycles.size() - 1);
                     backToStart = true;
                }
            }
        }
    }
    
    //Restore original values
    for (int i = 0; i < ArcsinSol.size(); i++){
        ArcsinSol[i].val = 0;
        ArcsinSol[i].selected = false;
    }
    
    return cycles;
}
bool PtivePriceCycle(IloNumArray& solpi, map<pair<int,int>,double>&Weights, vector<int>& cycle){
    double price = 0;
    
    for (int i = 0; i < cycle.size() - 1; i++)
        price+= Weights[make_pair(cycle[i], cycle[i + 1])] - solpi[cycle[i]];
    
    if (price > 0.0001){
        return true;
    }
    else{
        return false;
    }
}
void Problem::SetName2(IloNumVar& var, const char* prefix, IloInt i, IloInt j){
    string _prefix = prefix;
    string name = _prefix + "." + to_string(i) + "." + to_string(j);
    const char* varName = name.c_str();
    var.setName(varName);
}
bool Problem::isInOrderm(vector<pair<int,double>>::iterator& itVeci, int veci){
    for (auto it = Order.begin(); it != itVeci; it++){
        if (veci == vMDD[it->first].first[1].actual){
            return true;
        }
    }
    return false;
}
