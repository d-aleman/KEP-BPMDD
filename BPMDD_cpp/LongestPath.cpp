//
//  LongestPath.cpp
//
//  Created by Carolina Riascos Alvarez on 2019-12-16.
//  Copyright © 2019 Carolina Riascos Alvarez. All rights reserved.
//

#include "LongestPath.hpp"

tuple<int,int,double,int,double> Problem::maxMenosOcho(vector<state>& vLosNew, state& estoy, int i, bool cadena, int CuantosPahtsQueres){
    tuple<int,int,double,int,double> tCicloFinAux; double peso = 0; double mmax = -1000000.0;//-10000000.0;
    int k = estoy.vPrede[i];
    for (auto it : estoy.mActual[vLosNew[k]]){
        auto pp = make_pair(vLosNew[k].actual, it);
        peso = vLosNew[k].MaxWeight + Weights[pp] - solpi[vLosNew[k].actual]; //modified weight
        if (insideBB == true){
            if (BranchingMethod == 1){
                //Go over the node map ZERO
                auto pos = Tree[VoyBB.first][VoyBB.second].NodeGroupZero.find(posMDD);
                if (pos != Tree[VoyBB.first][VoyBB.second].NodeGroupZero.end()){
                    for (auto it2 = Tree[VoyBB.first][VoyBB.second].NodeGroupZero[posMDD].begin(); it2!= Tree[VoyBB.first][VoyBB.second].NodeGroupZero[posMDD].end(); it2++){
                        if (vLosNew[estoy.vPrede[i]].actual == *it2){
                            peso = -1000;
                            break;//found!
                        }
                    }
                }
            }else{
                //Go over the arc map ZERO
                for (auto it2 = Tree[VoyBB.first][VoyBB.second].ArcGroupZero.begin(); it2!= Tree[VoyBB.first][VoyBB.second].ArcGroupZero.end(); it2++){
                    if (vLosNew[k].actual == it2->first && it == it2->second){
                        peso = -1000;
                        break; //found!
                    }
                }
                //Go over the arc map ONE
                bool allowed = false;
                bool ifTail = false;
                for (auto it2 = Tree[VoyBB.first][VoyBB.second].ArcGroupOne.begin(); it2!= Tree[VoyBB.first][VoyBB.second].ArcGroupOne.end(); it2++){
                    if (vLosNew[k].actual == it2->first) ifTail = true;
                    if (vLosNew[k].actual == it2->first && it == it2->second){
                        allowed = true;
                        break;//found!
                    }
                }
                if (allowed == false && ifTail == true){
                    peso = -1000;
                }
            }
        }
        if (cadena == false){
            pp = make_pair(it, vLosNew[estoy.vSuce[0]].actual);
            double auxPeso = peso + Weights[pp] - solpi[it]; //modified weight
            if (insideBB == true){
                if (BranchingMethod == 1){
                    //Go over the node map ZERO
                    auto pos = Tree[VoyBB.first][VoyBB.second].NodeGroupZero.find(posMDD);
                    if (pos != Tree[VoyBB.first][VoyBB.second].NodeGroupZero.end()){
                        for (auto it2 = Tree[VoyBB.first][VoyBB.second].NodeGroupZero[posMDD].begin(); it2!= Tree[VoyBB.first][VoyBB.second].NodeGroupZero[posMDD].end(); it2++){
                            if (it == *it2){
                                auxPeso = -1000;
                                break;//found!
                            }
                        }
                    }
                }else{
                    //Go over the arc map ZERO
                    for (auto it2 = Tree[VoyBB.first][VoyBB.second].ArcGroupZero.begin(); it2!= Tree[VoyBB.first][VoyBB.second].ArcGroupZero.end(); it2++){
                        if (it == it2->first && vLosNew[estoy.vSuce[0]].actual == it2->second){
                            auxPeso = -1000;
                            break; //found!
                        }
                    }
                    //Go over the arc map ONE
                    bool allowed = false;
                    bool ifTail = false;
                    for (auto it2 = Tree[VoyBB.first][VoyBB.second].ArcGroupOne.begin(); it2!= Tree[VoyBB.first][VoyBB.second].ArcGroupOne.end(); it2++){
                        if (it == it2->first) ifTail = true;
                        if (it == it2->first && vLosNew[estoy.vSuce[0]].actual == it2->second){
                            allowed = true;
                            break;//found!
                        }
                    }
                    if (allowed == false && ifTail == true){
                        auxPeso = -1000;
                    }
                }
            }
            if (auxPeso > mmax) {
                tCicloFinAux = make_tuple(k, it, peso, 2, auxPeso);
                mmax = auxPeso;
                vLosNew[k].mActualCual = it;
            }
        }else{ // if it is a chain
            peso = vLosNew[k].MaxWeight + Weights[pp] - solpi[it];
            if (peso > mmax){
                mmax = peso;
                tCicloFinAux = make_tuple(k, it, 0, 0, peso);
                vLosNew[k].mActualCual = it;
                if (estoy.vTupleVarios.size() < CuantosPahtsQueres){
                    if (estoy.vTupleVarios.empty() == true){
                        estoy.vTupleVarios.push_back(make_tuple(k,it,peso));
                    } else {
                        if (peso > get<2>(estoy.vTupleVarios[0])){
                            estoy.vTupleVarios.clear();
                            estoy.vTupleVarios.push_back(make_tuple(k,it,peso));
                        }
                    }
                }
            } else if (estoy.vTupleVarios.size() < CuantosPahtsQueres and estoy.vTupleVarios.empty() == false){
                if (peso == get<2>(estoy.vTupleVarios[0])) estoy.vTupleVarios.push_back(make_tuple(k,it,peso));
            }
        }
    }
    return tCicloFinAux;
}

pair<int, double> Problem::predMax(vector<state>& vLosNew, int indexestoy, bool cadena, int CuantosPahtsQueres){
    pair<int, double> pM = make_pair(-1,   -10000000);//-100000000
    state estoy = vLosNew[indexestoy]; tCicloFin = make_tuple(-9, -9,   -10000000, -9,   -10000000);
    for (auto i = 0; i < estoy.vPrede.size(); i++){ // prede of estoy
        if (estoy.id == -8) { //when estoy es the last state
            tuple<int,int,double,int,double> tCicloFin2 = make_tuple(-1, -1,   -10000000, -1,   -10000000);
            tCicloFin2 = maxMenosOcho(vLosNew, estoy, i, cadena, CuantosPahtsQueres);
            if (get<4>(tCicloFin2) > get<4>(tCicloFin)){
                tCicloFin = tCicloFin2;
            }
        } else if (vLosNew[estoy.vPrede[i]].id == -1){ // if prede is origin
            pM.second = 0; //weight
            pM.first =  estoy.vPrede[i];
        } else if (estoy.id == -10){
            cout << "-10 You should not be here" << endl;
        } else {
            pair<int,int> pp = make_pair(vLosNew[estoy.vPrede[i]].actual, estoy.actual);
            double peso = -1000;
            if (cadena == false){
                peso = vLosNew[estoy.vPrede[i]].MaxWeight + Weights[pp] - solpi[vLosNew[estoy.vPrede[i]].actual]; // modified weight
            } else if (vLosNew[estoy.vPrede[i]].actual > Pairs - 1){
                peso = vLosNew[estoy.vPrede[i]].MaxWeight + Weights[pp] - solpi[vLosNew[estoy.vPrede[i]].actual] - solpi[estoy.actual]; // modified weight
            } else{
                peso = vLosNew[estoy.vPrede[i]].MaxWeight + Weights[pp] - solpi[estoy.actual]; // modified weight
            }
            if (insideBB == true){
                if (BranchingMethod == 1){
                    //Go over the node map ZERO
                    auto pos = Tree[VoyBB.first][VoyBB.second].NodeGroupZero.find(posMDD);
                    if (pos != Tree[VoyBB.first][VoyBB.second].NodeGroupZero.end()){
                        for (auto it2 = Tree[VoyBB.first][VoyBB.second].NodeGroupZero[posMDD].begin(); it2!= Tree[VoyBB.first][VoyBB.second].NodeGroupZero[posMDD].end(); it2++){
                            if (vLosNew[estoy.vPrede[i]].actual == *it2){
                                peso = -1000;
                                break;//found!
                            }
                        }
                    }
                }else{
                    //Go over the arc map ZERO
                    for (auto it2 = Tree[VoyBB.first][VoyBB.second].ArcGroupZero.begin(); it2!= Tree[VoyBB.first][VoyBB.second].ArcGroupZero.end(); it2++){
                        if (vLosNew[estoy.vPrede[i]].actual == it2->first && estoy.actual == it2->second){
                            peso = -1000;
                            break;//found!
                        }
                    }
                    //Go over the arc map ONE
                    bool allowed = false;
                    bool ifTail = false;
                    for (auto it2 = Tree[VoyBB.first][VoyBB.second].ArcGroupOne.begin(); it2!= Tree[VoyBB.first][VoyBB.second].ArcGroupOne.end(); it2++){
                        if (vLosNew[estoy.vPrede[i]].actual == it2->first) ifTail = true;
                        if (vLosNew[estoy.vPrede[i]].actual == it2->first && estoy.actual == it2->second){
                            allowed = true;
                            break;//found!
                        }
                    }
                    if (allowed == false && ifTail == true){
                        peso = -1000;
                    }
                }
            }
            if (peso > pM.second){
                pM.second = peso;
                pM.first =  estoy.vPrede[i];
                if (cadena == true && estoy.vTupleVarios.size() < CuantosPahtsQueres){
                    if (estoy.vTupleVarios.empty() == true){ // vector<tuple<int,int,double>>& vTupleVarios pos, actual, valor;
                        estoy.vTupleVarios.push_back(make_tuple(estoy.vPrede[i],vLosNew[estoy.vPrede[i]].actual,peso));
                    } else {
                        if (peso > get<2>(estoy.vTupleVarios[0])){
                            estoy.vTupleVarios.clear();
                            estoy.vTupleVarios.push_back(make_tuple(estoy.vPrede[i],vLosNew[estoy.vPrede[i]].actual,peso));
                        }
                    }
                }
            } else if (cadena == true && estoy.vTupleVarios.size() < CuantosPahtsQueres and estoy.vTupleVarios.empty() == false){
                if (peso == get<2>(estoy.vTupleVarios[0])) estoy.vTupleVarios.push_back(make_tuple(estoy.vPrede[i],vLosNew[estoy.vPrede[i]].actual,peso));
            }
        }
    }
    vLosNew[indexestoy] = estoy;
    return pM;
}

void Problem::FindLongestPath(){
    clock_t tStart = clock();
    tStart = clock();
    bool cadena = false;
    int ThisManyPaths = 15;
    cadena = ((vMDD[posMDD].first[1].actual) > Pairs - 1)? true : false;
    if (cadena == true){
        //cout << "chains...";
    }
    double max = -10000000; bool salir = false; int cadeLonPaMax = -1;
    auto itend = vMDD[posMDD].second.end(); //vLvEstados.end();
    //tuple<int,int,double,int,double> tCicloFin; //prede aux double prede double Only ciclo
    // 0: Position of predecessor state giving longest path of state -8
    // 1: "actual" of state -8 given by the predecessor at get<0>(tCicloFin)
    // 2: Maximum weight given by the max predecessor of -8
    // 3: Predecessor giving longest path of state -10
    // 4: Maximum weight given by the max predecessor state of -10
    //vector<tuple<int,int,double>> vTupleVarios; //pos actual double
    
    for (auto it = vMDD[posMDD].second.begin(); it != itend; it++){ // levels
        for (auto it2 = it->begin(); it2 != it->end(); it2++){ //states pos
            if (vMDD[posMDD].first[*it2].id != -88){
                vMDD[posMDD].first[*it2].vTupleVarios.clear();
                pair<int, double> pPredMax = predMax(vMDD[posMDD].first, *it2, cadena, ThisManyPaths);//predMax(vLosNew, int estoy, cadena);
                if (vMDD[posMDD].first[*it2].id == -8 && cadena == false){
                    vMDD[posMDD].first[2].MaxWeight = get<2>(tCicloFin);
                    vMDD[posMDD].first[2].predLongest = get<0>(tCicloFin); //position of state -8 is 2, -8 is the state before the last layer in a cycle
                    vMDD[posMDD].first[3].MaxWeight = get<4>(tCicloFin); //position of state -10 is 3, -10 is the last state at the last layer in a cycle
                    vMDD[posMDD].first[3].predLongest = get<3>(tCicloFin);
                    salir = true;
                    break;
                } else {    // cadena true and id != -8
                    if (vMDD[posMDD].first[*it2].id == -8){
                        vMDD[posMDD].first[*it2].MaxWeight = get<4>(tCicloFin);
                        vMDD[posMDD].first[*it2].predLongest = get<0>(tCicloFin);
                        pPredMax.second = get<4>(tCicloFin);
                        pPredMax.first = get<0>(tCicloFin);
                    }else{
                        vMDD[posMDD].first[*it2].MaxWeight = pPredMax.second;
                        vMDD[posMDD].first[*it2].predLongest = pPredMax.first;
                    }
                    if (pPredMax.second > max) {
                        max = pPredMax.second;
                        cadeLonPaMax = *it2;
                    }
                }
            }
        } //end for
        if (salir == true) break;
    }
    int maxDeOch = get<1>(tCicloFin);
    state it = vMDD[posMDD].first[cadeLonPaMax]; // end of chain
    if (cadena == true){
        if (it.MaxWeight <= 0) ThisManyPaths = 1;
    }else{
        ThisManyPaths = 1;
        it = vMDD[posMDD].first[3];
    }
    
    //once found, retrieve from bottom to top
    for (int cc = 0; cc < ThisManyPaths; cc++){
        vector<pair<int,double>> vLongest;
        if (cadena == true){
            it = vMDD[posMDD].first[cadeLonPaMax];
            if (cc >= it.vTupleVarios.size()) break;
            it.predLongest = get<0>(it.vTupleVarios[cc]);
            it.actual = get<1>(it.vTupleVarios[cc]);
            maxDeOch = it.actual;
            it.MaxWeight = get<2>(it.vTupleVarios[cc]);
        }
        while (it.id != -1){
            if (it.id == -8){
                auto it2 = it;
                it = vMDD[posMDD].first[it.predLongest]; // .first[it.predLongest];
                vLongest.push_back(make_pair(maxDeOch, it2.MaxWeight));
            } else {
                auto it2 = it;
                vLongest.push_back(make_pair(it.actual, it.MaxWeight));
                it = vMDD[posMDD].first[it.predLongest];
            }
        }
        reverse(vLongest.begin(), vLongest.end());
        LongestPath.push_back(vLongest);
        Verify();
    }
    LongestTime = (clock() - tStart)/double(CLOCKS_PER_SEC);
}

