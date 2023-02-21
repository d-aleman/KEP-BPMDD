//
//  AugmentedNetwork.cpp
//
//  Created by Carolina Riascos Alvarez on 2019-12-10.
//  Copyright © 2019 Carolina Riascos Alvarez. All rights reserved.

#include "AugmentedNetwork.hpp"
#include "Class_Problem_VFS.hpp"

random_device rd2;
mt19937 gen(rd2());

vector<int>::iterator nextVeciCade(state& esta, bool sillego, int lv, int k, map<int,bool>& mOrden, double porcentaje, vector<int>& vVeci){
    vector<int>::iterator resul;
    if (lv >= k or esta.vamosVeci/double(esta.totalVeci) >= porcentaje){
        return esta.ItEnd;
    }
    if (sillego == false && lv >= k){
        return esta.ItEnd;
    }else{
        if (esta.firstTime == true){
            esta.firstTime = false;
            esta.vamosVeci++;
            resul = esta.Itestoy ;
        } else {
            esta.vamosVeci++;
            resul = esta.Itestoy + 1;
        }
    }
    return resul;
}
vector<int>::iterator nextVeci(state& esta, bool sillego, int lv, int k, map<int,bool>& mOrden){
    vector<int>::iterator resul;
    if (lv >= k){
        return esta.ItEnd;
    }
    if (sillego == false && lv >= k){
        return esta.ItEnd;
    }else{
        if (esta.firstTime == true){
            esta.firstTime = false;
            resul = esta.Itestoy ;
        } else {
            resul = esta.Itestoy + 1;
        }
    }
    return resul;
}
pair<bool,bool> siVaFin(vector<int>::iterator itVeci, int origen, map<int,vector<int>>& Comp2, map<int,pair<bool,bool>>& mapaFin, int voyPeque, int cuantosMDDsPeque, int lv){
    pair<bool, bool> pareja;
    bool sisi = true;
    if (sisi == true){
        if (mapaFin.count(*itVeci) == 0){
            for (vector<int>::iterator it = Comp2[*itVeci].begin(); it != Comp2[*itVeci].end(); it++) {
                if (*it == origen) {
                    pareja.first = true;
                    if (lv == 4-1){
                        pareja.second = false;
                        return pareja;
                    }
                    if (Comp2[*itVeci].size() > 1) pareja.second = true;
                    mapaFin[*itVeci] = pareja;
                    return pareja;
                }
            }
            pareja.first = false;
            if (lv == 4-1){
                pareja.second = false;
                return pareja;
            }
            if (Comp2[*itVeci].size() > 0) pareja.second = true;
            mapaFin[*itVeci] = pareja;
            return pareja;
        } else {
            pair<bool,bool> aux = mapaFin[*itVeci];
            if (lv == 4-1){
                aux.second = false;
            }
            return aux;
        }
    }else{
        pareja.first = false;
        pareja.second = false;
        return pareja;
    }
    
}
pair<bool,bool> siVaFinCade(vector<int>::iterator itVeci, map<int,vector<int>>& Comp2, map<int,pair<bool,bool>>& mapaFin, int lv, int ChainLength){  //si va a origen, y si tiene más
    pair<bool, bool> pareja;
    if (lv == ChainLength){
        pareja.first = true;
        pareja.second = false;
    } else if (Comp2[*itVeci].size() > 0){
        pareja.second = true;
        pareja.first = false;
    } else {
        pareja.first = false;
        pareja.second = false;
    }
    return pareja;
}
bool yaVisitado (state& estoy, int candidato){
    bool resp;
    if (estoy.vVisited.count(candidato) == 0){
        resp = false;
    }else{
        resp = true;
    }
    return resp;
}
bool suceIguales (state& estoy, state& otro, vector<state>& vLosNew){ // comparación por posición
    vector<int> vEstoy, vOtro;
    for (auto it : estoy.vSuce){
        if (vLosNew[it].id == -8){
            if (vLosNew[it].mActual[estoy].size() == vLosNew[it].mActual[otro].size()){
                for (auto it2 : vLosNew[it].mActual[estoy]){
                    vEstoy.push_back(it2);
                }
            }else{
                return false;
            }
        } else {
            vEstoy.push_back(vLosNew[it].actual);
        }
    }
    for (auto it : otro.vSuce){
        if (vLosNew[it].id == -8){
            for (auto it2 : vLosNew[it].mActual[otro]){
                vOtro.push_back(it2);
            }
        } else {
            vOtro.push_back(vLosNew[it].actual);
        }
    }
    sort(vEstoy.begin(), vEstoy.end());
    sort(vOtro.begin(), vOtro.end());
    
    bool resp = vEstoy == vOtro;
    
    return resp;
}

//chains /////////////////////////////////////////
pair<vector<state>, vector<vector<int>>> Problem::BuildChainMDD(int iori, bool& hay, map<int,bool>& mOrden){
    int CuantosRegresosCero = 0;
    int origen = iori;
    map<int,vector<int>> Comp2 = Comp; // update comp for each copy
    
    vector<state> vLosNew;
    vector<vector<int>> vLvEstados;
    int lv = 1, id = 0;
    vLosNew.push_back(state(-1, -88, Comp2));
    vLosNew.push_back(state(id, origen, Comp2));
    //vLosNew.back().sizeBytes();
    int itNuevos = int(vLosNew.size()) - 1;
    vLosNew.back().vVisited[origen] = true;
    vLosNew.back().vPrede.push_back(0);
    vLosNew[0].vSuce.push_back(1);
    vLosNew.push_back(state(-8, -88, Comp2));
    int itDelFin = int(vLosNew.size()) - 1;
    
    vLvEstados.push_back(vector<int>());
    vLvEstados.back().push_back(itNuevos);
    map<int,pair<bool,bool>> mapaVaFin;
    map<int,bool> mapYaprede;
    int iLv = 0, jEs = 0; int jEsMer = 0;
    vector<vector<int>> vLvAfecta = vLvEstados; bool llegoMerge = false;
    
    double porceDeVeci = 0.18; // in percentage
    if (Nodes > 1500) porceDeVeci = 0.12;
    if (NDDs > 250) porceDeVeci = 0.08;
    int copyChainLength = 3;
    
    while (lv != 0){
        bool siLlego = false;
        auto it = nextVeciCade(vLosNew[vLvEstados[iLv][jEs]], siLlego, lv-1, int(copyChainLength), mOrden, porceDeVeci, Comp2[vLosNew[vLvEstados[iLv][jEs]].actual]);
        vLosNew[vLvEstados[iLv][jEs]].Itestoy = it;
        while(it != vLosNew[vLvEstados[iLv][jEs]].ItEnd){ // while there are neighbors
            if (yaVisitado(vLosNew[vLvEstados[iLv][jEs]], *it) == false) {
                pair<bool, bool> vaOriMas = siVaFinCade(it, Comp2, mapaVaFin, lv, int(copyChainLength));
                bool directo = false;
                if (vaOriMas.first == true){
                    vLosNew[vLvEstados[iLv][jEs]].cadeRepe = true;
                    vLosNew[itDelFin].mActual[vLosNew[vLvEstados[iLv][jEs]]].push_back(*it);
                    if (vLosNew[vLvEstados[iLv][jEs]].meConectoDos == false){
                        vLosNew[vLvEstados[iLv][jEs]].meConectoDos = true;
                        vLosNew[vLvEstados[iLv][jEs]].vSuce.push_back(itDelFin);
                    }
                    if (mapYaprede.count(vLvEstados[iLv][jEs]) == 0) {
                        mapYaprede[vLvEstados[iLv][jEs]] = true;
                        vLosNew[itDelFin].vPrede.push_back(vLvEstados[iLv][jEs]);
                    }
                    directo = true; siLlego = true; llegoMerge = true;
                    if (hay == false) hay = true;
                }
                if (vaOriMas.second == true && lv < copyChainLength) {
                    lv++; id++;
                    vLosNew.push_back(state(id, *it, Comp2));  // new node
                    itNuevos = int(vLosNew.size()) - 1;
                    vLosNew[itNuevos].vPrede.push_back(vLvEstados[iLv][jEs]);
                    vLosNew[itNuevos].vVisited =vLosNew[vLvEstados[iLv][jEs]].vVisited;
                    vLosNew[itNuevos].vVisited[vLosNew[itNuevos].actual] = true;
                    vLosNew[vLvEstados[iLv][jEs]].vSuce.push_back(itNuevos);
                    if (vLvEstados.size() < lv){
                        vLvEstados.push_back(vector<int>()); // create lv
                        vLvAfecta.push_back(vector<int>());
                    }
                    iLv++;
                    if (copyChainLength >= 3 && iLv >= 2 && iLv <= copyChainLength - 1){
                       JustVisitedByState[iLv][vLosNew[itNuevos].actual].push_back(itNuevos);
                    }
                    vLvEstados[iLv].push_back(itNuevos);
                    vLvAfecta[iLv].push_back(itNuevos);
                    jEs = int(vLvEstados[iLv].size()) - 1;
                    jEsMer = int(vLvAfecta[iLv].size() - 1);
                }
            } else {
                if (lv == int(copyChainLength) && it == vLosNew[vLvEstados[iLv][jEs]].ItEnd -1){
                    if (vLosNew[vLvEstados[iLv][jEs]].cadeRepe == false){
                        bool directo = false;
                        auto aaux = vLosNew[vLosNew[vLvEstados[iLv][jEs]].vPrede[0]];
                        vLosNew[itDelFin].mActual[aaux].push_back(vLosNew[vLvEstados[iLv][jEs]].actual);
                        if (vLosNew[vLvEstados[iLv][jEs]].meConectoDos == false){
                            vLosNew[vLvEstados[iLv][jEs]].meConectoDos = true;
                            vLosNew[vLvEstados[iLv][jEs]].vSuce.push_back(itDelFin);
                        }
                        if (mapYaprede.count(vLvEstados[iLv][jEs]) == 0) {
                            mapYaprede[vLvEstados[iLv][jEs]] = true;
                            vLosNew[vLvEstados[iLv][jEs]].meConectoDos = true;
                            vLosNew[itDelFin].vPrede.push_back(vLvEstados[iLv][jEs]);
                        }
                        directo = true; siLlego = false;
                        vLosNew[vLvEstados[iLv][jEs]].pathLlego = false;
                        if (hay == false) hay = true;
                        break;
                    }
                }
            }
            it = nextVeciCade(vLosNew[vLvEstados[iLv][jEs]], siLlego, lv-1, int(copyChainLength), mOrden, porceDeVeci, Comp2[vLosNew[vLvEstados[iLv][jEs]].actual]);
            vLosNew[vLvEstados[iLv][jEs]].Itestoy = it;
        }
        if (vLosNew[vLvEstados[iLv][jEs]].vSuce.size() == 0 && id != 0) {
            vLosNew[vLosNew[vLvEstados[iLv][jEs]].vPrede.back()].vSuce.pop_back();
            vLosNew.erase(vLosNew.begin() + vLvEstados[iLv][jEs]);
            vLvEstados[iLv].erase(vLvEstados[iLv].begin() + jEs);
            vLvAfecta[iLv].erase(vLvAfecta[iLv].begin() + jEsMer);
            id--;
        }
        lv--;
        iLv--;
        jEs = int(vLvEstados[iLv].size()) - 1;
        jEsMer = int(vLvAfecta[iLv].size() - 1);
    }
    vLvEstados.push_back(vector<int>(1,itDelFin));
    for (auto &i : vLosNew) i.eliminar88(vLosNew);
    return make_pair(vLosNew, vLvEstados);
}

//cycles /////////////////////////////////////////
pair<vector<state>, vector<vector<int>>> Problem::BuildCycleMDD(int iori, bool& hay, map<int,bool>& mOrden, int voyPeque, int cuantosMDDsPeque){
    int CuantosRegresosCero = 0;
    int origen = iori;
    map<int,vector<int>> Comp2 = Comp; // update comp for each copy, without the previous origin vertices
    if (mOrden.size() != 0){
        for (auto mi = 0; mi != Comp2.size(); mi++){
            auto mi2 = Comp2[mi].begin();
            for (; mi2 != Comp2[mi].end(); ) {
                if (mOrden.count(*mi2) == 1) {
                    mi2 = Comp2[mi].erase(mi2);
                } else {
                    mi2++;
                }
            }
        }
    }

    vector<state> vLosNew;
    vector<vector<int>> vLvEstados;
    int lv = 1, id = 0;
    vLosNew.push_back(state(-1, -88, Comp2));
    vLosNew.push_back(state(id, origen, Comp2));
    //vLosNew.back().sizeBytes();
    int itNuevos = int(vLosNew.size()) - 1;
    vLosNew.back().vVisited[origen] = true;
    vLosNew.back().vPrede.push_back(0);
    vLosNew[0].vSuce.push_back(1);  // suce origin
    vLosNew.push_back(state(-8, -88, Comp2));  // the end state
    int itDelFin = int(vLosNew.size()) - 1;  //iterator of the end state (now int)
    vLosNew.push_back(state(-10, origen, Comp2));
    int itDelFin2 = int(vLosNew.size()) - 1;
    vLosNew[int(vLosNew.size()) - 2].vSuce.push_back(int(vLosNew.size()) - 1);
    vLosNew[int(vLosNew.size()) - 1].vPrede.push_back(int(vLosNew.size()) - 2);
    
    vLvEstados.push_back(vector<int>());
    vLvEstados.back().push_back(itNuevos);
    map<int,pair<bool,bool>> mapaVaFin;
    map<int,bool> mapYaprede;
    int iLv = 0, jEs = 0; int jEsMer = 0; //by levels and states
    vector<vector<int>> vLvAfecta = vLvEstados; bool llegoMerge = false;
    while (lv != 0){
        bool siLlego = false;
        auto it = nextVeci(vLosNew[vLvEstados[iLv][jEs]], siLlego, lv, int(CycleLength), mOrden);
        vLosNew[vLvEstados[iLv][jEs]].Itestoy = it;
        while(it != vLosNew[vLvEstados[iLv][jEs]].ItEnd){ // while there are neighbors
            if (yaVisitado(vLosNew[vLvEstados[iLv][jEs]], *it) == false) {
                pair<bool, bool> vaOriMas = siVaFin(it, origen, Comp2, mapaVaFin, voyPeque, cuantosMDDsPeque, lv);
                bool directo = false;
                if (vaOriMas.first == true){
                    TotalCycles++;
                    vLosNew[itDelFin].mActual[vLosNew[vLvEstados[iLv][jEs]]].push_back(*it);
                    auto ut = NonPresentNodes.find(*it);
                    if (ut != NonPresentNodes.end()) NonPresentNodes.erase(*it);
                    if (vLosNew[vLvEstados[iLv][jEs]].meConectoDos == false){
                        vLosNew[vLvEstados[iLv][jEs]].meConectoDos = true;
                        vLosNew[vLvEstados[iLv][jEs]].vSuce.push_back(itDelFin);
                    }
                    if (mapYaprede.count(vLvEstados[iLv][jEs]) == 0) {
                        mapYaprede[vLvEstados[iLv][jEs]] = true;
                        vLosNew[itDelFin].vPrede.push_back(vLvEstados[iLv][jEs]);
                    }
                    directo = true; siLlego = true; llegoMerge = true;
                    if (hay == false) hay = true;
                }
                if (vaOriMas.second == true && lv < CycleLength - 1) {
                    lv++; id++;
                    vLosNew.push_back(state(id, *it, Comp2));  // new node
                    itNuevos = int(vLosNew.size()) - 1;  //iterator new node
                    vLosNew[itNuevos].vPrede.push_back(vLvEstados[iLv][jEs]);  //prede new state
                    vLosNew[itNuevos].vVisited =vLosNew[vLvEstados[iLv][jEs]].vVisited; // copy visited map
                    vLosNew[itNuevos].vVisited[vLosNew[itNuevos].actual] = true; // new state visited
                    vLosNew[vLvEstados[iLv][jEs]].vSuce.push_back(itNuevos);    // suce of the old state
                    if (vLvEstados.size() < lv){
                        vLvEstados.push_back(vector<int>()); // create lv
                        vLvAfecta.push_back(vector<int>());
                    }
                    iLv++;
                    vLvEstados[iLv].push_back(itNuevos);
                    vLvAfecta[iLv].push_back(itNuevos);
                    jEs = int(vLvEstados[iLv].size()) - 1;
                    jEsMer = int(vLvAfecta[iLv].size()) - 1;
                }
            }
            it = nextVeci(vLosNew[vLvEstados[iLv][jEs]], siLlego, lv, int(CycleLength), mOrden);
            vLosNew[vLvEstados[iLv][jEs]].Itestoy = it;
            while(it != vLosNew[vLvEstados[iLv][jEs]].ItEnd && *it == origen){
                it = nextVeci(vLosNew[vLvEstados[iLv][jEs]], siLlego, lv, int(CycleLength), mOrden);
                vLosNew[vLvEstados[iLv][jEs]].Itestoy = it;
            }
        }
        if (vLosNew[vLvEstados[iLv][jEs]].vSuce.size() == 0 && id != 0) {
            vLosNew[vLosNew[vLvEstados[iLv][jEs]].vPrede.back()].vSuce.pop_back();
            vLosNew.erase(vLosNew.begin() + vLvEstados[iLv][jEs]);
            vLvEstados[iLv].erase(vLvEstados[iLv].begin() + jEs);
            vLvAfecta[iLv].erase(vLvAfecta[iLv].begin() + jEsMer);
            id--;
        }
        lv--;   //go down one level
        iLv--;  //take iterator down
        jEs = int(vLvEstados[iLv].size()) - 1;
        jEsMer = int(vLvAfecta[iLv].size() - 1);
    }
    vLvEstados.push_back(vector<int>(1,itDelFin));
    vLvEstados.push_back(vector<int>(1,itDelFin2));
    for (auto &i : vLosNew) i.eliminar88(vLosNew);
    return make_pair(vLosNew, vLvEstados);
}

void Problem::BuildMDDs(){
    llenarComp();
    ordenFor = ordenDegre();
    clock_t tStart = clock();
    tStart = clock();
    map<int,bool> mOrden; bool cade = false;
    pair<vector<state>, vector<vector<int>>> aux;  bool hay = false; // aux is a copy
    int conta = 0;
    int cuantosMDDsPeque = 0;
    int contaCade = 0;
    int cuantosMDDsCade = NDDs;
    for (auto i = 0; i < ordenFor.size(); i++){//
        hay = false;
        if (ordenFor[i].first > Pairs - 1 && ChainLength > 0){ // altruistic node
            if (cade == false){
                cade = true;
            }
            if (ChainLength >= 3) JustVisitedByState = vector<map<int,vector<int>>>(ChainLength);// Three is fixed
            if (contaCade < cuantosMDDsCade){
                aux = BuildChainMDD(ordenFor[i].first, hay, mOrden);
                contaCade++;
                if(contaCade % 10 == 0) cout << "Chain number: " << contaCade << endl;
            }
        } else {  // cycle
            aux = BuildCycleMDD(ordenFor[i].first, hay, mOrden, conta, cuantosMDDsPeque);
        }
        if (hay == true){
            vMDD.push_back(aux);
            nstatesMDD+= aux.first.size();
            conta++;
        }
        if(i % 200 == 0) cout << "Number: " << i << endl;
        mOrden[ordenFor[i].first] = true;
    }
    MDDTime = (clock() - tStart)/double(CLOCKS_PER_SEC);
}
