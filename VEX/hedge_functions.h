// ### Prims Hedge - primshedge()
// Находит hedge между двумя примитивами. Выдает hedge который принадледит prim

function int primshedge(int geo; int prim; int primnei){
    int hedge = primhedge(geo, prim);
    if (hedge == -1) return -1;
    
    int starthedge = hedge;
    int hedgeresult = -1;
    do{
        int equiv = hedge_nextequiv(geo, hedge); 
        int p = hedge_prim(geo, equiv);
        if (p == primnei){
            hedgeresult = hedge;
            return hedgeresult;
        }
        hedge = hedge_next(geo, hedge);
    }while(hedge != starthedge && hedgeresult == -1);
    
    return hedgeresult;
    
}