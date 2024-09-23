// Near Edges to Edge - nearedgestoedge()
// Finds all edges that are close to the specified edge. 
// The function nearedgetoedge() returns a single edge.

function int[] nearedgestoedge(int geo; string edgegroup; int pt01, pt02; float dist){
    if (edgegroup == '*') edgegroup = '';

    vector pos01 = point(geo, "P", pt01);
    vector pos02 = point(geo, "P", pt02);
    int pts01[] = pcfind(geo, "P", pos01, dist, 1e15);
    int pts02[] = pcfind(geo, "P", pos02, dist, 1e15);

    removevalue(pts01, pt01);
    removevalue(pts01, pt02);

    removevalue(pts02, pt01);
    removevalue(pts02, pt02);
    
    int nearedges[] = array();
    foreach(int p01; pts01){
        foreach(int p02; pts02){
            if (inedgegroup(geo, edgegroup, p01, p02) && p01 > p02){
                append(nearedges, p01);
                append(nearedges, p02);
            }
        }
    }
    return nearedges;

}


function int[] nearedgestoedge(int geo; int pt01, pt02; float dist){
    string edgegroup = '';
    int edges[] = nearedgestoedge(geo, edgegroup, pt01, pt02, dist);
    
    return edges;
}

function int[] nearedgetoedge(int geo; string edgegroup; int pt01, pt02; float dist){
    int edge[] = nearedgestoedge(geo, edgegroup, pt01, pt02, dist);
    if (len(edge) > 0)
        return array(edge[0], edge[1]);
    else
        return array();
}

function int[] nearedgetoedge(int geo; int pt01, pt02; float dist){
    string edgegroup = '';
    int edge[] = nearedgestoedge(geo, edgegroup, pt01, pt02, dist);
    if (len(edge) > 0)
        return array(edge[0], edge[1]);
    else
        return array();
}


////////////////////////////////////
// Near Edges to Point - nearedgestopoint()
// Finds all edges that are close to a specified position.

function int[] nearedgestopoint(int geo; string edgegroup; vector ptpos; float ds[]){
    int prim = -1;
    vector primuv = 0;
    xyzdist(geo, ptpos, prim, primuv);
    if (edgegroup == '*') edgegroup = '';
    
    vector pos = primuv(geo, "P", prim, primuv);
    int pts[] = primpoints(geo, prim);
    
    float dists[] = array();
    int targetpts[] = array();
    
    foreach(int pt; pts){
        vector pos01 = point(geo, "P", pt);
        foreach(int nei; neighbours(geo, pt)){
            if (find(pts, nei) < 0) continue;
            if (pt < nei) continue;
            
            if (inedgegroup(geo, edgegroup, pt, nei)){
                vector pos02 = point(geo, "P", nei);
                float d = ptlined(pos01, pos02, pos);
                append(dists, d);
                append(targetpts, array(pt, nei));
            }
        }
    }
    int indices[] = argsort(dists);
    int targetpts_sort[] = array();
    
    if (len(dists) > 0){
        ds = reorder(dists, indices);
        for(int i = 0; i < len(pts); i++){
            int index = indices[i];
            append(targetpts_sort, array(targetpts[index*2], targetpts[index*2+1]));
        }
        return targetpts_sort;
    }else{ 
        return array();
    }
}


function int[] nearedgestopoint(int geo; vector ptpos; float dists[]){
    string edgegrp = '';
    int pts[] = nearedgestopoint(geo, edgegrp, ptpos, dists);
    return pts;
}

function int[] nearedgestopoint(int geo; vector ptpos){
    string edgegrp = '';
    float dists[] = array();
    int pts[] = nearedgestopoint(geo, edgegrp, ptpos, dists);

    return pts;
}

function int[] nearedgestopoint(int geo; string edgegroup; vector ptpos){
    float dists[] = array();
    int pts[] = nearedgestopoint(geo, edgegroup, ptpos, dists);
    
    return pts;
}

function int nearedgetopoint(int geo; string edgegroup; vector ptpos){
    float dists[] = array();
    int pts[] = nearedgestopoint(geo, edgegroup, ptpos, dists);
    if (len(pts) > 0)
	    return pts[0];
	 else
		 return -1;
}


//////////////////////////////////////////////////////////

// Returns all primitives that are connected to a specified edge.
function int[] edgeprims(int geo; int pt01, pt02){
    int prims[] = array();
    
    int hedge = pointedge(geo, pt01, pt02);
    int starthedge = hedge;
    
    if (hedge == -1) return prims;
    do{
        int prim = hedge_prim(geo, hedge);
        if (find(prims, prim) < 0)
            append(prims, prim);
            
        hedge = hedge_nextequiv(geo, hedge);
    }while(hedge > -1 && hedge != starthedge);
    
    return prims;
    
}



///////////////////////////////////////////////
// Returns all edges associated with a primitive.
// Each edge is defined by two points.
function int[] primedges(int geo; int prim){
    
    int pts[] = array();
    
    int hedge = primhedge(geo, prim);
    if (hedge == -1) return pts;
    int hedges[] = array();
    
    int starthedge = hedge;
    do{
        int pt01 = hedge_srcpoint(geo, hedge);
        int pt02 = hedge_dstpoint(geo, hedge);
        
        if (find(hedges, hedge) < 0){
            append(hedges, hedge);
            append(pts, pt01);
            append(pts, pt02);
        }
            
        hedge = hedge_next(geo, hedge);
    }while (hedge != starthedge && hedge > -1);
    
    return pts;
}



///////////////////////////////////////
// Point Edges - pointedges()
// Returns the edges that are connected to a specified point.

function int[] pointedges(int geo; int ptn){
    int pts[] = neighbours(geo, ptn);
    int edges[] = array();
    foreach(int pt; pts){
        append(edges, ptn);
        append(edges, pt);
    }
    return edges;
}



////////////////////////////////////////////
// Is Edge In Prim
// Checks if an edge is part of a primitive's edge group.
// Returns the two points that define the edge.

function int isedgeinprim(int geo; string edgegroup; int prim; int pts[]){
    int hedge = primhedge(geo, prim);
    if (hedge == -1) return 0;
    
    int starthedge = hedge;
    do{
        int pt01 = hedge_srcpoint(geo, hedge);
        int pt02 = hedge_dstpoint(geo, hedge);
        
        if (inedgegroup(geo, edgegroup, pt01, pt02)){
            append(pts, pt01);
            append(pts, pt02);
            return 1;
        }
            
        hedge = hedge_next(geo, hedge);
    }while (hedge != starthedge && hedge > -1);
    
    return 0;
}

function int isedgeinprim(int geo; string edgegroup; int prim){
    int pts[] = array();
    int res = isedgeinprim(geo, edgegroup, prim, pts);
    return res;
}





