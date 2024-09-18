// functions for working with lines (open polys)

//////////////////////////////////////////////////
// findprims - finding primitive lines by two points
// findprim - finding one primitive line by two points

function int[] findprims(int geo; string primgrp; int pt01; int pt02){    
    int prims[];
    
    int hedge = pointhedge(geo, pt01, pt02); //pts[0];
    if (hedge == -1)
        hedge = pointhedge(geo, pt02, pt01); 
    
    if (hedge == -1) return prims;
    
    int prim = hedge_prim(geo, hedge);
    if (prim >= 0) append(prims, prim);
    
    for (int nh = hedge_nextequiv(geo, hedge); nh != hedge; nh = hedge_nextequiv(geo, nh)){
        prim = hedge_prim(geo, nh);
        if (prim >= 0)
            append(prims, prim);
        
    }

    return prims;
}


// getting prims without primgroup
function int[] findprims(int geo; int pt01; int pt02){
    
    string primgrp = '';
    int prims[] = findprims(geo, primgrp, pt01, pt02);
    return prims;
}

// finding one prim with primgroup
function int findprim(int geo; string primgrp; int pt01; int pt02){

    int prims[] = findprims(geo, primgrp, pt01, pt02);
    
    if (len(prims) == 0)
        return -1;
    else
        return prims[0];
}

// finding one prim without primgroup
function int findprim(int geo; int pt01; int pt02){
    string primgrp = '';
    int prims[] = findprims(geo, primgrp, pt01, pt02);
    if (len(prims) == 0)
        return -1;
    else
        return prims[0];
}



/////////////////////////////////////////
// Removing Duplicated Lines. Iterating over primitives, using findprims()

function void removeduplines(int geo; int primnum){

    int pts[] = primpoints(geo, primnum);
    int prims[] = findprims(geo, pts[0], pts[1]);

    prims = sort(prims);
    if (primnum != prims[0]) return;

    foreach(int prim; prims){
        if (prim == primnum) continue;
        removeprim(geo, prim, 0);
    }
}


//////////////////////////////////////////////
// ### Near Line Prims to Prim - nearprims()
// Находит все линии, близкие к искомой линии
// мб заюзать лучше pcsegment()
// find near line prims 

function int[] nearprims(int geo; int prim; float rad){
    int nearprims[] = array();
    int pts[] = primpoints(geo, prim);
    vector pos01 = point(geo, "P", pts[0]);
    vector pos02 = point(geo, "P", pts[1]);

    int num = 1e15;
    int pts01[] = pcfind(geo, "P", pos01, rad, num);
    int pts02[] = pcfind(geo, "P", pos02, rad, num);
    
    foreach(int pt01; pts01){
        // if (pt01 == pts[0]) continue;
        foreach(int pt02; pts02){
            // if (pt02 == pts[1]) continue;
            int nprims[] = findprims(geo, pt01, pt02);
            foreach(int nprim; nprims){
                if (nprim >= 0 && prim != nprim){
                    append(nearprims, nprim);
                }
            }
            
        }
    }
    return nearprims;
}


// VERTICES
////////////////////////////////////////
// Vertex by Point and Prm - pointprimvertex()

function int pointprimvertex(int geo; int ptnum; int primnum){
    int vtxs01[] = pointvertices(geo, ptnum);
    int vtxs02[] = primvertices(geo, primnum);

    foreach(int vtx; vtxs01){
        if (find(vtxs02, vtx) >= 0){
            return vtx;
        }
    }
    return -1;
}


///////////////////////////////////////////////
// Near Vertices by UV - nearverticesbyuv()
// 
function int[] nearverticesbyuv(int geo; float targetuv; int primnum; float maxuvdist){
    int vtxs[] = primvertices(geo, primnum);

    int linearvtxs[] = array();
    float dists[] = array();
    
    foreach(int vtx; vtxs){
        int linearvtx = primvertex(geo, primnum, vtx);
        float uv = vertexcurveparam(geo, linearvtx);
        float dist = abs(targetuv - uv);
        
        if (dist <= maxuvdist){
            append(linearvtxs, linearvtx);
            append(dists, dist);
        }
    }
    int idxs[] = argsort(dists);

    return reorder(linearvtxs, idxs);
}

function int nearvertexbyuv(int geo; float targetuv; int primnum; float maxuvdist){
    int vtxs[] = nearverticesbyuv(geo, targetuv, primnum, maxuvdist);
    if (len(vtxs) > 0)
        return vtxs[0];
    else
        return -1;
}

function int nearvertexbyuv(int geo; float targetuv; int primnum){
    float maxuv = 2.0;
    int vtx = nearvertexbyuv(geo, targetuv, primnum, maxuv);
    return vtx;
}



///////////////////////////////////////////
// Points array from UV - pointsfromuv()
// GET ARRAY OF POINTS FROM U parametric values (START/END)

function void pointsfromuv(int geo; int primnum; vector2 uvw; int pts[]){
    int vtxs[] = primvertices(geo, primnum);
    float uvs[] = array(clamp(uvw[0], 0, 1), clamp(uvw[1], 0, 1));
    uvs = sort(uvs);
    foreach(int vtx; vtxs){
        float uv = vertexcurveparam(geo, vtx);
        int pt = vertexpoint(geo, vtx);
        if (uv >= uvs[0] && uv <= uvs[1] && find(pts, pt) < 0){
            append(pts, pt);
        }
    }
}

// SEPARATE UVW (TWO FLOATS)
function void pointsfromuv(int geo; int primnum; float start; float end; int pts[]){
    vector2 uv = set(start, end);
    pointsfromuv(geo, primnum, uv, pts);
}

// RETURN ARRAY OF POINTS
function int[] pointsfromuv(int geo; int primnum; vector2 uvw){
    int pts[];
    pointsfromuv(geo, primnum, uvw, pts);
    return pts;
}

function int[] pointsfromuv(int geo; int primnum; float start; float end){
    int pts[];
    vector2 uvw = set(start, end);
    pointsfromuv(geo, primnum, uvw, pts);
    return pts;
}



///////////////////////////////////
// ### UV From Points - uvfrompoint()
// Берет uv по номеру примитива и точки

function float uvfrompoint(int geo; int ptnum; int primnum){
    int vtx = pointprimvertex(geo, ptnum, primnum);
    return vertexcurveparam(geo, vtx);
}



////////////////////////////////////
// ### UVs From Points - uvsfrompoints()
// Выдает массив всех диапазонов uv, которые лежат между точками


function void uvsfrompoints_all(int geo; int pts[]; int primnum; vector2 uvranges[]; string ptgrps[]){
    int allpts[] = primpoints(geo, primnum);
    if (allpts[-1] == allpts[0]) // remove duplicated index from end
        removeindex(allpts, -1);
    float uvstart = 2.0;
    float uvend = -1.0;
    string ptgrp = '';
    foreach(int idx; int pt; allpts){
        float uv = uvfrompoint(geo, pt, primnum);
        if (find(pts, pt) >= 0){
            if (uvstart > uv || uvend < uv){
                ptgrp += itoa(pt) + ',';
                if (uvstart > uv) uvstart = uv;
                if (uvend < uv) uvend = uv;
            }
        }else{
            if (uvstart != 2.0 && uvend != -1.0){
                append(uvranges, set(uvstart, uvend));
                
                if (ptgrp[-1] == ',')
                    ptgrp = ptgrp[:-1];
                append(ptgrps, ptgrp);
            }
            uvstart = 2.0;
            uvend = -1.0;
            ptgrp = '';
        }
        
    }
    if (uvstart != 2.0 && uvend != -1.0)
        append(uvranges, set(uvstart, uvend));
    
}

// GET ONE OR TWO UV RANGE FROM POINTS ARRAY - находит самые крайние uv. СДелает два uv, если упрется в край (uvend = 1.0)
function void uvsfrompoints_main(int geo; int pts[]; int primnum; vector2 uvranges[]; string ptgrps[]){
    int allpts[] = primpoints(geo, primnum);
    if (allpts[-1] == allpts[0]) // remove duplicated index from end
        removeindex(allpts, -1);
    float uvstart = 2.0;
    float uvend = -1.0;
    string ptgrp = '';
   
    foreach(int idx; int pt; allpts){
        float uv = uvfrompoint(geo, pt, primnum);
        if (find(pts, pt) >= 0){
            if (uvstart > uv || uvend < uv){
                ptgrp += itoa(pt) + ',';
                if (uvstart > uv) uvstart = uv;
                if (uvend < uv) uvend = uv;
            }
        }
        
        if (uvend == 1.0){
            append(uvranges, set(uvstart, uvend));
            
            if (ptgrp[-1] == ',')
                ptgrp = ptgrp[:-1];
            append(ptgrps, ptgrp);
            
            uvstart = 2.0;
            uvend = -1.0;
            ptgrp = '';
        }
            
    }
    if (uvstart != 2.0 && uvend != -1.0)
        append(uvranges, set(uvstart, uvend));
    
}

// NO MODE
function void uvsfrompoints(int geo; int pts[]; int primnum; vector2 uvranges[]; string ptgrps[]){
    uvsfrompoints_main(geo, pts, primnum, uvranges, ptgrps);
}

// WITH MODE
function void uvsfrompoints(int geo; int pts[]; int primnum; vector2 uvranges[]; string ptgrps[]; string mode){
    if (mode == 'all')
        uvsfrompoints_all(geo, pts, primnum, uvranges, ptgrps);
    else
        uvsfrompoints_main(geo, pts, primnum, uvranges, ptgrps);
}

////// RETURN UVRANGES[]
function vector2[] uvsfrompoints(int geo; int pts[]; int primnum; string ptgrps[]; string mode){
    vector2 uvranges[];
    if (mode == 'all')
        uvsfrompoints_all(geo, pts, primnum, uvranges, ptgrps);
    else
        uvsfrompoints_main(geo, pts, primnum, uvranges, ptgrps);
    return uvranges;
}

// NO MODE
function vector2[] uvsfrompoints(int geo; int pts[]; int primnum; string ptgrps[]){
    vector2 uvranges[];
    string mode = 'main';
    if (mode == 'all')
        uvsfrompoints_all(geo, pts, primnum, uvranges, ptgrps);
    else
        uvsfrompoints_main(geo, pts, primnum, uvranges, ptgrps);
    return uvranges;
}

/////// RETURN UVRANGES[] WITHOUT PTGRPS[]
function vector2[] uvsfrompoints(int geo; int pts[]; int primnum; string mode){
    vector2 uvranges[];
    string ptgrps[];
    if (mode == 'all')
        uvsfrompoints_all(geo, pts, primnum, uvranges, ptgrps);
    else
        uvsfrompoints_main(geo, pts, primnum, uvranges, ptgrps);
    return uvranges;
}

// NO MODE
function vector2[] uvsfrompoints(int geo; int pts[]; int primnum){
    vector2 uvranges[];
    string ptgrps[];
    string mode = 'main';
    if (mode == 'all')
        uvsfrompoints_all(geo, pts, primnum, uvranges, ptgrps);
    else
        uvsfrompoints_main(geo, pts, primnum, uvranges, ptgrps);
    return uvranges;
}