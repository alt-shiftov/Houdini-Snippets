// functions for working with lines (open polys)

//////////////////////////////////////////////////
// findprims() - Returns the primitive lines that connect two points.
// findprim() - Returns a single primitive line that connects two points.

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
// Removing Duplicated Lines (Lines that connect the same two points). 
// This is done by iterating over primitives and using findprims().

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


// Get neighboring primitives for a polyline (a line with two points).
function int[] polyline_neighbours(int geo; int primnum){
    int pts[] = primpoints(geo, primnum);
    
    int neis0[] = neighbours(geo, pts[0]);
    int neis1[] = neighbours(geo, pts[1]);
    removevalue(neis0, pts[1]);
    removevalue(neis1, pts[0]);
    
    int neiprims[];
    foreach(int nei; neis0){
        int edge = pointedge(geo, pts[0], nei);
        if (edge == -1) continue;
        
        int prim = hedge_prim(geo, edge);
        if (find(neiprims, prim) < 0){
            append(neiprims, prim);
        }
    }
    
    foreach(int nei; neis1){
        int edge = pointedge(geo, pts[1], nei);
        if (edge == -1) continue;
        
        int prim = hedge_prim(geo, edge);
        if (find(neiprims, prim) < 0){
            append(neiprims, prim);
        }
    }
    removevalue(neiprims, primnum);
    return neiprims;
    
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// Simple connectivity check for polylines with a constraint point group.
// Finds connected polylines "split" at constraint points.
// Works at the Detail level.
// Constraint points break the original line into segments.

// Uses polyline_neighbours() to find connected lines, 
// ignoring lines where in between is constraint points.

// Can be faster by using the Connectivity node to create initial classes,
// then iterating over the pieces instead of all geometry at the Detail level.

function int[] polyline_neighbours(int geo; int primnum; string constrgrp){
    int pts[] = primpoints(geo, primnum);
    
    int neis0[] = neighbours(geo, pts[0]);
    int neis1[] = neighbours(geo, pts[1]);
    removevalue(neis0, pts[1]);
    removevalue(neis1, pts[0]);
    
    if (inpointgroup(0, constrgrp, pts[0]))
        neis0 = array();
    
    if (inpointgroup(0, constrgrp, pts[1]))
        neis1 = array();
    
    int neiprims[];
    foreach(int nei; neis0){
        int edge = pointedge(geo, pts[0], nei);
        if (edge == -1) continue;
        
        int prim = hedge_prim(geo, edge);
        if (find(neiprims, prim) < 0){
            append(neiprims, prim);
        }
    }
    
    foreach(int nei; neis1){
        int edge = pointedge(geo, pts[1], nei);
        if (edge == -1) continue;
        
        int prim = hedge_prim(geo, edge);
        if (find(neiprims, prim) < 0){
            append(neiprims, prim);
        }
    }
    removevalue(neiprims, primnum);
    return neiprims;
    
}

int primCount = nprimitives(0);
int componentID = 0;  // Unique ID for each connected component

// Create an array to store component assignments for points
int ids[] = array();
resize(ids, primCount, -1);

// storing constraint group
string constrgroup = chs('constrgroup');

// Flood fill function to assign component ID
function void floodFill(int startPrim, compID; int ids[]; string constrptgroup) {
    
    int stack[] = array(startPrim);  // Stack to keep track of points to visit
    
    while(len(stack) > 0) {
        int curr = pop(stack);  // Get the last point from the stack
        if (ids[curr] != -1) continue;  // Skip if already visited

        ids[curr] = compID;  // Assign component ID to current point
        
        // Find neighboring points connected by edges
        int neighbors[] = polyline_neighbours(0, curr, constrptgroup);
        foreach(int nb; neighbors) {
            if (ids[nb] == -1) {  // Only visit unvisited neighbors
                append(stack, nb);  // Add neighbor to stack for future visits
            }
        }
    }
}


// Iterate over all points and apply flood fill where necessary
for (int prim = 0; prim < primCount; prim++) {
    if (ids[prim] == -1) {
        floodFill(prim, componentID, ids, constrgroup);  // Pass ids as an argument
        componentID++;  // Increment component ID for the next connected piece
    }
}

// Finally, assign the IDs back to points
i[]@ids = ids;


////////////////////////////////////////////
// After the Detail Wrangle, use a Primitive Wrangle to assign the class attribute
int ids[] = detail(0, "ids");
i@class = ids[@primnum];



//////////////////////////////////////////////
//////////////////////////////////////////////
// ### Near Line Primitives to Primitive - nearprims()
// Finds the nearest line primitives to a given primitive.
// This should be updated to utilize the pcsegment() function.

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
// Get curve uv by primitive and point

function float uvfrompoint(int geo; int ptnum; int primnum){
    int vtx = pointprimvertex(geo, ptnum, primnum);
    return vertexcurveparam(geo, vtx);
}



////////////////////////////////////
// ### UVs From Points - uvsfrompoints()
// Returns an array of all uv ranges that lie between the points

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

// RETURN ONE OR TWO UV RANGE FROM POINTS ARRAY
// return corner uvs. Will make two uvs, if it hits the edge (uvend = 1.0)

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