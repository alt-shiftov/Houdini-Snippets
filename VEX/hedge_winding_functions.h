// Hedges Functions:

// ------ iterating over hedges (3 vars)
// ------ iterating over hedge equiv's (2 vars)

// ------ compare_hedge_winding() - compare primitives winding by 2 equiv hedges
// ------ primshedges() - getting 2 hedges by 2 neighbour prims.

// compare_hedge_winding() + primedges() - This combination can be used to compare the winding (facing direction) of two prims.

// Winding Functions:
// ------ reverse_winding() - analog Reverse SOP node

// ------ computeWindingNumber2D() - Determines if a point is inside or outside a planar polygon.
// The polygon can be arbitrarily oriented in 3D space; it does not need to be perpendicular to any axis.

// ------ computeWindingNumber3D() - Determines if a point is inside or outside a solid 3D geometry.
// Developed for research purposes. I recommend to use official windingnumber() function.



///// HEDGE FUNCTIONS ////////////////////////////
///////////////////////////////////////////////////
// ### Iterating Over Hedges

// Variant 1
int hedge = primhedge(0, @primnum); //initialize hedge
int starthedge = hedge; // stash hedge
do{
    // ... do something
    hedge = hedge_next(0, hedge); // go to new hedge
}while(hedge != -1 && hedge != starthedge);


// Variant 2
int hedge = primhedge(0, @primnum); //initialize hedge
int starthedge = hedge; // stash hedge

while (hedge != -1){
    // ... do something
    
    hedge = hedge_next(0, hedge);
    if (hedge == starthedge) break; // exit while loop
}


// Variant 3
int start = primhedge(0, @primnum);

for (int hedge = start; hedge != -1; ){
    // ... do something
    
    hedge = hedge_next(0, hedge);
    if (hedge == start)
        break;
}



////////////////////////////////////////////////////
// ### Iterate Over Hedges And Equivs
// Iterate over hedges and all their equivs

// Variant 1
int starthedge = primhedge(0, @primnum); //pts[0];

for (int hedge = starthedge; hedge != -1; ){

    for (int nh = hedge_nextequiv(0, hedge); nh != hedge; nh = hedge_nextequiv(0, nh)){
        // ... do something
        
    }
    hedge = hedge_next(0, hedge);
    if (hedge == starthedge)
        break;
}



// Variant 2

int hedge = primhedge(0, @primnum); //pts[0]
int starthedge = hedge;

do{
    int equiv = hedge_nextequiv(0, hedge);
    int startequiv = equiv;
    
    for(int i = 0; i < hedge_equivcount(0, hedge); i++){
        if (equiv == -1 || equiv == hedge) break;
        
        // ... do something
        equiv = hedge_nextequiv(0, equiv);
    }

    hedge = hedge_next(0, hedge);
    
}while(hedge != -1 && hedge != starthedge);




// Compares the facing direction of neighboring polygons.
// Useful for determining whether the target polygon (neighbor) has the same orientation as the source polygon or if it is reversed.

// Returns:
// 1 - The winding is similar (polygons are facing in the same direction).
// 0 - The winding is different (polygons are facing in opposite directions).

function int compare_hedge_winding(int geo; int hedge0; int hedge1){

    int dst0 = hedge_dstpoint(geo, hedge0);
    int src0 = hedge_srcpoint(geo, hedge0);
    
    int dst1 = hedge_dstpoint(geo, hedge1);
    int src1 = hedge_srcpoint(geo, hedge1);
    
    if (dst0 == dst1 && src0 == src1){
        return 0;
    }
    return 1;
}


/////////////////////////////////////////////////////////

// Finds two hedges between two primitives
// If edge between prims does not exist - return [-1; -1]

// hedge[0] - prim0's hedge
// hedge[1] - prim1's hedge

function int[] primshedges(int geo; int prim0, prim1){
    int hedges[] = array();
    resize(hedges, 2, -1);
    int starthedge = primhedge(geo, prim0);

    for (int hedge = starthedge; hedge != -1; ){
        for (int nh = hedge_nextequiv(geo, hedge); nh != hedge; nh = hedge_nextequiv(geo, nh)){
            int neiprim = hedge_prim(geo, nh);
            if (neiprim == prim1){ // find prim
                hedges[0] = hedge;
                hedges[1] = nh;
                return hedges;
            }
            
        }
        hedge = hedge_next(geo, hedge);
        if (hedge == starthedge)
            break;
    }
    
    return hedges;
    
}




/////////////////////////////////////

// Reverse Facing of a polygon (like Reverse SOP)

function void reverse_winding(int geo; int prims){
    
    int primverts[] = primvertices(geo, prims);
    int primpts[] = reverse(primpoints(geo, prims)); //reverse array of points
    
    for(int i = 0; i < len(primverts); i++){
        //relink vertices to points (reverse array)
        setvertexpoint(geo, -1, primverts[i], primpts[i]);
    }
}

///////////////////////////////////////////////////////////////////
// Compute Winding for a Point Position on 2D Polygon

// Determines if a point is inside or outside a planar polygon


// by point and vertices of the polygon
function int computeWindingNumber2D(vector pt, vertices[]) {

    int windingNumber = 0;
    int numVertices = len(vertices);
    
    // Calculate the normal of the polygon's plane using the first three vertices
    vector normal = normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));

    // Choose two vectors on the plane to form a local 2D coordinate system (u, v)
    vector u = normalize(vertices[1] - vertices[0]);
    vector v = normalize(cross(normal, u));

    // Project the point and vertices onto the 2D plane defined by (u, v)
    vector2 proj_pt = set(dot(pt - vertices[0], u), dot(pt - vertices[0], v));
    vector2 proj_vertices[];

    foreach(vector vert; vertices) {
        vector2 vtxpos = set(dot(vert - vertices[0], u), dot(vert - vertices[0], v));
        push(proj_vertices, vtxpos);
    }

    // 2D winding number algorithm
    for (int i = 0; i < numVertices; i++) {
        vector2 v1 = proj_vertices[i];
        vector2 v2 = proj_vertices[(i + 1) % numVertices];

        if (proj_pt == v1 || proj_pt == v2) {
            return 0; // Point lies on the edge, winding number is 0
        }

        float dy = v2.y - v1.y;

        if ((v1.y <= proj_pt.y && proj_pt.y < v2.y) || (v2.y <= proj_pt.y && proj_pt.y < v1.y)) {
            float intersectX = v1.x + (proj_pt.y - v1.y) * (v2.x - v1.x) / dy;

            if (proj_pt.x < intersectX) {
                windingNumber += (v2.y > v1.y) ? 1 : -1;
            }
        }
    }

    return windingNumber;
}


// by point and primitve number
function int computeWindingNumber2D(int geo, prim; vector pt) {
    int pts[] = primpoints(geo, prim);
    vector vertices[];
    foreach(int p; pts){
        vector pos = point(geo, "P", p);
        append(vertices, pos);
    }
    
    return computeWindingNumber2D(pt, vertices);
    
}


// Example
i@winding = computeWindingNumber2D(1, 0, v@P);




//////////////////////////////////////////////////
// Compute Winding For Point inside Solid Geometry

// It is a very obscure method. Made it for studying. Required trianglulated solid geometry
// Better to use windingnumber() function

// Formula is from: The Solid Angle of a Plane Triangle
// https://ieeexplore.ieee.org/document/4121581/authors#authors


// calculate the solid angle of the triangle
float solidAngle(vector P, A, B, C) {
    vector R1 = A - P;
    vector R2 = B - P;
    vector R3 = C - P;
    
    float R1m = length(R1);
    float R2m = length(R2);
    float R3m = length(R3);
    
    float numerator = dot(R1, cross(R2,R3));
    
    float denominator = R1m*R2m*R3m 
                        + dot(R1,R2)*R3m 
                        + dot(R1,R3)*R2m 
                        + dot(R2,R3)*R1m;

    return 2 * atan2(numerator, denominator);

}

// getting point positions of every triangle from solid geo
function float computeWindingNumber3D(int geo; vector pt) {

    vector faces[];
    for(int i = 0; i < nprimitives(geo); i++){
        int pts[] = primpoints(geo, i);
        foreach(int p; pts){
            vector pos = point(geo, "P", p);
            append(faces, pos);
        }
    }

    float totalSolidAngle = 0.0;

    for(int i = 0; i < len(faces); i+=3){
        // Assume each face is a triangle; if not, you need to triangulate it
        vector A = faces[i];
        vector B = faces[i+1];
        vector C = faces[i+2];
        
        // Sum the solid angles
        totalSolidAngle += solidAngle(pt, A, B, C);
    }

    // Calculate the winding number. Reversing, because Houdini is using a clocwise order
    float windingNumber = - (totalSolidAngle) / (4 * $PI);
    return windingNumber;
}

// using function
f@winding = computeWindingNumber3D(1, v@P);




