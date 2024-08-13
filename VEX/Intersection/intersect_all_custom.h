// Attempting to create an analog of the intersect_all() function.
// This is a very simple implementation.
// Function identifies the nearest primitives to a segment using primfind().
// It then iterates over these primitives and applies the appropriate intersection function based on the primitive type.
// The function only works triangles, segments, and curves.


// intersection between segments
int intersect_segment_segment(vector p1, p2, p3, p4, P; float parameter[], EPS)
{

    // Early out: Check for parallel or co-linear segments
    vector dir1 = p2 - p1;
    vector dir2 = p4 - p3;
    float dot_prod = dot(dir1, dir2);
    float len_prod = length(dir1) * length(dir2);
    if (abs(dot_prod - len_prod) < EPS || abs(dot_prod + len_prod) < EPS) {
        return 0; // Parallel or co-linear
    }


   resize(parameter, 2);
   
   vector p13,p43,p21;
   float d1343,d4321,d1321,d4343,d2121;
   float numer,denom;
   
   vector pa, pb;
   
   p13.x = p1.x - p3.x;
   p13.y = p1.y - p3.y;
   p13.z = p1.z - p3.z;
   p43.x = p4.x - p3.x;
   p43.y = p4.y - p3.y;
   p43.z = p4.z - p3.z;
   
   if (abs(p43.x) < EPS && abs(p43.y) < EPS && abs(p43.z) < EPS)
      return 0;
   
   p21.x = p2.x - p1.x;
   p21.y = p2.y - p1.y;
   p21.z = p2.z - p1.z;
   
   if (abs(p21.x) < EPS && abs(p21.y) < EPS && abs(p21.z) < EPS)
      return 0;

   d1343 = p13.x * p43.x + p13.y * p43.y + p13.z * p43.z;
   d4321 = p43.x * p21.x + p43.y * p21.y + p43.z * p21.z;
   d1321 = p13.x * p21.x + p13.y * p21.y + p13.z * p21.z;
   d4343 = p43.x * p43.x + p43.y * p43.y + p43.z * p43.z;
   d2121 = p21.x * p21.x + p21.y * p21.y + p21.z * p21.z;

   denom = d2121 * d4343 - d4321 * d4321;
   if (abs(denom) < EPS)
      return 0;
   numer = d1343 * d4321 - d1321 * d4343;
    
   float mua = numer / denom;
   float mub = (d1343 + d4321 * (mua)) / d4343;
   
   if (mua < 0 || mua > 1 || mub < 0 || mub > 1) {
      pa = set(0,0,0);
      pb = set(0,0,0);
      return 0;
   }
   
   
   pa.x = p1.x + mua * p21.x;
   pa.y = p1.y + mua * p21.y;
   pa.z = p1.z + mua * p21.z;
   pb.x = p3.x + mub * p43.x;
   pb.y = p3.y + mub * p43.y;
   pb.z = p3.z + mub * p43.z;
   
   parameter[0] = mua;
   parameter[1] = mub;
   
   if (distance(pa, pb) < 0.001){ // if the distance is small - segments are intersecting
        P = avg(pa, pb);
        return 1;
    }else
        return 0;
}


int boxes_intersect(vector A_min, A_max, B_min, B_max)
{
    // Check overlap on x-axis
    if (A_max.x < B_min.x || B_max.x < A_min.x) {
        return 0; // No overlap
    }

    // Check overlap on y-axis
    if (A_max.y < B_min.y || B_max.y < A_min.y) {
        return 0; // No overlap
    }

    // // Check overlap on z-axis
    // if (A_max.z < B_min.z || B_max.z < A_min.z) {
    //     return 0; // No overlap
    // }

    // Overlap on all axes
    return 1; // Intersection
}

// segment and curve intersection
int intersect_segment_curve(int geo; int prim; vector p0, p1; vector positions[], uvws[]; float parameters[], eps){
    int pts[] = primpoints(geo, prim); // get curve points
    
    vector4 quat = dihedral(p1-p0, set(1,0,0));
    
    vector A_min = qrotate(quat, p0);
    vector A_max = qrotate(quat, p1);
    
    for(int i = 1; i < len(pts); i++){
        vector q0 = point(geo, "P", pts[i-1]);
        vector q1 = point(geo, "P", pts[i]);
       
        vector q0_rot = qrotate(quat, q0);
        vector q1_rot = qrotate(quat, q1);
        
        vector B_max = max(q0_rot, q1_rot);
        vector B_min = min(q0_rot, q1_rot);
        
        int isinside = boxes_intersect(A_min, A_max, B_min, B_max);
        if (!isinside) continue;
        
        vector P; float params[];
        int hasinters = intersect_segment_segment(p0, p1, q0, q1, P, params, eps);
        if (hasinters){
            append(positions, P);
            
            // calc u for segment points
            float u0 = (i-1) / (len(pts)*1.0-1);
            float u1 = i / (len(pts)*1.0-1);
            // use param for lerp between u0 and u1
            float u = lerp(u0, u1, params[1]);
            append(uvws, set(u, 0, 0));
            append(parameters, params[0]);
        }
        
    }
    return len(positions);
    
}



// segment and triangle intersection
int intersect_segment_triangle(int geo; int prim; vector P0, P1, P[], uvw[]; float parameter[], eps){

    int pts[] = primpoints(geo, prim);
    vector A = point(geo, 'P', pts[0]);
    vector B = point(geo, 'P', pts[1]);
    vector C = point(geo, 'P', pts[2]);
    
    int hasinters = 0;
    resize(P, 1);
    resize(uvw, 1);
    resize(parameter, 1);
    
    // Edge vectors for triangle
    vector AB = B - A;
    vector AC = C - A;
    
    // Plane normal
    vector normal = cross(AB, AC);
    
    // Line direction
    vector dir = P1 - P0;
    
    // Check if line and plane are parallel
    float dot_nd = dot(normal, dir);
    if (abs(dot_nd) < eps) {
        return 0; // No intersection
    }
    
    // Compute distance from P0 to plane
    float d = dot(normal, A - P0) / dot_nd;
    parameter[0] = d;
    
     // Check if intersection point is on the line segment
    if (d < 0.0 || d > 1.0) {
        return 0; // No intersection on the line segment
    }
    
    // Intersection point on line
    P[0] = P0 + d * dir;
    
    // Barycentric coordinates
    vector v0 = B - A;
    vector v1 = C - A;
    vector v2 = P[0] - A;
    
    float dot00 = dot(v0, v0);
    float dot01 = dot(v0, v1);
    float dot02 = dot(v0, v2);
    float dot11 = dot(v1, v1);
    float dot12 = dot(v1, v2);
    
    // Compute barycentric coordinates (u, v, w)
    float denom = dot00 * dot11 - dot01 * dot01;
    if (abs(denom) < eps) {
        return 0; // No intersection or degenerate triangle
    }
    
    float u = (dot11 * dot02 - dot01 * dot12) / denom;
    if (u < 0.0 || u > 1.0) {
        return 0; // Intersection outside triangle
    }
    
    float v = (dot00 * dot12 - dot01 * dot02) / denom;
    if (v < 0.0 || u + v > 1.0) {
        return 0; // Intersection outside triangle
    }
    uvw[0] = set(u, v, 0);
    
    // Intersection point is within the triangle
    return 1;
    
}

//////////////////////////////////////////////////////////////////////////////
// intersection with line (described by two points) and prim (triangle or line)
int get_intersection(int geo; int prim; vector P0, P1, P[], uvw[]; float parameter[]){
    int hasinters = 0;
    float eps = 1e-15;
    
    int closed = primintrinsic(geo, "closed", prim);
    if (closed){ // is a triangle
        hasinters = intersect_segment_triangle(geo, prim, P0, P1, P, uvw, parameter, eps);

    }else{ // is a segment or curve
        float params[];
        vector pos[];
        hasinters = intersect_segment_curve(geo, prim, P0, P1, P, uvw, parameter, eps);
    }
    return hasinters;
}


int intersect_all_custom(int geo; vector orig, dir; vector pos_res[]; int prim_res[]; vector uvw_res[]; float parameter_res[]){

    vector origins[] = array();
    vector dirs[] = array();
    
    vector P0 = orig;
    vector P1 = orig+dir;
    
    vector min = min(P0, P1);
    vector max = max(P0, P1);
    
    int potentialprims[] = primfind(geo, min, max);
    
    foreach(int pr; potentialprims){
        vector pos[], uvw[];
        float parameter[];
        int hasinters = get_intersection(geo, pr, P0, P1, pos, uvw, parameter);
        
        if (hasinters){
            append(prim_res, pr);
            append(pos_res, pos);
            append(uvw_res, uvw);
            append(parameter_res, parameter);
        }
    }
    
    return len(prim_res);

}

// version without parameter;
int intersect_all_custom(int geo; vector orig, dir; vector pos[]; int prims[]; vector uvws[]){
    float params[];
    return intersect_all_custom(geo, orig, dir, pos, prims, uvws, params);
}


