// resamplebylength() - resample function, result is similar to the Resample SOP.
// The simple version generates the @sourceprim and @sourceprimuv attributes, 
// which are useful for interpolating attributes from the original curve to the new resampled points.

// resamplebylength_pin() - This version allows resampling a curve while keeping certain points fixed (pinned).
// Pinned points are defined by an integer array of point indices (pinpts[]).
// Pinned points will remain in place during resampling.

// uv01 - starting UV parameter on the curve from which the resampling begins.
// uv02 - ending UV parameter on the curve where the resampling stops.
// prim - primitive curve number.
// minlen - minimum segment length
// newprim - output primitive number
// newpts[] - output array of new point indices
// uvs[] - output UV coordinates of new points


// BASE function
function void resamplebylength(int geo; float uv01, uv02; int prim; float minlen; int newprim; int newpts[]; float uvs[]){
    int divs = 4;
    float dist = primarclen(geo, set(uv01, 0), set(uv02, 0), prim, divs, PRIMUV_UNIT_TO_REAL);
    int num = int(dist / (minlen*1.0));
    if (minlen > dist) num = 1;
    
    for(int n = 0; n <= num; n++){
        float uv = fit(n*1.0, 0.0, num, uv01, uv02);
        if (find(uvs, uv) >= 0) continue;
        
        vector pos = primuv(geo, "P", prim, set(uv, 0, 0));
        int pt = addpoint(0, pos);
        addvertex(0, newprim, pt);
        append(newpts, pt);
        append(uvs, uv);
    }
}


// RESAMPLE WITH PIN UVs. VOID
function void resamplebylength_pin(int geo; float uv01, uv02; int prim; float minlen; int newprim; float pinuvs[]; int newpts[]; float newuvs[]){
    // before pin points
    resamplebylength(geo, uv01, pinuvs[0], prim, minlen, newprim, newpts, newuvs);

    // between pin points
    for(int i = 0; i < len(pinuvs)-1; i++){
        int nexti = (i+1) % len(pinuvs);
        resamplebylength(geo, pinuvs[i], pinuvs[nexti], prim, minlen, newprim, newpts, newuvs);
    }
    // after pin points
    resamplebylength(geo, pinuvs[-1], uv02, prim, minlen, newprim, newpts, newuvs);
}

// RESAMPLE WITH PIN UVs. RETURN NEW PRIM
function int resamplebylength_pin(int geo; float uv01, uv02; int prim; float minlen; float pinuvs[]; int newpts[]; float newuvs[]){
    int newprim = addprim(geo, 'polyline');
    resamplebylength_pin(geo, uv01, uv02, prim, minlen, newprim, pinuvs, newpts, newuvs);
    return newprim;
}

// get vertex by prim and point
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

// RESAMPLE WITH PIN POINTS. RETURN NEW PRIM
function int resamplebylength_pin(int geo; float uv01, uv02; int prim; float minlen; int pinpts[]; int newpts[]; float newuvs[]){
    int newprim = addprim(geo, 'polyline');

    float pinuvs[];
    foreach(int pt; pinpts){
        int vtx = pointprimvertex(geo, pt, prim);
        float uv = vertexcurveparam(geo, vtx);
        append(pinuvs, uv);
    }

    resamplebylength_pin(geo, uv01, uv02, prim, minlen, newprim, pinuvs, newpts, newuvs);
    return newprim;
}


// This function assigns the @sourceprim and @sourceprimuvw attributes to the newly generated points.
// It is useful for transferring attributes from the original primitive to the new one.
function void setuvtopoints(int geo; int pts[]; int primnum; float uvs[]; string primchannel, uvchannel){
    foreach(int idx; int pt; pts){
        float u = uvs[idx];
        setpointattrib(geo, primchannel, pt, primnum);
        setpointattrib(geo, uvchannel, pt, set(u, 0, 0));
    }
}

///////////////////////////////////////////////////////////////////////////////////

// SIMPLE VERSION OF RESAMPLE
function int resamplebylength(int geo; int prim; float minlen){
    int newprim = addprim(geo, 'polyline');
    int newpts[] = array();
    float newuvs[] = array();
    resamplebylength(geo, 0.0, 1.0, prim, minlen, newprim, newpts, newuvs);
    setuvtopoints(geo, newpts, prim, newuvs, 'sourceprim', 'sourceprimuv');
    
    return newprim;
}


// SIMPLE VERSION OF RESAMPLE WITH PIN POINTS
function int resamplebylength_pin(int geo; float uv01, uv02; int prim; float minlen; int pinpts[]){
    int newprim = addprim(geo, 'polyline');

    float pinuvs[]; 
    int newpts[];
    float newuvs[];
    foreach(int pt; pinpts){
        int vtx = pointprimvertex(geo, pt, prim);
        float uv = vertexcurveparam(geo, vtx);
        append(pinuvs, uv);
    }

    resamplebylength_pin(geo, uv01, uv02, prim, minlen, pinuvs, newpts, newuvs);
    setuvtopoints(geo, newpts, prim, newuvs, 'sourceprim', 'sourceprimuv');

    return newprim;

}