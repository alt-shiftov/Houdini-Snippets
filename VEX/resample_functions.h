// ### Resample by length - resamplebylength()
// Функция для ресампла линии. с pin points
// Работает не совсем корректно, если много линий, то может хуйню выдавать. 

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

// VECTOR2 UVW

function void resamplebylength(int geo; vector2 uvw; int prim; float minlen; int newprim; int newpts[]; float newuvs[]){
    resamplebylength(geo, uvw[0], uvw[1], prim, minlen, newprim, newpts, newuvs);
}

// RESAMPLE WITH PIN POINTS
function void resamplebylength(int geo; float uv01, uv02; int prim; float minlen; int newprim; int pinpts[]; float pinuvs[]; int newpts[]; float newuvs[]){
    // before pin points
    resamplebylength(geo, uv01, pinuvs[0], prim, minlen, newprim, newpts, newuvs);

    // between pin points
    for(int i = 0; i < len(pinpts)-1; i++){
        int nexti = (i+1) % len(pinpts);
        resamplebylength(geo, pinuvs[i], pinuvs[nexti], prim, minlen, newprim, newpts, newuvs);
    }
    // after pin points
    resamplebylength(geo, pinuvs[-1], uv02, prim, minlen, newprim, newpts, newuvs);
}

// RETURN NEW PRIM (WITH PIN POINTS)
function int resamplebylength(int geo; float uv01, uv02; int prim; float minlen; int pinpts[]; float pinuvs[]; int newpts[]; float newuvs[]){
    int newprim = addprim(geo, 'polyline');
    resamplebylength(geo, uv01, uv02, prim, minlen, newprim, pinpts, pinuvs, newpts, newuvs);
    return newprim;
}

function int resamplebylength(int geo; vector2 uvw; int prim; float minlen; int pinpts[]; float pinuvs[]; int newpts[]; float newuvs[]){
    int newprim = addprim(geo, 'polyline');
    resamplebylength(geo, uvw[0], uvw[1], prim, minlen, newprim, pinpts, pinuvs, newpts, newuvs);
    return newprim;
}


function void setuvtopoints(int geo; int pts[]; int primnum; float uvs[]; string primchannel; string uvchannel){
    foreach(int idx; int pt; pts){
        float u = uvs[idx];
        setpointattrib(geo, primchannel, pt, primnum);
        setpointattrib(geo, uvchannel, pt, set(u, 0, 0));
    }
}


// RETURN NEW PRIM
function int resamplebylength(int geo; float uv01, uv02; int prim; float minlen; int newpts[]; float newuvs[]){
    int newprim = addprim(geo, 'polyline');
    resamplebylength(geo, uv01, uv02, prim, minlen, newprim, newpts, newuvs);
    return newprim;
}

function int resamplebylength(int geo; vector2 uvw; int prim; float minlen; int newpts[]; float newuvs[]){
    int newprim = addprim(geo, 'polyline');
    resamplebylength(geo, uvw[0], uvw[1], prim, minlen, newprim, newpts, newuvs);
    return newprim;
}

// SIMPLE VERSION OF RESAMPLE
function int resamplebylength(int geo; int prim; float minlen){
    int newprim = addprim(geo, 'polyline');
    int newpts[] = array();
    float newuvs[] = array();
    resamplebylength(geo, 0.0, 1.0, prim, minlen, newprim, newpts, newuvs);
    setuvtopoints(geo, newpts, prim, newuvs, 'sourceprim', 'sourceprimuv');
    return newprim;
}