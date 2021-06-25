#include "Assets/Materials/Specluar.hlsl"

//#define TOON_COLOR

//Timing functions
float easeOutElastic(float x){
    float c4 = (2 * 3.14) / 3;

    return x == 0
    ? 0
    : x == 1
    ? 1
    : pow(2, -10 * x) * sin((x * 10 - 0.75) * c4) + 1;
}

float easeOutBack(float x){
    float c1 = 1.70158;
    float c3 = c1 + 1;

    return 1 + c3 * pow(x - 1, 3) + c1 * pow(x - 1, 2);
}

float easeOutCubic(float x) {
    return 1 - pow(1 - x, 3);
}

float easeInOutCubic(float x) {
    return x < 0.5 ? 4 * x * x * x : 1 - pow(-2 * x + 2, 3) / 2;
}
//End Timing functions

//Operations

float opUnion( float d1, float d2 ) { return min(d1,d2); }

float opSubtraction( float d1, float d2 ) { return max(-d1,d2); }

float opIntersection( float d1, float d2 ) { return max(d1,d2); }

float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return lerp( d2, d1, h ) - k*h*(1.0-h); }

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return lerp( d2, -d1, h ) + k*h*(1.0-h); }

float opSmoothIntersection( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return lerp( d2, d1, h ) + k*h*(1.0-h); }

//End Operations

/*// Repeat around the origin by a fixed angle.
// For easier use, num of repetitions is use to specify the angle.
float pModPolar(inout float2 p, float repetitions) {
	float angle = 2*PI/repetitions;
	float a = atan(p.y / p.x) + angle/2.;
	float r = length(p);
	float c = floor(a/angle);
	a = (a % angle) - angle/2.;
	p = float2(cos(a), sin(a))*r;
	// For an odd number of repetitions, fix cell index of the cell in -x direction
	// (cell index would be e.g. -5 and 5 in the two halves of the cell):
	if (abs(c) >= (repetitions/2)) c = abs(c);
	return c;
}*/

float mod(float x, float y){
    return x - y * floor(x/y);
}

float2 pModPolar(inout float2 p, float repetitions, float fix) {
	float angle = 2*PI/repetitions;
	float a = atan(p.y / p.x) + angle/2.;
	float r = length(p);
	float c = floor(a/angle);
	a = mod(a,angle) - (angle/2.)*fix;
	p = float2(cos(a), sin(a))*r;

	return p;
}

// Same, but mirror every second cell so all boundaries match
float2 pModMirror2(inout float2 p, float2 size) {
	float2 halfsize = size*0.5;
	float2 c = floor((p + halfsize)/size);
	p = mod(p + halfsize, size) - halfsize;
	p *= mod(c,float2(2,2))*2 - float2(1,1);
	return c;
}

float3 opTwist( in float3 p, float k )
{
    float c = cos(k*p.y);
    float s = sin(k*p.y);
    float2x2 m = {c,-s,
                  s,c};
    float3 q = float3(mul(m, float2(p.x, p.z)),p.y);
    return q;
}


/*float opTwist( in sdf3d primitive, in vec3 p )
{
    const float k = 10.0; // or some other amount
    float c = cos(k*p.y);
    float s = sin(k*p.y);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(m*p.xz,p.y);
    return primitive(q);
}*/

float sdCappedCone(float3 p, float3 a, float3 b, float ra, float rb)
{
    float rba  = rb-ra;
    float baba = dot(b-a,b-a);
    float papa = dot(p-a,p-a);
    float paba = dot(p-a,b-a)/baba;
    float x = sqrt( papa - paba*paba*baba );
    float cax = max(0.0,x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;
    float k = rba*rba + baba;
    float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );
    float cbx = x-ra - f*rba;
    float cby = paba - f;
    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(cax*cax + cay*cay*baba,
                       cbx*cbx + cby*cby*baba) );
}

float sdCappedCylinder(float3 p, float3 a, float3 b, float r)
{
    float3 ba = b - a;
    float3 pa = p - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);
    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
}

float sphere_SDF(float3 p, float3 center, float radius){
    return length(p - center) - radius;
}

float box_SDF( float3 p, float3 b )
{
    float3 q = abs(p) - b;
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float boxOffset_SDF(float3 p, float3 center, float3 bounds){
    return box_SDF(p - center, bounds);
}

float sdCapsule_changed( float3 p, float3 a, float3 b, float r )
{
    float3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    float t = -((h * 2)-1)*((h * 2)-1) + 1;
    return length( pa - ba*h ) - r * t;
}

float daggerCapsule_changed( float3 p, float3 a, float3 b, float r )
{
    float3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.3 );
    float t = -h*h+1;
    return length( pa - ba*h ) - (r * t);
}

float grestnerWave(float scaledD, float offset, float amp1, float amp2){
    float scaledDir1 = scaledD;
    float sectionOffset = floor(scaledDir1) % 2;
    float offset1 = offset; 
    float sections = floor(scaledDir1 - (sectionOffset * offset1)) % 2;

    float norm1 = 1 + offset1;
    float2 tSections = float2(1,1);
    if(sections == 0){
        tSections = float2((scaledDir1 % 2) / (norm1), amp1);
    }else{
        tSections = float2((scaledDir1 % 2 - norm1) * (1 / (1 - offset1)), amp2);
    }
    float disp = sin(tSections.x * 3.14) * tSections.y;
    return disp;
}

float plane_SDF(in float3 p, float3 normal, float h){
    return dot(p,normal) + h;
}

float liquid_SDF(in float3 p, in float time, in UnityTexture2D noiseTex){
    float3 dir1 = normalize(float3(1,0.1,1));
    float alongDir1 = dot(dir1, p);
    float scale1 = 0.3;
    float scaledDir1 = alongDir1 * scale1 + time * 0.3;

    float3 dir2 = normalize(float3(0.5, 0, 1.0));
    float alongDir2 = dot(dir2, p);
    float scale2 = 0.5;
    float scaledDir2 = alongDir2 * scale2 + time * 2 + 0.6;

    float3 dir3 = normalize(float3(2.0, 0, 0.2));
    float alongDir3 = dot(dir3, p);
    float scale3 = 0.4;
    float scaledDir3 = alongDir3 * scale3 + time + 2;

    float disp = grestnerWave(scaledDir1, 0.3, 0.1, -0.2) + grestnerWave(scaledDir2, 0.8, 0.7, -0.08) + grestnerWave(scaledDir3, 0.5, 0.3, -0.2);
    float noiseVal = tex2D(noiseTex, float2(p.x, p.z) * 0.02);
    float3 pDisp = p + float3(0, disp, 0);
    return dot(pDisp,float3(0,1,0)) + 0;// boxOffset_SDF(pDisp, float3(0,0,0), float3(20,1,20)) - 0.8;
}

float3 liquid_normal_SDF(in float3 p, in float time, in UnityTexture2D noiseTex){
    const float h = 0.1; // replace by an appropriate value
    const float2 k = float2(1,-1);
    return normalize( k.xyy*liquid_SDF( p + k.xyy*h, time, noiseTex) + 
                      k.yyx*liquid_SDF( p + k.yyx*h, time, noiseTex ) + 
                      k.yxy*liquid_SDF( p + k.yxy*h, time, noiseTex ) + 
                      k.xxx*liquid_SDF( p + k.xxx*h, time, noiseTex ) );
}

float flaskTopOpening(in float3 p,in float3 flaskA,in float3 flaskB){
    float3 dirFlask = normalize(flaskB - flaskA);
    float3 A = flaskB - dirFlask * 0.1;
    float3 B = flaskB + dirFlask * 0.1;
    float3 StartToPoint = p - A;
    float3 AB = B - A;
    float unclampedT = dot(StartToPoint, AB) / (length(AB) * length(AB));
    float t = min(max(unclampedT, 0), 1);
    float3 PointOnLine = t * AB;
    return length(StartToPoint - PointOnLine) - 0.6;
}

float flaskMainBody_SDF(in float3 p, float3 A, float3 B, float radMin, float radMax){
    //This is a cylinder that is only capped at the top, and uses a smooth step function to change the radius
    float3 StartToPoint = p - A;
    float3 AB = B - A;
    float unclampedT = dot(StartToPoint, AB) / (length(AB) * length(AB));
    float t = min(max(unclampedT, 0), 1);
    float3 PointOnLine = t * AB;
    float oneMinusT = 1-t;
    float rad = smoothstep(0.3, 1, oneMinusT) * (radMax - radMin) + radMin; // mid rad is 0.3 and max is 0.3 + 0.8 // oneMinusT*oneMinusT * (radMax - radMin) + radMin;

    float3 topNormal = normalize(AB);
    float3 BP = p-B;
    float3 pPlane = BP - dot(BP, topNormal) * topNormal;
    float3 clamped_pPlane = normalize(pPlane) * min(rad, length(pPlane));

    float stepT = step(1, unclampedT); //This is used to change between the normal capsuled and the cap at the top
    float flaskForm = length(StartToPoint - PointOnLine) - rad;
    float flaskTop = length(BP - clamped_pPlane);

    return lerp(flaskForm, flaskTop, stepT);
}

float flask_SDF(in float3 p){
    float3 A = float3(0, -1, 0);
    float3 B = float3(0, 1.9, 0);

    float flaskMain = flaskMainBody_SDF(p, A, B, 0.36, 1.2);
    float flaskInner = flaskMainBody_SDF(p, A, B, 0.1, 1);

    //float flaskBody = opSmoothSubtraction(flaskInner, flaskMain, 0.1);

    float topThing = sdCappedCylinder(p, float3(0,1.65,0), float3(0,2,0), 0.55) - 0.05;
    float topHole = sdCappedCylinder(p, float3(0,1,0), float3(0,2.5,0), 0.35);
    //topHole = opSmoothUnion(flaskInner, topHole, 0.5);

    return opSmoothSubtraction(topHole, opSmoothUnion(flaskMain - 0.1, topThing, 0.1), 0.1);
}

float3 flask_normal_SDF(in float3 p){
    const float h = 0.01; // replace by an appropriate value
    const float2 k = float2(1,-1);
    return normalize( k.xyy*flask_SDF( p + k.xyy*h ) + 
                      k.yyx*flask_SDF( p + k.yyx*h ) + 
                      k.yxy*flask_SDF( p + k.yxy*h ) + 
                      k.xxx*flask_SDF( p + k.xxx*h ) );
}

float3 rotateAroundY(float3 p, float angle){
    float c = cos(angle);
    float s = sin(angle);
    float2x2 m = {c,-s,
                  s,c};
    float2 rotated = mul(m, float2(p.x, p.z));
    return float3(rotated.x,p.y,rotated.y);
}

float innerFlask_SDF(in float3 p, in float time){
    float3 A = float3(0, 1 * abs(sin(time)) - 1, 0);
    float3 B = float3(0, 3 * abs(sin(time)), 0);

    float3 pCopy = p;
    float capsule = sdCapsule_changed(pCopy, A, B, 0.3);
    float lowerCapsule = opSmoothSubtraction(sphere_SDF(pCopy, float3(0,-2,0), 0.8), capsule, 0.1);
    
    float dist = sphere_SDF(p, float3(0,-1,0), 0.9);
    float mult = 8;
    float speed = 6;
    float displacement = sin(mult*p.x)*sin(mult*p.y - time * speed)*sin(mult*p.z);
    float liquid = (dist * 0.7 + displacement * 0.2) * 0.5;
    return opSmoothUnion(lowerCapsule + displacement * 0.02, liquid, 0.8);//sdCapsule(q, A, B, 0.15) * 0.05; //opSmoothIntersection(innerFlask, plane, 0.4);
}

float dagger_SDF(float3 p, out int matSelect){
    matSelect = 0;

    float3 pMirror = p;
    pMirror.x = abs(pMirror.x) + 0.01;
    float3 wavePos = float3(pMirror.x, pMirror.y, pMirror.z + sin((pMirror.y - 4) * 3) * 0.2);
    float daggerCapsule = daggerCapsule_changed(float3(wavePos.x*1.6, wavePos.y, wavePos.z*0.5), float3(0, 4, 0), float3(0,10, 0), 0.8);
    float middelCut = daggerCapsule_changed(wavePos * float3(1,1,0.5), float3(0.55, 2, 0), float3(0.55, 12, 0), 0.4);
    float plane = plane_SDF(pMirror, float3(-1,0,0), 0);
    float plane2 = plane_SDF(pMirror, float3(-1,0,0), 0.3);
    float plane3 = plane_SDF(pMirror, float3(0,-1,0), 3.5);

    float handleDetail = sphere_SDF(float3(p.x, p.y, abs(p.z)), float3(0, 3.5, 3), 0.7);
    float handleDetailSub = sdCapsule_changed(float3(p.x, p.y, abs(p.z)), float3(2, 3.5, 3), float3(-2, 3.5, 3), 0.4);
    handleDetail = opSmoothSubtraction(handleDetailSub, handleDetail, 0.4);
    float handle = opSmoothUnion(handleDetail, sdCapsule_changed(p, float3(0,3.5, -3), float3(0,3.5, 3), 0.5), 0.2);
    float blade = opSmoothIntersection(plane3, opSmoothSubtraction(middelCut, opSmoothSubtraction(plane2, opSmoothIntersection(daggerCapsule, plane, 0.01), 0.01), 0.001), 0.05);

    if(handle < blade){
        matSelect = 1;
    }

    return opUnion(handle, blade);
}

float liquidAndDagger_SDF(float3 p, float time, float animTime, out int matSelect, out float testColor){
    matSelect = 0;
    float createLiquidDaggerTime = clamp(animTime, 0, 80)/80;
    createLiquidDaggerTime = easeInOutCubic(createLiquidDaggerTime);
    float createSwordBodyTime = (clamp(animTime, 30, 100) - 30)/70;
    createSwordBodyTime = easeInOutCubic(createSwordBodyTime);

    float displacementTime = (clamp(animTime, 40, 80) - 40)/40;

    float3 A = float3(0, -2 + 3 * createLiquidDaggerTime , 0);
    float3 B = float3(0, 1 + 9 * createLiquidDaggerTime, 0);

    float3 daggerCapsuleA = float3(0, -0.5 * createLiquidDaggerTime , 0);
    float3 daggerCapsuleCreateA = float3(0, -1 , 0);
    float3 daggerCapsuleCreateB = float3(0, -0.5 + 12 * createSwordBodyTime , 0);

    float mult = 6 - 3 * createLiquidDaggerTime;
    float speed = 6;
    float displacement = sin(mult*p.x)*sin(mult*p.y + time * speed)*sin(mult*p.z);
    
    float dist = sphere_SDF(p, float3(0,-1 + 2.4 * createLiquidDaggerTime,0), 1 * (1 - createLiquidDaggerTime));
    float capsule = sdCapsule_changed(p, A, B, 0.3 + 0.2 * createLiquidDaggerTime);
    float swordCapsule = sdCapsule_changed(p, daggerCapsuleA, B, 0.4 + 4.5 * createLiquidDaggerTime);
    
    float daggerCapsuleCreate = sdCapsule_changed(p, daggerCapsuleCreateA, daggerCapsuleCreateB, 0 + 5 * createSwordBodyTime);

    float liquid = opSmoothUnion(dist, capsule, 0.4);

    float dagger = dagger_SDF(p, matSelect);

    float intersect = opSmoothIntersection(dagger, swordCapsule, 0.5);
    float allLiquid = opSmoothUnion(intersect, liquid, 0.4);
    allLiquid = opSmoothSubtraction(daggerCapsuleCreate, allLiquid, 1.3);

    float cretingDagger = opSmoothIntersection(daggerCapsuleCreate, dagger, 1);
    testColor = smoothstep(0, 0.4, saturate(exp(-allLiquid) - 0.3));
    float liquidDisplaced = (allLiquid + displacement * (0.1 + 0.2 * displacementTime));
    

    return opSmoothUnion(liquidDisplaced, cretingDagger, 1.5)  * 0.5;//opUnion(daggerCapsuleCreate, (allLiquid + displacement * 0.2)) * 0.5;
}

float3 liquidAndDagger_normal_SDF(float3 p, float time, float animTime){
    const float h = 0.01; // replace by an appropriate value
    const float2 k = float2(1,-1);
    int matSelect;
    float testColor;
    return normalize( k.xyy*liquidAndDagger_SDF( p + k.xyy*h, time, animTime, matSelect, testColor) + 
                      k.yyx*liquidAndDagger_SDF( p + k.yyx*h, time, animTime, matSelect, testColor) + 
                      k.yxy*liquidAndDagger_SDF( p + k.yxy*h, time, animTime, matSelect, testColor) + 
                      k.xxx*liquidAndDagger_SDF( p + k.xxx*h, time, animTime, matSelect, testColor) );
}

float3 innerFlask_normal_SDF(in float3 p, in float time){
    const float h = 0.01; // replace by an appropriate value
    const float2 k = float2(1,-1);
    return normalize( k.xyy*innerFlask_SDF( p + k.xyy*h, time ) + 
                      k.yyx*innerFlask_SDF( p + k.yyx*h, time ) + 
                      k.yxy*innerFlask_SDF( p + k.yxy*h, time ) + 
                      k.xxx*innerFlask_SDF( p + k.xxx*h, time ) );
}

float flaskDetail(in float3 p, in float time){

    float3 pCopy = p;

    float3 pOffset = p + float3(0, 0, 0);
    float positionAngle = atan(pOffset.x / pOffset.z);
    pOffset.y += (sin((positionAngle + time) * 8) + sin((positionAngle) * -12)) * 0.1;
    float mainBase = sdCappedCone(pOffset, float3(0, -3.3, 0), float3(0, -2.5, 0), 5, 4);
    float subBase = sdCappedCone(pOffset, float3(0, -3.3, 0), float3(0, -2.5, 0), 4, 5);
    float intersectBase = sdCappedCone(pOffset, float3(0, -3.6, 0), float3(0, -2.5, 0), 1, 9);
    float base = opSmoothIntersection(opSmoothSubtraction(subBase, mainBase, 0.05), intersectBase, 0.02);

    pCopy.xz = pModPolar(pCopy.xz, 4, 1);
    
    float3 moveDirection = normalize(float3(1,1,0));
    float3 pFloating = (pCopy) + moveDirection * (sin(time * PI) * 0.5);

    float floatingSphere = sphere_SDF(pFloating, float3(1.5, 0.5, 0), 0.2);
    
    pCopy.y *= 0.75;

    float scale = lerp(0.59, 1.3, smoothstep(-0.5, 0.5, pCopy.y));
    pCopy.xz *= scale;
    float s1 = sphere_SDF(pCopy, float3(0.5,0,0), 0.45);
    float s2 = sphere_SDF(pCopy, float3(1,0,0), 0.3);
    float plane = plane_SDF(pCopy, float3(-1,0,0), 0.73);
    float3 pCopyOffset = pCopy + float3(0,0.5,0);
    float outerThing = sdCapsule_changed(pCopyOffset, float3(0.25,-1,0), float3(1.2,-1,0), 0.4);
    outerThing = opSmoothSubtraction(sphere_SDF(pCopyOffset, float3(0,-1,0),0.8), outerThing, 0.05);
    outerThing = opSmoothIntersection(plane_SDF(p, float3(0,-1,0), -2), outerThing, 0.3);
    
    float mainBody = opSmoothIntersection(opSmoothSubtraction(s2, s1, 0.05), plane, 0.1);
    return opSmoothUnion(base, opSmoothUnion(outerThing, opSmoothUnion(floatingSphere, mainBody, 0.6), 0.1), 0.001);
}

float3 flaskDetail_normal_SDF(in float3 p, in float time){
    const float h = 0.01; // replace by an appropriate value
    const float2 k = float2(1,-1);
    return normalize( k.xyy*flaskDetail( p + k.xyy*h, time ) + 
                      k.yyx*flaskDetail( p + k.yyx*h, time ) + 
                      k.yxy*flaskDetail( p + k.yxy*h, time ) + 
                      k.xxx*flaskDetail( p + k.xxx*h, time ) );
}

float tarDagger_SDF(in float3 p, in float time, out int matSelect){
    matSelect = 0;

    float3 pMirror = p;
    pMirror.x = abs(pMirror.x) + 0.1;
    float3 wavePos = float3(pMirror.x, pMirror.y, pMirror.z + sin((pMirror.y - 4) * 3) * 0.2);
    float daggerCapsule = daggerCapsule_changed(float3(wavePos.x*1.6, wavePos.y, wavePos.z*0.5), float3(0, 4, 0), float3(0,10, 0), 0.8);
    float middelCut = daggerCapsule_changed(wavePos * float3(1,1,0.5), float3(0.55, 2, 0), float3(0.55, 12, 0), 0.4);
    float plane = plane_SDF(pMirror, float3(-1,0,0), 0);
    float plane2 = plane_SDF(pMirror, float3(-1,0,0), 0.3);
    float plane3 = plane_SDF(pMirror, float3(0,-1,0), 3.5);

    float handleDetail = sphere_SDF(float3(p.x, p.y, abs(p.z)), float3(0, 3.5, 3), 0.7);
    float handleDetailSub = sdCapsule_changed(float3(p.x, p.y, abs(p.z)), float3(2, 3.5, 3), float3(-2, 3.5, 3), 0.4);
    handleDetail = opSmoothSubtraction(handleDetailSub, handleDetail, 0.4);
    float handle = opSmoothUnion(handleDetail, sdCapsule_changed(p, float3(0,3.5, -3), float3(0,3.5, 3), 0.5), 0.2);
    float blade = opSmoothIntersection(plane3, opSmoothSubtraction(middelCut, opSmoothSubtraction(plane2, opSmoothIntersection(daggerCapsule, plane, 0.01), 0.01), 0.001), 0.05);

    if(handle < blade){
        matSelect = 1;
    }

    float completeBlade = opUnion(handle, blade);
    float tarMult = 4;
    float tarSpeed = 10;
    float speed = 0.15;

    float upperLimit = abs(sin(time * speed)) * 11;
    float tarDisplacement = sin(tarMult*p.x)*sin(tarMult*p.y - time * tarSpeed)*sin(tarMult*p.z) * (0.3 * smoothstep( upperLimit + 1, upperLimit, p.y )-0.1);
    

    float tarMask = box_SDF(p + float3(0,- (abs(sin(time * speed)) * 10.7),0), float3(1, 1, 8));
    float normalMask = box_SDF(p + float3(0,-9 - (abs(sin(time * speed)) * 10.7),0), float3(2,8, 8));

    float tarDagger = opSmoothIntersection(tarMask, completeBlade + tarDisplacement, 0.8);
    float normalDagger = opIntersection(normalMask, completeBlade);

    if(tarDagger < normalDagger){
        matSelect = 2;
    }

    return opSmoothUnion(tarDagger, normalDagger, 1.3);
}

float3 tarDagger_normal_SDF(in float3 p, in float time){
    const float h = 0.01; // replace by an appropriate value
    const float2 k = float2(1,-1);
    int dumpMat;
    return normalize( k.xyy*tarDagger_SDF( p + k.xyy*h, time, dumpMat) + 
                      k.yyx*tarDagger_SDF( p + k.yyx*h, time, dumpMat) + 
                      k.yxy*tarDagger_SDF( p + k.yxy*h, time, dumpMat) + 
                      k.xxx*tarDagger_SDF( p + k.xxx*h, time, dumpMat) );
}

void PotionSwordSDF_float(float3 startPosition, float3 directionNorm, float3 positionOffset, float time, UnityTexture2D noiseTex, float animTime, out float alpha, out float3 color){
    float minDist = 0.01;
    //startPosition -= positionOffset;
    alpha = 0;
    color = float3(0,0,0);

    float distanceTraveled = 0;
    float outerLightData = 0;
    float outerFlaskDepth = 100000;

    float flaskLightR = 0;
    float flaskLightG = 0;
    float flaskLightB = 0;
    
    float flaskSpecularR = 0;
    float flaskSpecularG = 0;
    float flaskSpecularB = 0;
    
    [loop]
    for(int i = 0; i < 50; i++){
        float3 testPoint = startPosition + directionNorm * distanceTraveled;
        float dist = flask_SDF(testPoint);
        if(abs(dist) < minDist){
            alpha = 1;
            float3 normal = flask_normal_SDF(testPoint);
            float vdotn = dot(normal, -directionNorm);
            float3 lightDir = normalize(float3(1,1,1));
            
            float s = specular(lightDir, normal, directionNorm, 8, 1, float3(1,2,1));
            flaskSpecularB = s;
            flaskSpecularR = specular(lightDir, normal, directionNorm, 7, 1, float3(1,2,0.9)); 
            flaskSpecularG = specular(lightDir, normal, directionNorm, 6, 1, float3(1,2,1.1));
            //s = smoothstep(0, 0.6, s);

            flaskLightB = exp(-vdotn * 3.5);
            flaskLightR = exp(-vdotn * 4.5);
            flaskLightG = exp(-vdotn * 5);

            #ifdef TOON_COLOR
                flaskSpecularB = smoothstep(0.6, 0.6, flaskSpecularB);
                flaskSpecularR = smoothstep(0.6, 0.6, flaskSpecularR);
                flaskSpecularG = smoothstep(0.6, 0.6, flaskSpecularG);

                flaskLightB = smoothstep(0.22,0.22, flaskLightB);
                flaskLightG = smoothstep(0.25,0.25, flaskLightG);
                flaskLightR = smoothstep(0.3,0.3, flaskLightR);
            #else
                flaskSpecularB = smoothstep(0.3, 0.6, flaskSpecularB);
                flaskSpecularR = smoothstep(0.3, 0.6, flaskSpecularR);
                flaskSpecularG = smoothstep(0.3, 0.6, flaskSpecularG);

                flaskLightB = smoothstep(0,0.5, flaskLightB);
                flaskLightG = smoothstep(0,0.6, flaskLightG);
                flaskLightR = smoothstep(0,0.7, flaskLightR);
            #endif

            //light = smoothstep(0,0.3, light);

            outerLightData = saturate(( flaskLightB + flaskLightR + flaskLightG + flaskSpecularB + flaskSpecularR + flaskSpecularG)) * 0.6;
            outerFlaskDepth = distanceTraveled;
            //color = lerp(float3(0,0.1,0.3),float3(0.3,0.9,1), (exp(-vdotn * 8) + s));
            break;
        }
        distanceTraveled += dist * 0.7;
    }

    alpha = outerLightData;

    //Liquid
    bool liquidRendered = false;
    float liquidDepth = 1000000;
    distanceTraveled = 0;
    [loop]
    for(int i = 0; i < 70; i++){
        float3 testPoint = startPosition + directionNorm * distanceTraveled;
        int matSelect;
        float testColor;
        float dist = liquidAndDagger_SDF(testPoint, time, animTime, matSelect, testColor);
        if(abs(dist) < minDist){
            alpha = 1;
            liquidRendered = true;
            float3 normal = liquidAndDagger_normal_SDF(testPoint, time, animTime);
            float vdotn = dot(normal, -directionNorm);
            float3 lightDir = normalize(float3(1,1,1));
            float s = specular(lightDir, normal, directionNorm, 6, 1, float3(1,2,1));
            float fresnel = exp(-vdotn * 8);
            fresnel = clamp(fresnel, -1, 4);
            
            float sMetal = specular(lightDir, normal, directionNorm, 10, 2, float3(1,2,1));
            float fresnelMetal = (exp(-vdotn * 6));
            fresnelMetal = clamp(fresnelMetal, -1, 1.5);

            float lightData = (fresnel + s);
            float lightDataMetal = fresnelMetal + sMetal;

            #ifdef TOON_COLOR
                lightData = smoothstep(0.1, 0.1, lightData);
            #endif

            if(matSelect == 0){
                color = lerp(float3(0.1,0.1,0.1), float3(0.6, 0.6, 0.6), lightDataMetal);
            }else if(matSelect == 1){
                color = lerp(float3(0.5,0.1,0), float3(1, 0.5, 0), lightDataMetal);
            }

            float3 tarColor = lerp(float3(0,0,0),float3(0,0.5,0.5), lightData);

            color = lerp(color, tarColor, testColor);

            //color = lerp(float3(0,0,0),float3(0,0.5,0.5), lightData);
            liquidDepth = distanceTraveled;
            break;
        }
        distanceTraveled += dist;
    }

    if(outerFlaskDepth < liquidDepth){
        float3 blueColor = float3(0.5,0.5,1) * (flaskLightB + flaskSpecularB);
        float3 redColor = float3(1, 0.2, 0.5) * (flaskLightR + flaskSpecularR);
        float3 greenColor = float3(0.2, 1, 0.5) * (flaskLightG + flaskSpecularG);
        color = lerp(color, saturate(blueColor + redColor + greenColor), outerLightData);
        //color = lerp(color, float3(flaskLightR + flaskSpecularR, flaskLightG + flaskSpecularG, flaskLightB + flaskSpecularB), outerLightData);
    }

    distanceTraveled = 0;
    float detailDepth = 100000;
    [loop]
    for(int i = 0; i < 120; i++){
        float3 testPoint = startPosition + directionNorm * distanceTraveled;
        float dist = flaskDetail(testPoint, time);
        if(abs(dist) < minDist){
            detailDepth = distanceTraveled;
            float3 prevColor = color;
            float3 normal = flaskDetail_normal_SDF(testPoint, time);
            float vdotn = dot(normal, -directionNorm);
            float3 lightDir = normalize(float3(1,1,1));
            float s = specular(lightDir, normal, directionNorm, 10, 4, float3(1,2,1));
            float fresnel = (exp(-vdotn * 6));
            fresnel = clamp(fresnel, -1, 4);

            #ifdef TOON_COLOR
                s = smoothstep(0.7, 1.2, s) * 3;
                fresnel = smoothstep(0.2, 0.21, fresnel);
            #endif
            float lightData = fresnel + s;
            //lightData = smoothstep(0, 0.1, lightData);
            if(liquidDepth > detailDepth){
                color = lerp(float3(0.5,0.1,0), float3(1, 0.5, 0), lightData);
                alpha = 1;
                //color = lerp(color, float3(0, 1, 0.5), lightData + 0.1);//lerp(float3(0,0.1,0.5), float3(0, 1, 0.5), lightData);
                if(outerFlaskDepth < distanceTraveled){
                    color = lerp(color, prevColor, saturate(outerLightData));
                }else{
                    color = color;
                }
            }
            break;
        }
        distanceTraveled += dist * 0.5;
    }

    /*distanceTraveled = 0; 
    [loop]
    for(int i = 0; i < 150; i++){
        float3 testPoint = startPosition + directionNorm * distanceTraveled;
        float dist = liquid_SDF(testPoint, time + 100, noiseTex);
        if(abs(dist) < minDist && distanceTraveled < liquidDepth && distanceTraveled < detailDepth && distanceTraveled < outerFlaskDepth){
            alpha = 1;
            float3 normal = liquid_normal_SDF(testPoint, time + 100, noiseTex);
            float vdotn = dot(normal, -directionNorm);
            float3 lightDir = normalize(float3(1,1,1));
            float s = specular(lightDir, normal, directionNorm, 6, 2);
            color = lerp(float3(0,0,0.5),float3(0,0.3,2), s);
            break;
        }
        distanceTraveled += dist;
    }*/
    /*
    distanceTraveled = 0;
    float daggerDepth = 100000;
    [loop]
    for(int i = 0; i < 120; i++){
        float3 testPoint = startPosition + directionNorm * distanceTraveled;
        int matSelector;
        float dist = tarDagger_SDF(testPoint, time, matSelector);
        if(abs(dist) < minDist){
            daggerDepth = distanceTraveled;
            float3 prevColor = color;
            float3 normal = tarDagger_normal_SDF(testPoint, time);
            float vdotn = dot(normal, -directionNorm);
            float3 lightDir = normalize(float3(1,1,1));
            float s = specular(lightDir, normal, directionNorm, 10, 2, float3(1,2,1));
            float fresnel = (exp(-vdotn * 6));
            fresnel = clamp(fresnel, -1, 1);
            #ifdef TOON_COLOR
                s = smoothstep(0.5, 1.6, s) *2.5;
                fresnel = smoothstep(0.2, 0.21, fresnel);
            #endif
            float lightData = fresnel + s;
            //lightData = smoothstep(0, 0.1, lightData);
            
            if(liquidDepth > daggerDepth && detailDepth > daggerDepth){
                alpha = 1;
                if(matSelector == 0){
                    color = lerp(float3(0.1,0.1,0.1), float3(0.6, 0.6, 0.6), lightData);
                }else if(matSelector == 1){
                    color = lerp(float3(0.5,0.1,0), float3(1, 0.5, 0), lightData);
                }else{
                    color = lerp(float3(0,0,0),float3(0,0.5,0.5), lightData);
                }

                if(outerFlaskDepth < distanceTraveled){
                    color = lerp(color, prevColor, saturate(outerLightData));
                }else{
                    color = color;
                }
                //color = lerp(float3(0.5,0.1,0), float3(1, 0.5, 0), lightData);
                //color = lerp(color, float3(0, 1, 0.5), lightData + 0.1);//lerp(float3(0,0.1,0.5), float3(0, 1, 0.5), lightData);
                /*if(outerFlaskDepth < distanceTraveled){
                    alpha = 1;
                    color = lerp(color, prevColor, saturate(outerLightData));
                }else{
                    alpha = 1;
                    color = color;
                }
            }
            break;
        }
        distanceTraveled += dist * 0.5;
    }*/

    /*distanceTraveled = 0; 
    [loop]
    for(int i = 0; i < 150; i++){
        float3 testPoint = startPosition + directionNorm * distanceTraveled;
        float dist = liquid_SDF(testPoint, time + 100, noiseTex);
        if(abs(dist) < minDist && distanceTraveled < liquidDepth && distanceTraveled < detailDepth && distanceTraveled < outerFlaskDepth){
            alpha = 1;
            float3 normal = liquid_normal_SDF(testPoint, time + 100, noiseTex);
            float vdotn = dot(normal, -directionNorm);
            float3 lightDir = normalize(float3(1,1,1));
            float s = specular(lightDir, normal, directionNorm, 6, 2);
            color = lerp(float3(0,0,0.5),float3(0,0.3,2), s);
            break;
        }
        distanceTraveled += dist;
    }*/
}