#version 430

layout(location = 1) uniform int pass;
layout(location = 3) uniform int mode = 0;
layout(location = 6) uniform float time;
layout(location = 7) uniform vec4 slider;
layout(location = 8) uniform int scene = 0;
layout (location = 10) uniform mat4 rotation;
layout (location = 11) uniform vec3 camPos;
layout (location = 12) uniform vec3 colSlider;
layout (location = 13) uniform int switchMtrl;
layout (location = 14) uniform int switchTerm;
layout (location = 15) uniform float mslider;
layout (location = 16) uniform bool stopMove;

layout(binding = 0) uniform sampler2D backfaces_tex;

layout(location = 0) out vec4 fragcolor;  
         
in vec3 vpos;  

//forward function declarations
vec4 raytracedcolor(vec3 rayStart, vec3 rayStop);
vec4 clear_color(vec3 rayDir);
vec4 lighting(vec3 pos, vec3 rayDir);
float distToShape(vec3 pos);
vec3 normal(vec3 pos);

const vec3 light_pos = vec3(0.8, 0.8, 0.8);

const vec4 La = vec4(0.75, 0.75, 0.75, 1.0);
const vec4 Ld = vec4(0.74, 0.74, 0.74, 1.0);
const vec4 Ls = vec4(1.0, 1.0, 0.74, 1.0);

const vec4 Ka = vec4(0.4, 0.4, 0.34, 1.0);
const vec4 Kd = vec4(1.0, 1.0, 0.73, 1.0);
const vec4 Ks = vec4(0.1, 0.1, 0.073, 1.0);

//vec2 sceneRes = vec2 (1.0f);
//given distance function==============================================================

//primitive combination--------------------------------------------------
float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); }

float opSmoothIntersection( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h); }

//primitive alteration-------------------------------------------------------
vec4 opElongate( in vec3 p, in vec3 h )
{
    vec3 q = abs(p)-h;
    return vec4( max(q,0.0), min(max(q.x,max(q.y,q.z)),0.0) );
}

float opOnion( in float sdf, in float thickness )
{
    return abs(sdf)-thickness;
}

vec3 opTwist( vec3 p, float t)
{
    float  c = cos(10.0*slider.x*p.y+10.0f*t);
    float  s = sin(10.0*slider.x*p.y+10.0f*t);
    mat2   m = mat2(c,-s,s,c);
    return vec3(m*p.xz,p.y);
}


// shape function definitions  ----------------------------------------------------------------
float sdFloor(vec3 p){
	return p.z; 
}
float sdCone( vec3 p, vec2 c )
{
    // c must be normalized
    float q = length(p.xy);
    return dot(c,vec2(q,p.z));
}
float sdHexPrism( vec3 p, vec2 h )
{
    const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec2 d = vec2(
       length(p.xy-vec2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
       p.z-h.y );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h ) - r;
}

float sdOctahedron( in vec3 p, in float s)
{
    p = abs(p);
    float m = p.x+p.y+p.z-s;
    vec3 q;
         if( 3.0*p.x < m ) q = p.xyz;
    else if( 3.0*p.y < m ) q = p.yzx;
    else if( 3.0*p.z < m ) q = p.zxy;
    else return m*0.57735027;
    
    float k = clamp(0.5*(q.z-q.y+s),0.0,s); 
    return length(vec3(q.x,q.y-s+k,q.z-k)); 
}

float sdEllipsoid( in vec3 p, in vec3 r )
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float sdSphere( vec3 p, float s )
{
	return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

//start here===========================================================================================
//compare result of sceneSDF, declare as global var to hold ID
vec2 res = vec2(0.0);

vec2 opU( vec2 dist1, vec2 dist2 )
{
	return (dist1.x<dist2.x) ? dist1 : dist2;
}

//scene dist function
float sceneSDF(vec3 pos)
{
	if(scene == 0)
	{
		const float radius = 0.4;
		vec3 offset = 2.0*slider.xyz;

		float d3 = sdCapsule(vec3(inverse(rotation)*vec4(pos, 1.0f)), vec3(0.01f+ 0.05f *abs(cos(time*0.8))), vec3(0.05f+0.04f *abs( sin(time*0.9))), 0.03f+0.5*abs(cos(time*0.01)));
		//float float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
		float d2 = sdOctahedron(opTwist(vec3(inverse(rotation)*vec4(pos, 1.0f)), 0.0f) - vec3(0.9* cos(time), 0.9 * sin(time), 0.0f), radius);

		float d4 = sdOctahedron(opTwist(vec3(inverse(rotation)*vec4(pos, 1.0f)), 0.0f) - vec3(0.9* sin(time-1.5), 0.9 * -cos(time- 1.5), 0.0f), radius);

		return min(min(d2, d3), d4);

	}

	else if(scene == 1)
	{
		//elongated torus
		vec3 q = pos- vec3(-0.5f,0.0f,0.0f);
		vec4 w = opElongate(q, vec3(0.0f, slider.z, slider.w));
		float d1 = sdTorus(w.xyz, vec2 (0.3f, 0.1f));
		float d2 = sdTorus(opTwist((pos-vec3(0.5f, 0.0f, 0.0f)), time*slider.y), vec2(0.3f, 0.1f));
		return min(d1,d2);
	}

	else if(scene == 2)
	{
		const float radius = 0.4;
		vec3 c = vec3(1.0, 1.0, 1.0) + 2.0*slider.xyz;
		vec3 q = mod(pos,c)-0.5*c;
		float d1 = sdSphere(q, radius);
		float d2 = sdCapsule( q, vec3 (0.0), vec3 (0.5),  radius );
		float d3 =  sdTorus(  q, vec2 (0.3, 0.1) );

		float d3On = opOnion( d3, slider.w);

		float d4 =  sdHexPrism(  q, vec2 (radius, radius/2.0f) );
		
		return opSmoothIntersection( d3On, d4,  0.3 ); 
	}
	else if(scene == 3)
	{
		float dMove = sdEllipsoid(vec3(inverse(rotation)*vec4(pos, 1.0f)), vec3(0.1 +0.05*abs(cos(time*0.6)), 0.2 + 0.04*abs(sin(time*0.9)), 0.4- 0.3*abs(cos(time*0.99))));
		float dStop = sdEllipsoid(pos, vec3(0.3, 0.4, 0.3));
		float d3 = 0.0f;
		if (!stopMove) {d3 = dMove;}
		else {d3=dStop;}
		//float d3 = mix (dMove, dStop, stopMove);

		float d2 = sdOctahedron(opTwist(vec3(inverse(rotation)*vec4(pos, 1.0f)), 0.0f) - vec3(0.9* cos(time), 0.9 * sin(time), 0.0f), 0.4);
		float d4 = sdOctahedron(opTwist(vec3(inverse(rotation)*vec4(pos, 1.0f)), 0.0f) - vec3(0.9* sin(time-1.5), 0.9 * -cos(time- 1.5), 0.0f), 0.4);
		float d5 = sdFloor(pos-vec3(-1.0f));

		//compare shortest result and assign ID to each dist
		res = vec2(d5, 5.0);
		res = opU(res, vec2(d2, 2.0));
		res = opU(res, vec2(d3, 3.0));
		res = opU(res, vec2(d4, 4.0));
		
	return res.x;
	}
}



//normal vector of the shape we are drawing.
//Estimated as the gradient of the signed distance function.
vec3 normal(vec3 pos)
{
	const float h = 0.001;
	const vec3 Xh = vec3(h, 0.0, 0.0);	
	const vec3 Yh = vec3(0.0, h, 0.0);	
	const vec3 Zh = vec3(0.0, 0.0, h);	

	return normalize(vec3(sceneSDF(pos+Xh)-sceneSDF(pos-Xh), sceneSDF(pos+Yh)-sceneSDF(pos-Yh), sceneSDF(pos+Zh)-sceneSDF(pos-Zh)));
}

//shadow function (improved)
float softshadow( in vec3 pos, in vec3 lightDir, float mint, float maxt, float k )
{
    float res = 1.0;
	float t = mint;

	float h = sceneSDF(pos + lightDir*t);

    float ph = 1e20;
    for( float t=mint; t < maxt;  )
    {
        h = sceneSDF(pos + lightDir*t);
        if( h< 1e-4) return 0.0;
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min( res, k*d/max(0.0,t-y));
        ph = h;
        t += h;
		if( res<0.0001 || t>maxt ) break;
    }
    return clamp(res, 0.0f, 1.0f);
}
//**PBR FDG Equation======================================================================
const float PI = 3.14159265359;

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;
	
    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
	
    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float nom   = NdotV;
    float denom = NdotV * (1.0 - k) + k;
	
    return nom / denom;
}
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);
	 
    return ggx1 * ggx2;
}
vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);

} 

float DistributionBeckmann(vec3 N, vec3 H, float roughness){
 float a = acos(max(dot(N,H), 0.0));
 float m = roughness;

 float num = exp(-  ((tan(a)/m)* (tan(a)/m)) );
 float denom = 4.0 * m*m * pow( cos(a) ,4.0);

return num/denom;
}

float GeometryCookTorrance(vec3 N, vec3 V, vec3 L, vec3 H)
{
 float NdotH = max(dot(N, H), 0.0);
 float NdotV = max(dot(N, V), 0.0);
 float VdotH = max(dot(V, H), 0.0);
 float NdotL = max(dot(N, L), 0.0);

 return min( min(1, 2.0*NdotH*NdotV/VdotH), 2.0* NdotH*NdotL/VdotH);
}

//use fresnel term from class notEqual
vec3 fresnelCookTorrance(float nDotV, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - nDotV, 5.0);
} 

//render PBR
//preset albedo color
vec3 alumi = vec3(0.96, 0.96, 0.97);
vec3 gold = vec3(1.00, 0.86, 0.57);
vec3 copper = vec3(0.95, 0.64, 0.54);
vec3 plastic = vec3(0.24, 0.24, 0.24);
//material params
vec3 albedo;
vec3 albedoKa = vec3(1.0, 1.0, 1.0);
vec3 albedoKd = vec3(1.0, 1.0, 1.0);
vec3 albedoKs = vec3(1.0, 1.0, 1.0);
float roughness = 0.1; 
float roughnessD = 0.1; //roughness only for Distribution
float metallic = 0.3;  
vec3 F0 = vec3(0.04);


vec3 renderPBR(vec3 sampleP){
//light_pos: const vec3 light_pos = vec3(5.0, 5.0, 5.0);
//light color: const vec4 Ls = vec4(1.0, 1.0, 0.74, 1.0);

//assign param to each dist
if(res.y == 3.0) 
 {
       // roughness = 0.4;
		 roughnessD = 0.5-mslider;
        metallic= 0.9;
		if (switchMtrl == 0) albedoKa = colSlider;
		if (switchMtrl == 1) albedoKd = colSlider;
		if (switchMtrl == 2) albedoKs = colSlider;
		if (switchMtrl == 3) F0 = colSlider;
		
        //F0 =  mix(F0, albedo, metallic);
    }
    else if(res.y == 2.0)
    {
        roughness =0.1;
		roughnessD = roughness;
        metallic = 0.9;

		albedoKa = copper;
		albedoKd = copper;
		albedoKs = copper;
        
       	F0 = mix(F0, albedo, metallic);
    }
    else if (res.y == 4.0)
    {
        roughness =0.1;
		roughnessD = roughness;
        metallic = 0.9;
        albedoKa = gold;
		albedoKd = gold;
		albedoKs = gold;
        F0 = mix(F0, albedo, metallic);
    }
    else if (res.y == 5.0){
       	roughness =0.5;
		roughnessD = roughness;
        metallic = 0.5;
     
		albedoKa = plastic;
		albedoKd = plastic;
		albedoKs = plastic;
        F0 = mix(F0, albedo, metallic); 
    }

	//a few precomputed variable
vec3 N = normal(sampleP);
vec3 V = normalize(camPos - sampleP);

//F0 = mix(F0, albedoKs, metallic); 
vec3 Lo = vec3(0.0); //declare final return color here

vec3 L = normalize(light_pos - sampleP);
vec3 H = normalize(V+L);

 float attenDist = length(light_pos-sampleP);
 float attenuation = 1.0/(attenDist*attenDist);        
vec3 radiance = vec3(Ls.xyz) * attenuation;

//float D = DistributionGGX(N, H, roughnessD);
//float G = GeometrySmith(N, V, L, roughnessD);
//vec3 F = fresnelSchlick(max(dot(H,V), 0.0), F0);
//
float D = DistributionBeckmann(N, H, roughnessD);
float G = GeometryCookTorrance(N, V, L, H);
vec3 F = fresnelCookTorrance(max(dot(N,V), 0.0), F0);


vec3 amb  = vec3 (La)*albedoKa* 0.1f;
vec3 diff = vec3 (Ld)*albedoKd* max(dot(N, L), 0.0) ;
vec3 spec = vec3 (Ls)*albedoKs*(F* D* G)/(PI*dot(N, V)+0.5f);
//spec = mix(diff, spec, smoothstep(0.2f, 1.0f, dot(N,L)));


if(res.y == 3.0) 
 {
	if (switchTerm ==1) {
		spec = vec3 (Ls)*albedoKs*F;
		amb = vec3 (0.0f);
		diff = vec3 (0.0f);
		}
	if (switchTerm ==2) {
		spec = vec3 (Ls)*albedoKs*D;
		amb = vec3 (0.0f);
		diff = vec3 (0.0f);
		}
	if (switchTerm ==3) {
		spec = vec3 (Ls)*albedoKs*G;
		amb = vec3 (0.0f);
		diff = vec3 (0.0f);
		}
}

//vec3 DFGUpper = D * G * F;
//float DFGLower = PI * max(dot(N, V), 0.0) * max(dot(N,L), 0.0)+ 0.001; 
// vec3 kS = F;
// vec3 kD = vec3(1.0) - kS;
// kD *= 1.0 - roughness;
//this one use the equation from learn openGL
//vec3 spec = (kD * albedo/PI + DFGUpper/DFGLower)*radiance* max(dot(N,L), 0.0); 

float shadow = softshadow(sampleP, L, 0.01, 5.0, 32);

//vec3 amb  = vec3 (La)*albedoKa* 0.1f;
//vec3 diff = vec3 (Ld)*albedoKd* max(dot(N, L), 0.0) ;
//vec3 spec = vec3 (Ls)*albedoKs*(F* D* G)/(PI*dot(N, V)+1.0f);


//if (max(dot(N,L), 0.3f) == 0.3f) spec = vec3( 0.0f);



 Lo = amb + (diff + spec) * shadow; 

 //**phong specular debug purpose
//vec3 R = normalize( reflect(-L, N));
 //Lo = amb + (diff + pow(dot(R, V), 8)) * shadow; 
 return Lo;
}

//******
//Compute lighting on the raycast surface using Phong lighting model
vec4 lighting(vec3 pos, vec3 rayDir)
{
	const vec3 light = normalize(light_pos-pos); //light direction from surface
	vec3 n = normal(pos);

	//vec4 La = clear_color(n);
	float diff = 0.0f;
	diff = 0.1* max(0.0, dot(n, light));
	vec3 Spec = vec3(0.0f);

	if(scene == 3) {
	
	return vec4(renderPBR(pos), 1.0);
	
	}
	else diff = max(0.0, dot(n, light));
	//float diff = max(0.0, dot(n, light));

	return La*Ka + Ld*Kd*diff ;	
}

//Clear Color Pattern==================================================================================
float map (vec3 p){
vec3 q = fract (p-2.0f)* 2.0f -1.0f;

return length(q) - 0.2f;

//return 0.0f;
}

float march(vec3 o, vec3 r){

float t = 0.0f;
for (int i =0; i< 32; ++i){
vec3 p = o + r*t;
float d = map(p);
t += d*0.5;
}

return t;
//return 0.0f;
}

vec4 clearCol(vec3 o, vec3 r){


return 0.0f;
}



vec4 clear_color(vec3 rayDir)
{	
	vec2 resolution = vec2(1280.0f, 720.0f);
	vec2 uv = gl_FragCoord.xy/ resolution;
	uv = 2.0f * uv - 1.0f;
	uv.x *= resolution.x/resolution.y;
	vec3 r = normalize(vec3(uv, 1.0f));
	vec3 o = vec3( 0.0f, time, time);

	float tt = time*0.25f;
	r.xz *= mat2(cos(tt), -sin(tt), sin(tt), cos(tt));

	float dist = march(o, r);
	
	float fog = 1.0f/ (1.0f + dist*dist*0.1f);

	if(fog > 0.3f && fog < 0.5f) fog = abs(sin(time));
	
//	return vec4(vec3(fog), 1.0f);
	
	vec3 col =  0.4f*cos( time+2.0f*uv.yxy+vec3(6.0f,4.0f,2.0f));
	return vec4( col, 1.0f);
}

//===============================================================================

vec4 raytracedcolor(vec3 rayStart, vec3 rayStop)
{
	const int MaxSamples = 1000000; //max number of steps along ray

	vec3 rayDir = normalize(rayStop-rayStart);	//ray direction unit vector
	float travel = distance(rayStop, rayStart);	
	float stepSize = travel/MaxSamples;	//initial raymarch step size
	vec3 pos = rayStart;				//position along the ray
	vec3 step = rayDir*stepSize;		//displacement vector along ray
	
	for (int i=0; i < MaxSamples && travel > 0.0; ++i, pos += step, travel -= stepSize)
	{
		float dist = sceneSDF(pos); //How far are we from the shape we are raycasting?

		//Distance tells us how far we can safely step along ray without intersecting surface
		stepSize =  0.3*dist;
		step = rayDir*stepSize;
		
		//Check distance, and if we are close then perform lighting
		const float eps = 1e-5;
		if(dist <= eps)
		{
			return lighting(pos, rayDir);
		}	
	}
	//If the ray never intersects the scene then output clear color
	return clear_color(rayDir);
}


//main=============================================================

void main(void)
{   
	if(pass == 1)
	{
		fragcolor = vec4((vpos), 1.0); //write cube positions to texture
	}
	else if(pass == 2) 
	{
		if(mode == 0) // for debugging: show backface colors
		{
			fragcolor = texelFetch(backfaces_tex, ivec2(gl_FragCoord), 0);
			return;
		}
		else if(mode == 1) // for debugging: show frontface colors
		{
			fragcolor = vec4((vpos), 1.0);
			return;
		}
		else // raycast
		{
			vec3 rayStart = vpos.xyz;
			vec3 rayStop = texelFetch(backfaces_tex, ivec2(gl_FragCoord.xy), 0).xyz;
			fragcolor = raytracedcolor(rayStart, rayStop);
		}
	}
}


