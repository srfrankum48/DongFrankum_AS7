uniform vec4 lightDir;

varying vec3  vNormal;
varying vec3  vLightVec;
varying vec3  vViewVec;

void main(void)
{
   // Output transformed vertex position:
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;

   // Compute the light vector (view space):
   vLightVec   = -lightDir.xyz;
   vLightVec.z = -vLightVec.z;    // notice that we need to flip the z here to let OGL match D3D

   // Transform vertex position into view space:
   vec3 Pview =  vec3(gl_ModelViewMatrix * gl_Vertex);

   // Transform normal into view space:        
   vNormal = normalize( gl_NormalMatrix * gl_Normal);

   // Compute view vector (view space):
   vViewVec = -normalize(Pview);  
   
   vTexCoord = vec2(gl_MultiTexCoord0);
}