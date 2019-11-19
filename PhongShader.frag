uniform float n_specular;
uniform float Ks;
uniform float Ka;
uniform float Kd;
uniform vec4 diffuse;
uniform vec4 specular;
uniform vec4 ambient;

varying vec3  vNormal;
varying vec3  vLightVec;
varying vec3  vViewVec;

void solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
   float discr = b * b - 4 * a * c;
   
   if (discr < 0) 
      return 0;
   else if (discr == 0) 
      x0 = x1 = -0.5 * b / a;
   else {
      float q = (b > 0) ?
         -0.5 * (b + sqrt(discr)) :
         -0.5 * (b - sqrt(discr));
      x0 = q / a;
      x1 = c / q;
   }
   if (x0 > x1) std::swap(x0, x1);

   return 1;
}

void sphereIntersect(const vec3 ray) 
{
   float t0, t1; // solutions for t if the ray intersects 
#if 0 
   // geometric solution
   vec3 L = center - orig;
   float tca = L.dotProduct(dir);
   // if (tca < 0) return false;
   float d2 = L.dotProduct(L) - tca * tca;
   if (d2 > radius2) 
      return false;
   float thc = sqrt(radius2 - d2);
   t0 = tca - thc;
   t1 = tca + thc;
#else 
      // analytic solution
   vec3 L = orig - center;
   float a = dir.dotProduct(dir);
   float b = 2 * dir.dotProduct(L);
   float c = L.dotProduct(L) - radius2;
   if (!solveQuadratic(a, b, c, t0, t1)) 
      return false;
#endif 
   if (t0 > t1) {
      float temp = t0;
      t0 = t1;
      t1 = temp;
   }
   if (t0 < 0) {
      t0 = t1; // if t0 is negative, let's use t1 instead 
      if (t0 < 0) return false; // both t0 and t1 are negative 
   }

   t = t0;

   return true;
}

void main(void)
{
  // Compute the reflection vector:
   vec3 vReflect = normalize( 2.0 * dot( vNormal, vLightVec) * vNormal - vLightVec );       

   // Compute ambient term:
   vec4 AmbientColor = ambient * Ka;

   // Compute diffuse term:
   vec4 DiffuseColor = diffuse * Kd * max( 0.0, dot( vNormal, vLightVec ));

   // Compute specular term:
   vec4 SpecularColor = specular * Ks * pow( max( 0.0, dot(vReflect, vViewVec)), n_specular );
   
   vec4 FinalColor = (AmbientColor + DiffuseColor) + SpecularColor;
   
   gl_FragColor = FinalColor;
}