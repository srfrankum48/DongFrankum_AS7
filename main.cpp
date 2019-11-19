#include <stdio.h>
#include <stdlib.h>
#include <gl\glut.h>
#include <math.h>

#include "frame_buffer.h"
#include "primitives.h"
#include "color.h"
#include <vector>
#include "glprocs.h"

#include <iostream>
#include <fstream>

#define ON 1
#define OFF 0

using namespace std;


// Global variables
char *shaderFileRead(char *fn);
int window_width, window_height;    // Window dimensions
Scene* pDisplayScene;
Camera* pDisplayCamera;
GLuint vertex_shader, fragment_shader, p;


const int INITIAL_RES = 400;

FrameBuffer* fb;

class point
{
public:
	double x,y,z,w;

	point(){ x = 0; y = 0; z = 0; w = 1;}
	point(double xa, double ya, double za)
	{
		x = xa; y = ya; z = za; w = 1.0;
	}
	point(double xa, double ya, double za, double wa)
	{
		x = xa; y = ya; z = za; w = wa;
	}
};

Vertex::Vertex()
{
	x = y = z = 0;
	h = 1;
}

void Vertex::Normalize()
{
	x = x / h;
	y = y / h;
	z = z / h;
	h = 1;
}

Object::Object()
{
	// Load the identity for the initial modeling matrix
	ModelMatrix[0] = ModelMatrix[5] = ModelMatrix[10] = ModelMatrix[15] = 1;
	ModelMatrix[1] = ModelMatrix[2] = ModelMatrix[3] = ModelMatrix[4] =
		ModelMatrix[6] = ModelMatrix[7] = ModelMatrix[8] = ModelMatrix[9] =
		ModelMatrix[11] = ModelMatrix[12] = ModelMatrix[13] = ModelMatrix[14] = 0;
}

Object::~Object()
{
	delete[] pVertexList;
	delete[] pFaceList;
}

void Object::Load(char* file, float s, float rx, float ry, float rz,
	float tx, float ty, float tz)
{
	FILE* pObjectFile = fopen(file, "r");
	if (!pObjectFile)
		cout << "Failed to load " << file << "." << endl;
	else
		cout << "Successfully loaded " << file << "." << endl;

	char DataType;
	float a, b, c;

	// Scan the file and count the faces and vertices
	VertexCount = FaceCount = 0;
	while (!feof(pObjectFile))
	{
		fscanf(pObjectFile, "%c %f %f %f\n", &DataType, &a, &b, &c);
		if (DataType == 'v')
			VertexCount++;
		else if (DataType == 'f')
			FaceCount++;
	}
	pVertexList = new Vertex[VertexCount];
	pFaceList = new Face[FaceCount];

	fseek(pObjectFile, 0L, SEEK_SET);

	cout << "Number of vertices: " << VertexCount << endl;
	cout << "Number of faces: " << FaceCount << endl;

	// Load and create the faces and vertices
	int CurrentVertex = 0, CurrentFace = 0;
	while (!feof(pObjectFile))
	{
		fscanf(pObjectFile, "%c %f %f %f\n", &DataType, &a, &b, &c);
		if (DataType == 'v')
		{
			pVertexList[CurrentVertex].x = a;
			pVertexList[CurrentVertex].y = b;
			pVertexList[CurrentVertex].z = c;
			CurrentVertex++;
		}
		else if (DataType == 'f')
		{
			// Convert to a zero-based index for convenience
			pFaceList[CurrentFace].v1 = (int)a - 1;
			pFaceList[CurrentFace].v2 = (int)b - 1;
			pFaceList[CurrentFace].v3 = (int)c - 1;
			CurrentFace++;
		}
	}

	// Apply the initial transformations in order
	LocalScale(s);
	WorldRotate((float)(M_PI*rx / 180.0), (float)(M_PI*ry / 180.0), (float)(M_PI*rz / 180.0));
	WorldTranslate(tx, ty, tz);
}

void Object::Load(float Sx, float Sy, float Sz, float radius, float SAmbientR, float SAmbientG, float SAmbientB,
	float SDiffuseR, float SDiffuseG, float SDiffuseB, float SSpecularR, float SSpecularG, float SSpecularB,
	float Skambient, float Skdiffuse, float floatSkspecular, float SkspecularExp, float SIndexRefraction,
	float SkReflective, float SkRefractive)
{
	char *file = "Sphere.obj";
	FILE* pObjectFile = fopen(file, "r");
	if (!pObjectFile)
		cout << "Failed to load " << file << "." << endl;
	else
		cout << "Successfully loaded " << file << "." << endl;

	char DataType;
	float a, b, c;

	// Scan the file and count the faces and vertices
	VertexCount = FaceCount = 0;
	while (!feof(pObjectFile))
	{
		fscanf(pObjectFile, "%c %f %f %f\n", &DataType, &a, &b, &c);
		if (DataType == 'v')
			VertexCount++;
		else if (DataType == 'f')
			FaceCount++;
	}
	pVertexList = new Vertex[VertexCount];
	pFaceList = new Face[FaceCount];

	fseek(pObjectFile, 0L, SEEK_SET);

	cout << "Number of vertices: " << VertexCount << endl;
	cout << "Number of faces: " << FaceCount << endl;

	// Load and create the faces and vertices
	int CurrentVertex = 0, CurrentFace = 0;
	while (!feof(pObjectFile))
	{
		fscanf(pObjectFile, "%c %f %f %f\n", &DataType, &a, &b, &c);
		if (DataType == 'v')
		{
			pVertexList[CurrentVertex].x = a;
			pVertexList[CurrentVertex].y = b;
			pVertexList[CurrentVertex].z = c;
			CurrentVertex++;
		}
		else if (DataType == 'f')
		{
			// Convert to a zero-based index for convenience
			pFaceList[CurrentFace].v1 = (int)a - 1;
			pFaceList[CurrentFace].v2 = (int)b - 1;
			pFaceList[CurrentFace].v3 = (int)c - 1;
			CurrentFace++;
		}
	}
}

// Do world-based translation
void Object::WorldTranslate(float tx, float ty, float tz)
{
	ModelMatrix[12] += tx;
	ModelMatrix[13] += ty;
	ModelMatrix[14] += tz;
}

// Perform world-based rotations in x,y,z order (intended for one-at-a-time use)
void Object::WorldRotate(float rx, float ry, float rz)
{
	float temp[16];

	if (rx != 0)
	{
		float cosx = cos(rx), sinx = sin(rx);
		for (int i = 0; i < 16; i++)
			temp[i] = ModelMatrix[i];
		ModelMatrix[1] = temp[1] * cosx - temp[2] * sinx;
		ModelMatrix[2] = temp[2] * cosx + temp[1] * sinx;
		ModelMatrix[5] = temp[5] * cosx - temp[6] * sinx;
		ModelMatrix[6] = temp[6] * cosx + temp[5] * sinx;
		ModelMatrix[9] = temp[9] * cosx - temp[10] * sinx;
		ModelMatrix[10] = temp[10] * cosx + temp[9] * sinx;
		ModelMatrix[13] = temp[13] * cosx - temp[14] * sinx;
		ModelMatrix[14] = temp[14] * cosx + temp[13] * sinx;
	}

	if (ry != 0)
	{
		float cosy = cos(ry), siny = sin(ry);
		for (int i = 0; i < 16; i++)
			temp[i] = ModelMatrix[i];
		ModelMatrix[0] = temp[0] * cosy + temp[2] * siny;
		ModelMatrix[2] = temp[2] * cosy - temp[0] * siny;
		ModelMatrix[4] = temp[4] * cosy + temp[6] * siny;
		ModelMatrix[6] = temp[6] * cosy - temp[4] * siny;
		ModelMatrix[8] = temp[8] * cosy + temp[10] * siny;
		ModelMatrix[10] = temp[10] * cosy - temp[8] * siny;
		ModelMatrix[12] = temp[12] * cosy + temp[14] * siny;
		ModelMatrix[14] = temp[14] * cosy - temp[12] * siny;
	}

	if (rz != 0)
	{
		float cosz = cos(rz), sinz = sin(rz);
		for (int i = 0; i < 16; i++)
			temp[i] = ModelMatrix[i];
		ModelMatrix[0] = temp[0] * cosz - temp[1] * sinz;
		ModelMatrix[1] = temp[1] * cosz + temp[0] * sinz;
		ModelMatrix[4] = temp[4] * cosz - temp[5] * sinz;
		ModelMatrix[5] = temp[5] * cosz + temp[4] * sinz;
		ModelMatrix[8] = temp[8] * cosz - temp[9] * sinz;
		ModelMatrix[9] = temp[9] * cosz + temp[8] * sinz;
		ModelMatrix[12] = temp[12] * cosz - temp[13] * sinz;
		ModelMatrix[13] = temp[13] * cosz + temp[12] * sinz;
	}
}

// Perform locally-based rotations in x,y,z order (intended for one-at-a-time use)
void Object::LocalRotate(float rx, float ry, float rz)
{
	float temp[16];

	if (rx != 0)
	{
		float cosx = cos(rx), sinx = sin(rx);
		for (int i = 0; i < 16; i++)
			temp[i] = ModelMatrix[i];
		ModelMatrix[4] = temp[4] * cosx + temp[8] * sinx;
		ModelMatrix[5] = temp[5] * cosx + temp[9] * sinx;
		ModelMatrix[6] = temp[6] * cosx + temp[10] * sinx;
		ModelMatrix[7] = temp[7] * cosx + temp[11] * sinx;
		ModelMatrix[8] = temp[8] * cosx - temp[4] * sinx;
		ModelMatrix[9] = temp[9] * cosx - temp[5] * sinx;
		ModelMatrix[10] = temp[10] * cosx - temp[6] * sinx;
		ModelMatrix[11] = temp[11] * cosx - temp[7] * sinx;
	}

	if (ry != 0)
	{
		float cosy = cos(ry), siny = sin(ry);
		for (int i = 0; i < 16; i++)
			temp[i] = ModelMatrix[i];
		ModelMatrix[0] = temp[0] * cosy - temp[8] * siny;
		ModelMatrix[1] = temp[1] * cosy - temp[9] * siny;
		ModelMatrix[2] = temp[2] * cosy - temp[10] * siny;
		ModelMatrix[3] = temp[3] * cosy - temp[11] * siny;
		ModelMatrix[8] = temp[8] * cosy + temp[0] * siny;
		ModelMatrix[9] = temp[9] * cosy + temp[1] * siny;
		ModelMatrix[10] = temp[10] * cosy + temp[2] * siny;
		ModelMatrix[11] = temp[11] * cosy + temp[3] * siny;
	}

	if (rz != 0)
	{
		float cosz = cos(rz), sinz = sin(rz);
		for (int i = 0; i < 16; i++)
			temp[i] = ModelMatrix[i];
		ModelMatrix[0] = temp[0] * cosz + temp[4] * sinz;
		ModelMatrix[1] = temp[1] * cosz + temp[5] * sinz;
		ModelMatrix[2] = temp[2] * cosz + temp[6] * sinz;
		ModelMatrix[3] = temp[3] * cosz + temp[7] * sinz;
		ModelMatrix[4] = temp[4] * cosz - temp[0] * sinz;
		ModelMatrix[5] = temp[5] * cosz - temp[1] * sinz;
		ModelMatrix[6] = temp[6] * cosz - temp[2] * sinz;
		ModelMatrix[7] = temp[7] * cosz - temp[3] * sinz;
	}
}

// Do locally-based uniform scaling
void Object::LocalScale(float s)
{
	for (int i = 0; i <= 11; i++)
		ModelMatrix[i] = s * ModelMatrix[i];
}

void Scene::Load(char* file)
{
	FILE* pSceneFile = fopen(file, "r");
	if (!pSceneFile)
		cout << "Failed to load " << file << "." << endl;
	else
		cout << "Successfully loaded " << file << "." << endl;

	char MeshFile[255];
	float NumLights, NumSpheres, NumTriangleMeshes;
	float LightType, Lx, Ly, Lz, Lr, Lg, Lb;
	float Sx, Sy, Sz, radius, SAmbientR, SAmbientG, SAmbientB,
		SDiffuseR, SDiffuseG, SDiffuseB, SSpecularR, SSpecularG, SSpecularB,
		Skambient, Skdiffuse, Skspecular, SkspecularExp, SIndexRefraction,
		SkReflective, SkRefractive;
	float Scaling, RotationX, RotationY, RotationZ,
		TranslationX, TranslationY, TranslationZ;

	// Step through the file and count the objects
	ObjectCount = 0;
	while (!feof(pSceneFile))
	{
		char letter;
		fscanf(pSceneFile, "%c", &letter);
		if (letter == 'L') {
			fscanf(pSceneFile, "%f %f %f %f %f %f %f\n", &LightType, &Lx, &Ly, &Lz,
				&Lr, &Lg, &Lb);
		}
		else if (letter == 'S') {
			fscanf(pSceneFile, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
				&Sx, &Sy, &Sz, &radius, &SAmbientR, &SAmbientG, &SAmbientB, &SDiffuseR, &SDiffuseG, &SDiffuseB,
				&SSpecularR, &SSpecularG, &SSpecularB, &Skambient, &Skdiffuse, &Skspecular,
				&SkspecularExp, &SIndexRefraction, &SkReflective, &SkRefractive);
			ObjectCount++;
		}
		else if (letter == 'M') {
			fscanf(pSceneFile, "%s %f %f %f %f %f %f %f\n", MeshFile, &Scaling,
				&RotationX, &RotationY, &RotationZ, &TranslationX, &TranslationY, &TranslationZ);
			ObjectCount++;
		}
		else {
			fscanf(pSceneFile, "%f %f %f\n", &NumLights, &NumSpheres, &NumTriangleMeshes);
		}
	}
	pObjectList = new Object[ObjectCount];

	fseek(pSceneFile, 0L, SEEK_SET);

	// Step through the file and create/load the objects
	for (int i = 0; i < ObjectCount; i++)
	{
		char letter= ' ';
		fscanf(pSceneFile, "%c", &letter);
		if (letter == 'S') {
			fscanf(pSceneFile, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
				&Sx, &Sy, &Sz, &radius, &SAmbientR, &SAmbientG, &SAmbientB, &SDiffuseR, &SDiffuseG, &SDiffuseB,
				&SSpecularR, &SSpecularG, &SSpecularB, &Skambient, &Skdiffuse, &Skspecular,
				&SkspecularExp, &SIndexRefraction, &SkReflective, &SkRefractive);
			pObjectList[i].Load(Sx, Sy, Sz, radius, SAmbientR, SAmbientG, SAmbientB, SDiffuseR, SDiffuseG, SDiffuseB,
				SSpecularR, SSpecularG, SSpecularB, Skambient, Skdiffuse, Skspecular,
				SkspecularExp, SIndexRefraction, SkReflective, SkRefractive);
		}
		else if (letter == 'M') {
			fscanf(pSceneFile, "%s %f %f %f %f %f %f %f\n", MeshFile, &Scaling,
				&RotationX, &RotationY, &RotationZ, &TranslationX, &TranslationY, &TranslationZ);
			pObjectList[i].Load(MeshFile, Scaling, RotationX, RotationY, RotationZ,
				TranslationX, TranslationY, TranslationZ);
		}
	}

	cout << "Number of Objects Loaded: " << ObjectCount << endl;
}


typedef struct _faceStruct {
  int v1,v2,v3;
  int n1,n2,n3;
} faceStruct;

int verts, faces, norms;    // Number of vertices, faces and normals in the system
point *vertList, *normList; // Vertex and Normal Lists
faceStruct *faceList;	    // Face List

// The mesh reader itself
// It can read *very* simple obj files
void meshReader (char *filename,int sign)
{
  float x,y,z,len;
  int i;
  char letter;
  point v1,v2,crossP;
  int ix,iy,iz;
  int *normCount;
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) { 
    printf("Cannot open %s\n!", filename);
    exit(0);
  }

  // Count the number of vertices and faces
  while(!feof(fp))
    {
      fscanf(fp,"%c %f %f %f\n",&letter,&x,&y,&z);
      if (letter == 'v')
	{
	  verts++;
	}
      else
	{
	  faces++;
	}
    }

  fclose(fp);

  printf("verts : %d\n", verts);
  printf("faces : %d\n", faces);

  // Dynamic allocation of vertex and face lists
  faceList = (faceStruct *)malloc(sizeof(faceStruct)*faces);
  vertList = (point *)malloc(sizeof(point)*verts);
  normList = (point *)malloc(sizeof(point)*verts);

  fp = fopen(filename, "r");

  // Read the veritces
  for(i = 0;i < verts;i++)
    {
      fscanf(fp,"%c %f %f %f\n",&letter,&x,&y,&z);
      vertList[i].x = x;
      vertList[i].y = y;
      vertList[i].z = z;
    }

  // Read the faces
  for(i = 0;i < faces;i++)
    {
      fscanf(fp,"%c %d %d %d\n",&letter,&ix,&iy,&iz);
      faceList[i].v1 = ix - 1;
      faceList[i].v2 = iy - 1;
      faceList[i].v3 = iz - 1;
    }
  fclose(fp);


  // The part below calculates the normals of each vertex
  normCount = (int *)malloc(sizeof(int)*verts);
  for (i = 0;i < verts;i++)
    {
      normList[i].x = normList[i].y = normList[i].z = 0.0;
      normCount[i] = 0;
    }

  for(i = 0;i < faces;i++)
    {
      v1.x = vertList[faceList[i].v2].x - vertList[faceList[i].v1].x;
      v1.y = vertList[faceList[i].v2].y - vertList[faceList[i].v1].y;
      v1.z = vertList[faceList[i].v2].z - vertList[faceList[i].v1].z;
      v2.x = vertList[faceList[i].v3].x - vertList[faceList[i].v2].x;
      v2.y = vertList[faceList[i].v3].y - vertList[faceList[i].v2].y;
      v2.z = vertList[faceList[i].v3].z - vertList[faceList[i].v2].z;

      crossP.x = v1.y*v2.z - v1.z*v2.y;
      crossP.y = v1.z*v2.x - v1.x*v2.z;
      crossP.z = v1.x*v2.y - v1.y*v2.x;

      len = sqrt(crossP.x*crossP.x + crossP.y*crossP.y + crossP.z*crossP.z);

      crossP.x = -crossP.x/len;
      crossP.y = -crossP.y/len;
      crossP.z = -crossP.z/len;

      normList[faceList[i].v1].x = normList[faceList[i].v1].x + crossP.x;
      normList[faceList[i].v1].y = normList[faceList[i].v1].y + crossP.y;
      normList[faceList[i].v1].z = normList[faceList[i].v1].z + crossP.z;
      normList[faceList[i].v2].x = normList[faceList[i].v2].x + crossP.x;
      normList[faceList[i].v2].y = normList[faceList[i].v2].y + crossP.y;
      normList[faceList[i].v2].z = normList[faceList[i].v2].z + crossP.z;
      normList[faceList[i].v3].x = normList[faceList[i].v3].x + crossP.x;
      normList[faceList[i].v3].y = normList[faceList[i].v3].y + crossP.y;
      normList[faceList[i].v3].z = normList[faceList[i].v3].z + crossP.z;
      normCount[faceList[i].v1]++;
      normCount[faceList[i].v2]++;
      normCount[faceList[i].v3]++;
    }
  for (i = 0;i < verts;i++)
    {
      normList[i].x = (float)sign*normList[i].x / (float)normCount[i];
      normList[i].y = (float)sign*normList[i].y / (float)normCount[i];
      normList[i].z = (float)sign*normList[i].z / (float)normCount[i];
    }

}


void setParameters(GLuint program) {
	//The parameters in quotes are the names of the corresponding variables
	ambient_loc = glGetUniformLocationARB(program, "AmbientContribution");
	glUniform3fvARB(ambient_loc, 1, ambient_cont);

	diffuse_loc = glGetUniformLocationARB(program, "DiffuseContribution");
	glUniform3fvARB(diffuse_loc, 1, diffuse_cont);

	specular_loc = glGetUniformLocationARB(program, "SpecularContribution");
	glUniform3fvARB(specular_loc, 1, specular_cont);

	exponent_loc = glGetUniformLocationARB(program, "exponent");
	glUniform1fARB(exponent_loc, exponent);

	//Access attributes in vertex shader
	tangent_loc = glGetAttribLocationARB(program, "tang");
	glVertexAttrib1fARB(tangent_loc, tangent);
}


void drawRect(double x, double y, double w, double h)
{
	glVertex2f(x,y);
	glVertex2f(x+w,y);
	glVertex2f(x+w,y+h);
	glVertex2f(x, y+h);
}

// Transform a point with an arbitrary matrix
Vertex Transform(float* matrix, Vertex& point)
{
	Vertex temp;
	temp.x = matrix[0] * point.x + matrix[4] * point.y + matrix[8] * point.z + matrix[12] * point.h;
	temp.y = matrix[1] * point.x + matrix[5] * point.y + matrix[9] * point.z + matrix[13] * point.h;
	temp.z = matrix[2] * point.x + matrix[6] * point.y + matrix[10] * point.z + matrix[14] * point.h;
	temp.h = matrix[3] * point.x + matrix[7] * point.y + matrix[11] * point.z + matrix[15] * point.h;
	return temp;
}

void setShaders(char v[], char f[]) {

	char *vs = NULL, *fs = NULL;

	//create the empty shader objects and get their handles
	vertex_shader = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
	fragment_shader = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);


	//read the shader files and store the strings in corresponding char. arrays.
	vs = shaderFileRead(v);
	fs = shaderFileRead(f);

	const char * vv = vs;
	const char * ff = fs;

	//set the shader's source code by using the strings read from the shader files.
	glShaderSourceARB(vertex_shader, 1, &vv, NULL);
	glShaderSourceARB(fragment_shader, 1, &ff, NULL);

	free(vs); free(fs);

	//Compile the shader objects
	glCompileShaderARB(vertex_shader);
	glCompileShaderARB(fragment_shader);


	//create an empty program object to attach the shader objects
	p = glCreateProgramObjectARB();

	//attach the shader objects to the program object
	glAttachObjectARB(p, vertex_shader);
	glAttachObjectARB(p, fragment_shader);

	/*
	**************
	Programming Tip:
	***************
	Delete the attached shader objects once they are attached.
	They will be flagged for removal and will be freed when they are no more used.
	*/
	glDeleteObjectARB(vertex_shader);
	glDeleteObjectARB(fragment_shader);

	//Link the created program.
	/*
	**************
	Programming Tip:
	***************
	You can trace the status of link operation by calling
	"glGetObjectParameterARB(p,GL_OBJECT_LINK_STATUS_ARB)"
	*/
	glLinkProgramARB(p);

	//Start to use the program object, which is the part of the current rendering state
	glUseProgramObjectARB(p);

}

// The display function. It is called whenever the window needs
// redrawing (ie: overlapping window moves, resize, maximize)
// You should redraw your polygons here
void	display(void)
{
    // Clear the background
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	double w = 10/double(fb->GetWidth());
	double h = 10/double(fb->GetHeight());

	Color cl;
	glColor3f(0,0,1);

	glBegin(GL_QUADS);

	printf("width %d, height %d\n", fb->GetWidth(), fb->GetHeight());

	for(int y = 0; y < fb->GetHeight(); y++)
	{
		for(int x = 0; x < fb->GetHeight(); x++)
		{
			cl = fb->buffer[x][y].color;
			glColor3f(cl.r, cl.g, cl.b);

			drawRect(w*x, h*y, w, h);
		}
	}

	Vertex* input;
	Vertex	temp, temp1, temp2, temp3;
	for (int i = 0; i < pDisplayScene->ObjectCount; i++)
	{
		for (int j = 0; j < pDisplayScene->pObjectList[i].FaceCount; j++)
		{
			input = new Vertex[3];
			input[0] = pDisplayScene->pObjectList[i].pVertexList[pDisplayScene->pObjectList[i].pFaceList[j].v1];
			input[1] = pDisplayScene->pObjectList[i].pVertexList[pDisplayScene->pObjectList[i].pFaceList[j].v2];
			input[2] = pDisplayScene->pObjectList[i].pVertexList[pDisplayScene->pObjectList[i].pFaceList[j].v3];


			for (int k = 0; k < 3; k++) {
				temp = Transform(pDisplayScene->pObjectList[i].ModelMatrix, input[k]);
				temp2 = Transform(pDisplayCamera->ViewingMatrix, temp);
				input[k] = Transform(pDisplayCamera->ProjectionMatrix, temp2);
			}

			glBegin(GL_POLYGON);
			for (int k = 0; k < 3; k++)
				glVertex2f(input[k].x / input[k].h, input[k].y / input[k].h);
			glEnd();

			delete[] input;
			input = NULL;
		}
	}

	glEnd();
    glutSwapBuffers();
}


// This function is called whenever the window is resized. 
// Parameters are the new dimentions of the window
void	resize(int x,int y)
{
    glViewport(0,0,x,y);
    window_width = x;
    window_height = y;
    
    printf("Resized to %d %d\n",x,y);
}


// This function is called whenever the mouse is pressed or released
// button is a number 0 to 2 designating the button
// state is 1 for release 0 for press event
// x and y are the location of the mouse (in window-relative coordinates)
void	mouseButton(int button,int state,int x,int y)
{
   ;
}


//This function is called whenever the mouse is moved with a mouse button held down.
// x and y are the location of the mouse (in window-relative coordinates)
void	mouseMotion(int x, int y)
{
	;
}

int new_x, new_y;

// This function is called whenever there is a keyboard input
// key is the ASCII value of the key pressed
// x and y are the location of the mouse
void	keyboard(unsigned char key, int x, int y)
{
    switch(key) {
    case 'q':                           /* Quit */
		exit(1);
		break;
	case '-':
		new_x = fb->GetWidth() / 2;
		new_y = fb->GetHeight() / 2;
		fb->Resize(new_y, new_x);
		resize(new_y, new_x);
		BresenhamLine(fb, fb->GetWidth()*0.1, fb->GetHeight()*0.1, fb->GetWidth()*0.9, fb->GetHeight()*0.9, Color(1,0,0));
		break;
	case '=':
		new_x = fb->GetWidth() * 2;
		new_y = fb->GetHeight() * 2;
		fb->Resize(new_y, new_x);
		resize(new_y, new_x);
		BresenhamLine(fb, fb->GetWidth()*0.1, fb->GetHeight()*0.1, fb->GetWidth()*0.9, fb->GetHeight()*0.9, Color(1,0,0));
		break;
	case ']': // move image plane farther to origin (zaxis)
		for (int i = 0; i < fb->x_res; i++) {
			for (int j = 0; j < fb->y_res; j++) {
				fb->SetPixel(i, j, fb->GetPixel(i, j).color, fb->GetPixel(i, j).z_value + 1);
			}
		}
		display();
		break;
	case '[': // move image plane closer to origin (zaxis)
		for (int i = 0; i < fb->GetWidth(); i++) {
			for (int j = 0; j < fb->GetHeight(); j++) {
				fb->SetPixel(i, j, fb->GetPixel(i, j).color, fb->GetPixel(i, j).z_value - 1);
			}
		}
		display();
		break;
	case '.':
		break;
	case ',':
		break;
	case 'r':
		setParams();
		setShaders("PhongShader.frag", "PhongShader.vert");
		glutPostRedisplay();
		break;
    default:
		break;
    }

    // Schedule a new display event
    glutPostRedisplay();
}


int main(int argc, char* argv[])
{    

	fb = new FrameBuffer(INITIAL_RES, INITIAL_RES);

	BresenhamLine(fb, fb->GetWidth()*0.1, fb->GetHeight()*0.1, fb->GetWidth()*0.9, fb->GetHeight()*0.9, Color(1,0,0));

    // Initialize GLUT
	pDisplayScene = new Scene;
	pDisplayScene->Load("red_sphere_and_teapot.rtl");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Raytracer");
    glutDisplayFunc(display);
    glutReshapeFunc(resize);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);
    glutKeyboardFunc(keyboard);

	setShaders("PhongShader.frag", "PhongShader.vert");


	glEnable(GL_LIGHTING);

	// ambient light
	glEnable(GL_LIGHT0);

	GLfloat light0_pos[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat light0_a[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	GLfloat light0_d[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	GLfloat light0_s[] = { 0.2f, 0.2f, 0.2f, 1.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_a);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_d);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_s);

	glEnable(GL_NORMALIZE);

    // Initialize GL
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,10,0,10,-10000,10000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glEnable(GL_DEPTH_TEST);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

    // Switch to main loop
    glutMainLoop();
    return 0;        
}

//Read the shader files, given as parameter.
char *shaderFileRead(char *fn) {


	FILE *fp = fopen(fn, "r");
	if (!fp)
	{
		cout << "Failed to load " << fn << endl;
		return " ";
	}
	else
	{
		cout << "Successfully loaded " << fn << endl;
	}

	char *content = NULL;

	int count = 0;

	if (fp != NULL)
	{
		fseek(fp, 0, SEEK_END);
		count = ftell(fp);
		rewind(fp);

		if (count > 0)
		{
			content = (char *)malloc(sizeof(char) * (count + 1));
			count = fread(content, sizeof(char), count, fp);
			content[count] = '\0';
		}
		fclose(fp);
	}
	return content;
}
