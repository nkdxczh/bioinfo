// Renderer.h: interface for the Renderer class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_GEOHASH_RENDERER_)
#define _GEOHASH_RENDERER_

class Renderer  
{
public:
	Renderer();
	virtual ~Renderer();
	static int main_window;		///Return pathway for GLUI to OpenGL Rendering Window

	////Glut Callbacks
	static void init(void);
	static void idle(void);
	static void display(void);
	static void reshape(int w, int h);

	//Keyboard input callbacks
	static void normalKeys(unsigned char key, int x, int y);
	static void processMouse(int button, int state, int x, int y);

	//Picking callbacks and functions
	static void pickObjects(int x, int y);
	static void processHits(GLint hits, GLuint buffer[]);

	//Rendering Functions and callback
	static void render(int mode);
	static void beginRendering(void);

	//GLUI activators and callback
	static void setUpGLUI(void);
	static void controlCB(int control);

	//GLUI windows
	static GLUI * glui_controls;

	//Live variables for GLUI and other User Interface
	static float * obj_pos;
	static float * obj_rot;
	static float aspectRatio;
};

#endif


