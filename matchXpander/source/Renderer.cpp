// Renderer.cpp: implementation of the Renderer class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int Renderer::main_window;		///referenced by GLUI to get back to the main window

//GLUI windows
GLUI * Renderer::glui_controls;
float Renderer::aspectRatio;


//Live Variables Used by GLUI and other UI
float * Renderer::obj_pos;
float * Renderer::obj_rot;


Renderer::Renderer()
{
	obj_pos = new float[3];
	obj_pos[0] = 0.0;
	obj_pos[1] = 0.0;
	obj_pos[2] = 0.0;

	obj_rot = new float[16];
	obj_rot[0] = 1.0;	obj_rot[1] = 0.0;	obj_rot[2] = 0.0;	obj_rot[3] = 0.0;
	obj_rot[4] = 0.0;	obj_rot[5] = 1.0;	obj_rot[6] = 0.0;	obj_rot[7] = 0.0;
	obj_rot[8] = 0.0;	obj_rot[9] = 0.0;	obj_rot[10] = 1.0;	obj_rot[11] = 0.0;
	obj_rot[12] = 0.0;	obj_rot[13] = 0.0;	obj_rot[14] = 0.0;	obj_rot[15] = 1.0;

}

Renderer::~Renderer()
{

	delete(obj_pos);
	delete(obj_rot);
//	delete(motifVar_GLUIsaveName);
//	delete(motifVar_GLUIloadName);


}

void Renderer::init(void)
{

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	glDepthFunc(GL_LEQUAL);									// The Type Of Depth Testing To Do	
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_TEXTURE_2D);								// Enable T
	GLfloat light_position[] = { 10.0, 10.0, 10.0, 0.0 };
    GLfloat model_ambient[] = { 0.8, 0.8, 0.8, 1.0 };
	glClearDepth(1.0f);										// Depth Buffer Setup
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, model_ambient);

    glClearColor(0.0, 0.0, 0.0, 0.0);
	glPointSize(POINTSIZE);
	glLineWidth(LINEWIDTH);

}

void Renderer::idle( void )
{
  /* According to the GLUT specification, the current window is 
     undefined during an idle callback.  So we need to explicitly change
     it if necessary */
	if ( glutGetWindow() != main_window ){
		glutSetWindow(main_window);  
	}

	glutPostRedisplay();

}

void Renderer::beginRendering(void)
{
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(200, 250);

	main_window = glutCreateWindow("Geometric Hashing for Proteins by Brian Chen");

	glutIgnoreKeyRepeat(1);
	glutDisplayFunc(display);
		////mouse commands
	glutMouseFunc(processMouse);
		////keyboard commands
	glutKeyboardFunc(normalKeys);
		//creates and registers GLUI windows
	setUpGLUI();

	init();
	glutMainLoop();
}

void Renderer::render(int renderMode){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30, aspectRatio, 0.1f, CLIPPINGVOL);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();

	///Transform for Basic rotation/translation
	///translation
	glTranslatef(0.0, 0.0, ZDISPLACEMENT);			///to make viewing easier
	glTranslatef(obj_pos[0], obj_pos[1], obj_pos[2]);
	///rotation
	glMultMatrixf(obj_rot);

	glPushMatrix();

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
	int i = 0;

	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(10.0f, 0.0f, 0.0f);
	glEnd();
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 10.0f, 0.0f);
	glEnd();
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 10.0f);
	glEnd();

//     AA 2.596000 -0.298000 -2.048000 68 CA HIS 158 HREQ
//     AA 1.693000 -10.406000 -8.850000 44 CA ARG 162 RVNEY
//     AA -3.488000 -2.836000 4.963000 21 CA TRP 181 WR
//     AA -6.807000 -2.564000 -2.386000 44 CA GLU 184 ESQG
	double * *source = new double*[3];
	source[0] = new double[3];
	source[0][0] = 2.596;
	source[0][1] = -.298;
	source[0][2] = -2.048;
	source[1] = new double[3];
	source[1][0] = 1.693;
	source[1][1] = -10.406;
	source[1][2] = -8.85;
	source[2] = new double[3];
	source[2][0] = -3.488;
	source[2][1] = -2.836;
	source[2][2] = 4.963;

//    AA 18.869000 83.968000 59.384000 -1 CA ARG 174
//    AA 19.705000 92.291000 65.207000 -1 CA GLU 171
//    AA 14.565000 82.128000 52.486000 -1 CA ARG 178
//    AA 18.089000 79.303000 56.159000 -1 CA GLY 176
	double ** target = new double*[3];
	target[0] = new double[3];
	target[0][0] = 18.869;
	target[0][1] = 83.968;
	target[0][2] = 59.384;
	target[1] = new double[3];
	target[1][0] = 19.705;
	target[1][1] = 92.291;
	target[1][2] = 65.207;
	target[2] = new double[3];
	target[2][0] = 14.565;
	target[2][1] = 82.128;
	target[2][2] = 52.386;

	double * centroid = new double[3];
	centroid[0] = -(target[0][0] + target[1][0] + target[2][0])/3;
	centroid[1] = -(target[0][1] + target[1][1] + target[2][1])/3;
	centroid[2] = -(target[0][2] + target[1][2] + target[2][2])/3;
	
	glTranslatef(centroid[0], centroid[1], centroid[2]);
	glPushMatrix();

//	double * newAlignment = alignThreePoints(source, target, 3);
//	double * newAlignment = alignThreePoints(source, target);
	double * newAlignment = NULL;	///to compile 

	double * tempVect;

	if(newAlignment[15] == 1){
		glColor3f(1.0, 0.0, 0.0);
	}
	else{
		glColor3f(0.0, 1.0, 0.0);
	}

	for(i = 0; i<3; i++){
		glPushMatrix();
		if(newAlignment[15] == -1){
			tempVect = source[i];	
		}
		else{
			tempVect = transformVector3x4(newAlignment, source[i]);
		}

//		glTranslatef( source[i][0], source[i][1], source[i][2] );
		glTranslatef( tempVect[0], tempVect[1], tempVect[2] );

		glutSolidSphere(1, 8, 8);

		//testMatrix(newAlignment);

		glPopMatrix();
		if(newAlignment[15] != -1){
			delete[](tempVect);
		}
	}
	glColor3f(0.0, 0.0, 1.0);
	for(i = 0; i<3; i++){
		glPushMatrix();
		glTranslatef( target[i][0], target[i][1], target[i][2] );
		glutSolidSphere(1, 8, 8);
		glutWireSphere(EXTMATCHTHRESHOLD, 8, 8);

		glPopMatrix();
	}
	glPopMatrix();

	delete[](newAlignment);
	for(i = 0; i<3; i++){
		delete[](source[i]);
		delete[](target[i]);	
	}
	delete[](source);
	delete[](target);
	delete[](centroid);
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

	glPopMatrix();

	glutSwapBuffers();

	glFlush ();

}

void Renderer::display(void)
{
	render(GL_RENDER);

}

////////////////////////////////////////////////
//////Used for mouse selection
void Renderer::pickObjects(int x, int y)
{

} 


void Renderer::processHits(GLint hits, GLuint buffer[])
{

}

void Renderer::reshape(int w, int h)
{
	if (h==0) h=1;					//prevent divide by zero
	aspectRatio = ( (GLfloat) w/(GLfloat) h );
	
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);

	glutPostRedisplay();

}

void Renderer::normalKeys(unsigned char key, int x, int y)
{

}

void Renderer::processMouse(int button, int state, int x, int y)
{

}

void Renderer::controlCB(int control){

}

void Renderer::setUpGLUI(void)
{
	////////////////////////////////////////////////////////////////////////////////
	///Register Callbacks with GLUI/////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	//glutReshapeFunc(reshape);
	GLUI_Master.set_glutReshapeFunc( reshape );

	//glutIdleFunc(display);	
	GLUI_Master.set_glutIdleFunc( idle );

	////////////////////////////////////////
	//Create GLUI window
	glui_controls = GLUI_Master.create_glui( "Controls", 0, 807, 200);
	glui_controls->set_main_gfx_window( main_window );
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	///Viewing Controls//////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////

	//Translation Arrows
	GLUI_Translation *view_trans = 
		glui_controls->add_translation( "XY Control", GLUI_TRANSLATION_XY, obj_pos );
		view_trans->set_speed( .1 );
		view_trans->set_alignment( GLUI_ALIGN_RIGHT );

	glui_controls->add_column( false );

	//Zoom Arrows
	GLUI_Translation *trans_z = 
		glui_controls->add_translation( "Z Control", GLUI_TRANSLATION_Z, &obj_pos[2] );
		trans_z->set_speed( 5 );
		trans_z->set_alignment( GLUI_ALIGN_RIGHT );

	glui_controls->add_column( false );

	//Rotation ArcBall
	GLUI_Rotation *view_rot = 
		glui_controls->add_rotation( "Rotation", obj_rot );
		view_rot->set_alignment( GLUI_ALIGN_RIGHT );


}






