#include "glcommon.h"
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

GLenum printGLError_helper(const char * file,int line){
  static GLenum lastErr=0;
  static int lastLine=0;
  static std::string lastFile;
  static int times=0;
  GLenum err=glGetError();
  if(err){
    if(line==lastLine && lastErr==err && lastFile==file)
      times++;
    else times=0;
#define ERROR_CASE(X) case X:cerr << #X;break;
    
    if(times==0||times%100==0){
      cerr << "glError==" << err << "(at " << 
	file << ":" << line << ")";
      switch(err){
        ERROR_CASE(GL_INVALID_ENUM);
        ERROR_CASE(GL_TABLE_TOO_LARGE);
        ERROR_CASE(GL_OUT_OF_MEMORY);
        ERROR_CASE(GL_STACK_UNDERFLOW);
        ERROR_CASE(GL_STACK_OVERFLOW);
        ERROR_CASE(GL_INVALID_OPERATION);
        ERROR_CASE(GL_INVALID_VALUE);
      default:
	cerr << "unknown";
        break;
      }
      if(times)cerr << " (" << times << " times)" << endl;
      else cerr << "\n";
    }
#undef ERROR_CASE
    if(times==0){
      times=0;
      lastErr=err;
      lastFile=file;
      lastLine=line;
    }
  }
  return err;
}

void cameraPerspectiveMatrix(double fx,double fy,
			     double px,double py,
			     int imw,int imh,
			     double minz,double maxz,
			     int origintl)
{
  double mat[16];
  memset(mat,0,sizeof(mat));

  mat[0]=fx*2/imw;
  
  //Aspect
  //mat[4]=intrinsics[4]*2/imw;
  mat[5]=fy*2/imh;

  mat[8]=px*2/imw-1;
  mat[9]=py*2/imh-1;
  mat[10]=(origintl?-1: 1)*(minz+maxz)/(minz-maxz);
  mat[11]=(origintl? 1:-1);
  
  mat[14]=2*maxz*minz/(minz-maxz);

  glMultMatrixd(mat);
}



time_t getModifiedTime(const std::string & fname){
  struct stat sbuf;
  if(!stat(fname.c_str(),&sbuf)){
    return  sbuf.st_mtime;
  }
  return 0;
}


string resolvePath(const std::string & fname,
		   const std::vector<std::string> & paths){
  struct stat sbuf;
  if(!stat(fname.c_str(),&sbuf)){
    return fname;
  }
  //paths assumed to have no "/"
  for(int i=0; i<(int)paths.size(); i++){
    string name=paths[i]+"/"+fname;
    if(!stat(name.c_str(),&sbuf)){
       return name;
    }
  }
  return "";
}

void printShaderInfoLog(GLuint shader,const char * name){
  GLint logLength=0;
  glGetShaderiv(shader,GL_INFO_LOG_LENGTH,&logLength);
  
  if(logLength){
    GLint sz=0;
    GLchar * info=new GLchar[logLength+1];
    glGetShaderInfoLog(shader,logLength,&sz,info);
    printf("Shader %s Info: %s\n",name,info);
    delete [] info;
  }
}

void printProgramInfoLog(GLuint handle,const char * name){
  GLint logLength=0;
  glGetProgramiv(handle,GL_INFO_LOG_LENGTH,&logLength);

  if(logLength){
    GLint sz=0;
    GLchar * info=new GLchar[logLength+1];
    glGetProgramInfoLog(handle,logLength,&sz,info);
    printf("Program %s Info: %s\n",name,info);
    delete [] info;
  }
}

time_t load_shader(const std::string & name,GLuint shader){
  FILE * file=fopen(name.c_str(),"r");
  if(!file)return 0;
  printf("opened shader file %s\n",name.c_str());

  struct stat sbuf;
  if(!fstat(fileno(file),&sbuf)){
    char * data=new char[sbuf.st_size];
    GLint len=fread(data,1,sbuf.st_size,file);
    if(len!=sbuf.st_size){
      printf("error loading entire file:%s\n",name.c_str());
    }
    printGLError();
    const GLchar * lines[1]={data};

    glShaderSource(shader,1,lines,&len);
    glCompileShader(shader);

    GLint compileStatus;
    glGetShaderiv(shader,GL_COMPILE_STATUS,&compileStatus);
    if(!compileStatus){
      printf("compilation of shader %s was unsuccessful\n",name.c_str());
    }
    printShaderInfoLog(shader,name.c_str());

    delete [] data;
  }
  else printf("error statting shader %s\n",name.c_str());
  fclose(file);
  return time(0);
}


//////////////////////////////////////////////////////////////////////
// Shader program resource
//////////////////////////////////////////////////////////////////////

void checkLinkStatus(GLuint prog,const char * progname,
		     const char * filename,const char * function,int lineno){
  GLint linkStatus;
  glGetProgramiv(prog,GL_LINK_STATUS,&linkStatus);
  if(!linkStatus){
    printf("Program %s failed to link (%s:%s:%d)\n",progname,filename,function,lineno);
  }
}

//////////////////////////////////////////////////////////////////////
// Arb program resources
/////////////////////////////////////////////////////////////////////

time_t load_arb_program(const std::string & name,GLenum target){
  FILE * file=fopen(name.c_str(),"r");
  if(!file)return 0;

  struct stat sbuf;
  if(!fstat(fileno(file),&sbuf)){
    char * data=new char[sbuf.st_size];
    int len=fread(data,1,sbuf.st_size,file);
    printGLError();
    
    glProgramStringARB(target,GL_PROGRAM_FORMAT_ASCII_ARB,len,data);

    if(!printGLError())
      printf("load successful %d bytes\n",len);
    else{
      GLint errPos;
      glGetIntegerv(GL_PROGRAM_ERROR_POSITION_ARB,&errPos);
      printf("error in program at %d %s\n",errPos,glGetString(GL_PROGRAM_ERROR_STRING_ARB));
    }
    delete [] data;
  }
  else printf("error statting fragment/vertex program file\n");

  fclose(file);
  return time(0);
}

/////////////////////////////////////////////////////////////////
// Quick and dirty loading/linking of programs.
/////////////////////////////////////////////////////////////////

void setProgramDefaults(GLuint prog){
  GLint current;

  glGetIntegerv(GL_CURRENT_PROGRAM, &current);

  GLint loc;
  glUseProgram(prog);
  for(int i=0;i<16;i++){
    stringstream ss;
    ss<<"tex"<<i;
    loc=glGetUniformLocation(prog,ss.str().c_str());
    if(loc>=0)glUniform1i(loc,i);
  }
  loc=glGetUniformLocation(prog,"textures");
  GLint texs[16]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  if(loc>=0)glUniform1iv(loc,16,texs);

  glUseProgram(current);
}


bool link_program(GLuint prog, GLuint s1, GLuint s2, GLuint s3, GLuint s4){
  if(s1)glAttachShader(prog, s1);
  if(s2)glAttachShader(prog, s2);
  if(s3)glAttachShader(prog, s3);
  if(s4)glAttachShader(prog, s4);

  glLinkProgram(prog);

  checkLinkStatus(prog,"vertName",__FILE__,__FUNCTION__,__LINE__);
  printGLError();

  setProgramDefaults(prog);

  GLint linkStatus;
  glGetProgramiv(prog,GL_LINK_STATUS,&linkStatus);
  return linkStatus != 0;
}

GLuint create_linked_program(GLuint s1, GLuint s2, GLuint s3, GLuint s4){
  GLuint prog = glCreateProgram();
  link_program(prog, s1, s2, s3, s4);
  return prog;
}

GLuint load_vsh(const std::string & name){
  GLuint shader = glCreateShader(GL_VERTEX_SHADER);
  load_shader(name, shader);
  return shader;
}


GLuint load_fsh(const std::string & name){
  GLuint shader = glCreateShader(GL_FRAGMENT_SHADER);
  load_shader(name, shader);
  return shader;
}

GLuint load_program(GLuint & vertShader,
		    GLuint & fragShader,
		    const std::string & vertName,
		    const std::string & fragName){
  GLuint prog = glCreateProgram();

  if(vertName.size())vertShader = load_vsh(vertName);
  
  if(fragName.size())fragShader = load_fsh(fragName);
  
  link_program(prog, vertShader, fragShader);
  return prog;
}
