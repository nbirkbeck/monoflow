#include <nimage/image.h>
#include <nmath/vec3.h>
#include <assert.h>
#include <nmath/matrix.h>

using namespace nacb;

static double sqr(double x){
  return x*x;
}

Image8 bilateralSilEdges(Image8 & image, int hwin){
  if(image.nchannels < 4)return image.copy();

  double sigmap = double(hwin)/1.7;
  double sigmaa = 20.0;
  Image8 result;
  result = image.copy();
  
  for(int y=hwin; y<result.h-hwin; y++){
    for(int x=hwin; x<result.w-hwin; x++){
      if(image(x,y,3)){
	double wtot = 0;
	double tot = 0;
	Vec3f col(image(x,y,0), image(x,y,1), image(x,y,2));
	
	for(int yy=y-hwin; yy<=y+hwin; yy++){
	  for(int xx=x-hwin; xx<=x+hwin; xx++){
	    double wx = exp(-sqr(xx-x)/(2.0*sigmap*sigmap));
	    double wy = exp(-sqr(yy-y)/(2.0*sigmap*sigmap));
	    Vec3f c(image(xx,yy,0), image(xx,yy,1), image(xx,yy,2));
	    double w = exp(-sqr((c-col).len())/(2.0*sigmaa*sigmaa))*wx*wy;
	    tot  += double(image(xx,yy,3))*w;
	    wtot += w;
	  }
	}
	result(x,y,3) = (unsigned char)std::max(0.0, std::min(255.0, tot/wtot)); 
      }
    }
  }
  return result;
}



template <class T>
Imagef grad(Image<T> & im,int chan){
  Imagef g(im.w,im.h,2);
  for(int y=1;y<im.h-1;y++){
    for(int x=1;x<im.w-1;x++){
      g(x,y,0)=(im(x+1,y,chan)-im(x-1,y,chan))/2.0;
      g(x,y,1)=(im(x,y+1,chan)-im(x,y-1,chan))/2.0;
    }
  }
  for(int y=0;y<im.h;y++){
    g(0,y,0)=im(1,y,chan)-im(0,y,chan);
    g(im.w-1,y,0)=im(im.w-1,y,chan)-im(im.w-2,y,chan);

    if(y==0){
      g(0,0,1)=im(0,1,chan)-im(0,0,chan);
      g(im.w-1,0,1)=im(im.w-1,1,chan)-im(im.w-1,0,chan);
    }
    else if(y==im.h-1){
      g(0,y,1)=im(0,im.h-1,chan)-im(0,im.h-2,chan);
      g(im.w-1,y,1)=im(im.w-1,im.h-1,chan)-im(im.w-1,im.h-2,chan);
    }
    else {
      g(0,y,1)=(im(0,y+1,chan)-im(0,y-1,chan))/2.0;
      g(im.w-1,y,1)=(im(im.w-1,y+1,chan)-im(im.w-1,y-1,chan))/2.0;
    }
  }
  for(int x=0;x<im.w;x++){
    g(x,0,1)=(im(x,1,chan)-im(x,0,chan));
    g(x,im.h-1,1)=(im(x,im.h-1,chan)-im(x,im.h-2,chan));

    if(x==0){
      g(0,0,0)=im(1,0,chan)-im(0,0,chan);
      g(0,im.h-1,0)=im(1,im.h-1,chan)-im(0,im.h-1,chan);
    }
    else if(x==im.w-1){
      g(x,0,0)=im(x,0,chan)-im(x-1,0,chan);
      g(x,im.h-1,0)=im(x,im.h-1,chan)-im(x-1,im.h-1,chan);
    }
    else{
      g(x,0,0)=(im(x+1,0,chan)-im(x-1,0,chan))/2.0;
      g(x,im.h-1,0)=(im(x+1,im.h-1,chan)-im(x-1,im.h-1,chan))/2.0;
    }
  }
  return g;
}
template Imagef grad(Image8 & im,int chan);
template Imagef grad(Imagef & im,int chan);

Imagef motion_flow(Image8 & im1,Image8 & im2,int hwin,int its){
  assert(im1.w==im2.w);
  assert(im1.h==im2.h);
  assert(im1.nchannels==im2.nchannels);
  
  Imagef delta(im1.w,im1.h,2);

  int nch=std::min(3,im1.nchannels);
  Imagef grads[nch];
  for(int k=0;k<nch;k++)
    grads[k]=grad(im2,k);
  printf("got gradient\n");

  int win=2*hwin+1;
  int ncon=(win*win*nch);
  Matrix A(ncon,2),b(ncon,1);

  //I2(x+dx,y+dy)-I1(x,y)=0
  //grad(I2)[dx,dy]'=I1-I2

  Imagef weights(win,win,1);
  //double sum=0;
  double sigma=double(hwin)*2.0/3.0;
  //don't really care about overall scale of weights
  for(int y=0;y<win;y++){
    for(int x=0;x<win;x++){
      double xf=x-hwin;
      double yf=y-hwin;
      weights(x,y)=exp(-(xf*xf/(2*sigma*sigma)-yf*yf/(2*sigma*sigma)));
    }
  }
  delta=0.0f;

  for(int it=0;it<1;it++){
    for(int y=0;y<im1.h;y++){
      for(int x=0;x<im1.w;x++){
	if(x>=hwin && y>=hwin &&
	   x<im1.w-hwin &&
	   y<im1.h-hwin && 
	   (im1.nchannels<=3 || im1(x,y,3)>128)){

	  int on=0;
	  for(int k=0;k<nch;k++){
	    for(int yi=0;yi<win;yi++){
	      int y2=yi+y-hwin;
	      for(int xi=0;xi<win;xi++){
		int x2=xi+x-hwin;
		/*Vec3f g,c2;
		bilinear(grads[k],delta(x,y,0)+float(x2),delta(x,y,1)+float(y2),g);
		bilinear(im2,delta(x,y,0)+float(x2),delta(x,y,1)+float(y2),g);

		A(k*win*win+yi*win+xi,0)=g.x*weights(xi,yi);
		A(k*win*win+yi*win+xi,1)=g.y*weights(xi,yi);
		b[k*win*win+yi*win+xi]=(double(im1(x2,y2,k))-c2.data[k])*weights(xi,yi);
		*/
		A(k*win*win+yi*win+xi,0)=grads[k](x2,y2,0)*weights(xi,yi);
		A(k*win*win+yi*win+xi,1)=grads[k](x2,y2,1)*weights(xi,yi);
		b[k*win*win+yi*win+xi]=(double(im1(x2,y2,k))-double(im2(x2,y2,k)))*weights(xi,yi);
		on+=(im1.nchannels<=3)?1:(im1(x2,y2,3)>128);
	      }
	    }
	  }
	  if(2*on>=ncon){
	    Matrix upd=Matrix::LlinLeastSq(A,b);
	    delta(x,y,0)+=std::min(std::max(upd[0],-1.0),1.0);
	    delta(x,y,1)+=std::min(std::max(upd[1],-1.0),1.0);
	  }
	}
      }
    }
  }
  return delta;
}


template<class T>
Image<T> box_filt(Image<T> & im,int hs){
  Image<T> xf(im.w,im.h,1);
  Image<T> yf(im.w,im.h,1);
  int ws=2*hs+1;
  for(int y=0;y<im.h;y++){
    double sum=0;
    for(int x=0;x<ws;x++){
      sum+=im(x,y);
      if(x>=hs)xf(x-hs,y)=(T)(sum/(x+1));
    }
    for(int x=ws;x<im.w;x++){
      xf(x-hs-1,y)=(T)(sum/ws);
      sum-=im(x-ws,y);
      sum+=im(x,y);
    }
    for(int x=im.w-hs-1;x<im.w;x++){
      xf(x,y)=(T)(sum/(im.w-x+hs));
      sum-=im(x-hs-1,y);
    }
  }
  for(int x=0;x<im.w;x++){
    double sum=0;
    for(int y=0;y<ws;y++){
      sum+=im(x,y);
      if(y>=hs)yf(x,y-hs)=(T)(sum/(y+1));
    }
    for(int y=ws;y<im.h;y++){
      yf(x,y-hs-1)=(T)(sum/ws);
      sum-=xf(x,y-ws);
      sum+=xf(x,y);
    }
    for(int y=im.h-hs-1;y<im.h;y++){
      yf(x,y)=(T)(sum/(im.h-y+hs));
      sum-=xf(x,y-hs-1);
    }
  }
  return yf;
}

//Template definitions
template Image8 box_filt(Image8 & im,int hs);
template Imagef box_filt(Imagef & im,int hs);

#ifdef FILTERS_MAIN

int main(int ac, char * av[]){
  for(int i=1; i<ac; i++){
    Image8 im(av[i]);
    im = bilateralSilEdges(im, 4);
    im.write(av[i]);
  }
}

#endif
