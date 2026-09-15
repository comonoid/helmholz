/* pview — быстрая картинка-показ: комната, зеркальные и диффузные шары.
 * Whitted: первичные лучи, зеркальные отражения, точечный свет с тенями.
 * Аналитическая сцена (плоскости+сферы), без меша. Вывод: img/room_view.ppm */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { double c[3], r; int mir; double col[3]; } S;
static const S SP[] = {
  {{-1.1,0.9,0.3},0.9,1,{0.95,0.95,1.0}},   /* зеркало */
  {{ 1.2,0.8,-0.9},0.8,1,{0.9,0.95,1.0}},   /* зеркало */
  {{ 0.2,0.55,1.4},0.55,0,{0.85,0.25,0.15}},/* диффуз красный */
  {{-0.3,0.45,-1.5},0.45,0,{0.15,0.35,0.9}},/* диффуз синий */
  {{ 2.3,0.35,0.6},0.35,0,{0.9,0.75,0.15}}  /* диффуз жёлтый */
};
#define NS 5
static const double L[3] = {0.0,3.05,0.0}; /* лампа под потолком */
static double iabs(const double *v){return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);}
static int hit(const double o[3],const double d[3],double tmin,double tmax,double *t,double n[3],double col[3],int *mir){
  double best=tmax; int hitid=-1;
  for(int i=0;i<NS;i++){double oc[3]={o[0]-SP[i].c[0],o[1]-SP[i].c[1],o[2]-SP[i].c[2]};
    double b=oc[0]*d[0]+oc[1]*d[1]+oc[2]*d[2], c=oc[0]*oc[0]+oc[1]*oc[1]+oc[2]*oc[2]-SP[i].r*SP[i].r;
    double disc=b*b-c; if(disc<0)continue; double t=-b-sqrt(disc);
    if(t>tmin&&t<best){best=t;hitid=100+i;for(int a=0;a<3;a++){n[a]=(o[a]+t*d[a]-SP[i].c[a])/SP[i].r;col[a]=SP[i].col[a];}*mir=SP[i].mir;}}
  /* пол y=0, потолок y=3.1, стены x=±4, z=±3.5 */
  if(fabs(d[1])>1e-9){double t=(0.0-o[1])/d[1]; if(t>tmin&&t<best){best=t;hitid=1;double nn[3]={0,1,0};for(int a=0;a<3;a++)n[a]=nn[a];col[0]=col[1]=col[2]=0.62;*mir=0;}}
  if(fabs(d[1])>1e-9){double t=(3.1-o[1])/d[1]; if(t>tmin&&t<best){best=t;hitid=2;double nn[3]={0,-1,0};for(int a=0;a<3;a++)n[a]=nn[a];col[0]=col[1]=col[2]=0.5;*mir=0;}}
  if(fabs(d[0])>1e-9){double t=(4.0-o[0])/d[0]; if(t>tmin&&t<best){best=t;hitid=3;double nn[3]={-1,0,0};for(int a=0;a<3;a++)n[a]=nn[a];col[0]=col[1]=col[2]=0.6;*mir=0;}}
  if(fabs(d[0])>1e-9){double t=(-4.0-o[0])/d[0]; if(t>tmin&&t<best){best=t;hitid=4;double nn[3]={1,0,0};for(int a=0;a<3;a++)n[a]=nn[a];col[0]=col[1]=col[2]=0.6;*mir=0;}}
  if(fabs(d[2])>1e-9){double t=(3.5-o[2])/d[2]; if(t>tmin&&t<best){best=t;hitid=5;double nn[3]={0,0,-1};for(int a=0;a<3;a++)n[a]=nn[a];col[0]=col[1]=col[2]=0.6;*mir=0;}}
  if(fabs(d[2])>1e-9){double t=(-3.5-o[2])/d[2]; if(t>tmin&&t<best){best=t;hitid=6;double nn[3]={0,0,1};for(int a=0;a<3;a++)n[a]=nn[a];col[0]=col[1]=col[2]=0.6;*mir=0;}}
  *t=best; return hitid;
}
static double shadow(const double p[3]){double d[3]={L[0]-p[0],L[1]-p[1],L[2]-p[2]};double dist=iabs(d);for(int a=0;a<3;a++)d[a]/=dist;
  for(int i=0;i<NS;i++){double oc[3]={p[0]-SP[i].c[0],p[1]-SP[i].c[1],p[2]-SP[i].c[2]};
    double b=oc[0]*d[0]+oc[1]*d[1]+oc[2]*d[2],c=oc[0]*oc[0]+oc[1]*oc[1]+oc[2]*oc[2]-SP[i].r*SP[i].r;
    double disc=b*b-c; if(disc<0)continue; double t=-b-sqrt(disc); if(t>1e-4&&t<dist) return 0.0;}
  return 1.0;}
static void trace(const double o[3],const double d[3],int depth,double out[3]){
  double t,n[3],col[3];int mir;int h=hit(o,d,1e-4,1e30,&t,n,col,&mir);
  if(h<0){out[0]=out[1]=out[2]=0.02;return;}
  double p[3]={o[0]+t*d[0],o[1]+t*d[1],o[2]+t*d[2]};
  double ld[3]={L[0]-p[0],L[1]-p[1],L[2]-p[2]};double dist=iabs(ld);for(int a=0;a<3;a++)ld[a]/=dist;
  double ndl=n[0]*ld[0]+n[1]*ld[1]+n[2]*ld[2]; if(ndl<0)ndl=0;
  double vis=shadow(p);
  double att=110.0/(dist*dist+1.0);
  double e[3]={col[0]*ndl*vis*att+0.06*n[1]>0?col[0]*(ndl*vis*att+0.14*(n[1]*0.5+0.5)):0,
               col[1]*ndl*vis*att+0.14*(n[1]*0.5+0.5),
               col[2]*ndl*vis*att+0.14*(n[1]*0.5+0.5)};
  out[0]=e[0];out[1]=e[1];out[2]=e[2];
  /* лампа — светящаяся сфера, видна напрямую */
  {double oc[3]={o[0]-L[0],o[1]-L[1],o[2]-L[2]};double b=oc[0]*d[0]+oc[1]*d[1]+oc[2]*d[2],c=oc[0]*oc[0]+oc[1]*oc[1]+oc[2]*oc[2]-0.09;
   double disc=b*b-c; if(disc>=0){double tt=-b-sqrt(disc); if(tt>1e-4&&tt<t){out[0]=out[1]=out[2]=40.0;return;}}}
  if(mir&&depth<4){double r[3]={d[0]-2*(d[0]*n[0]),d[1]-2*(d[1]*n[1]),d[2]-2*(d[2]*n[2])};
    double o2[3]={p[0]+1e-4*r[0],p[1]+1e-4*r[1],p[2]+1e-4*r[2]};double rc[3];
    trace(o2,r,depth+1,rc);out[0]+=0.85*rc[0];out[1]+=0.85*rc[1];out[2]+=0.85*rc[2];}
}
int main(void){
  int W=1920,H=1080;double cam[3]={-3.7,2.6,-3.4},la[3]={0.4,0.6,0.2};
  double f[3]={la[0]-cam[0],la[1]-cam[1],la[2]-cam[2]};double fl=iabs(f);for(int a=0;a<3;a++)f[a]/=fl;
  double up[3]={0,1,0},r[3]={f[1]*up[2]-f[2]*up[1],f[2]*up[0]-f[0]*up[2],f[0]*up[1]-f[1]*up[0]};
  double rl=iabs(r);for(int a=0;a<3;a++)r[a]/=rl;
  double u[3]={r[1]*f[2]-r[2]*f[1],r[2]*f[0]-r[0]*f[2],r[0]*f[1]-r[1]*f[0]};
  double tanh_=0.75;
  FILE*fp=fopen("img/room_view.ppm","wb");
  fprintf(fp,"P6\n%d %d\n255\n",W,H);
  unsigned char*row=malloc((size_t)W*3);
  for(int j=0;j<H;j++){for(int i=0;i<W;i++){
    double px=(2.0*(i+0.5)/W-1.0)*1.0667*tanh_, py=(1.0-2.0*(j+0.5)/H)*tanh_;
    double d[3]={f[0]+px*r[0]+py*u[0],f[1]+px*r[1]+py*u[1],f[2]+px*r[2]+py*u[2]};
    double dl=iabs(d);for(int a=0;a<3;a++)d[a]/=dl;
    double c[3];trace(cam,d,0,c);
    for(int a=0;a<3;a++){double v=c[a]/(1.0+c[a]);int q=(int)(255.0*pow(v,1.0/2.2));row[3*i+a]=q<0?0:(q>255?255:q);}
  }fwrite(row,1,(size_t)W*3,fp);}
  fclose(fp);free(row);return 0;
}
