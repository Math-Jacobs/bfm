#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "mex.h"





double compute_l2_ot(double *mu, double *nu, double *phi, double *dual, double totalMass, double sigma, int maxIters, int n1, int n2);



void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    
    
    double *mu=mxGetPr(prhs[0]);
    double *nu=mxGetPr(prhs[1]);
    int maxIters=(int) mxGetScalar(prhs[2]);
    double sigma =(double) mxGetScalar(prhs[3]);
    
	int n1=mxGetM(prhs[0]);
	int n2=mxGetN(prhs[0]);
    
    int pcount=n1*n2;
    
    plhs[0] = mxCreateDoubleMatrix(n1,n2,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n1,n2,mxREAL);
    
    double *phi=mxGetPr(plhs[0]);
    double *psi=mxGetPr(plhs[1]);
    
        
    double sum=0;

    for(int i=0;i<pcount;i++){
		if (mu[i]<0){
			mexErrMsgTxt("Initial density contains negative values");
		}
		if (nu[i]<0){
			mexErrMsgTxt("Final density contains negative values");
		}
        
		sum+=mu[i];
    }
    
       
    
    double totalMass=sum/(n1*n2*1.0);
	   
   
    double value=compute_l2_ot(mu, nu, phi, psi, totalMass, sigma, maxIters, n1, n2);
    
    

    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);
            
            phi[i*n1+j]=.5*(x*x+y*y)-phi[i*n1+j];
            psi[i*n1+j]=.5*(x*x+y*y)-psi[i*n1+j];            
        }
    }
    
}






typedef struct{
    fftw_plan dctIn;
    fftw_plan dctOut;
    double *kernel;
    double *workspace;
}poisson_solver;


double *create_negative_laplace_kernel(int n1, int n2){
    double *kernel=calloc(n1*n2,sizeof(double));
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            double x=M_PI*j/(n1*1.0);
            double y=M_PI*i/(n2*1.0);
            
            double negativeLaplacian=2*n1*n1*(1-cos(x))+2*n2*n2*(1-cos(y));
            
            kernel[i*n1+j]=negativeLaplacian;
                
            
            
        }
    }
    return kernel;
}


poisson_solver create_poisson_solver_workspace(int n1, int n2){
    clock_t b,e;
    b=clock();
    poisson_solver fftps;
    fftps.workspace=calloc(n1*n2,sizeof(double));
    fftps.kernel=create_negative_laplace_kernel(n1,n2);
    
    fftps.dctIn=fftw_plan_r2r_2d(n2, n1, fftps.workspace, fftps.workspace,
                                 FFTW_REDFT10, FFTW_REDFT10,
                                 FFTW_MEASURE);
    fftps.dctOut=fftw_plan_r2r_2d(n2, n1, fftps.workspace, fftps.workspace,
                                  FFTW_REDFT01, FFTW_REDFT01,
                                  FFTW_MEASURE);
    
    e=clock();
    
    mexPrintf("FFT setup time: %.2fs\n", (e-b)/(CLOCKS_PER_SEC*1.0));
    mexEvalString("pause(.001);");

    return fftps;
}



void destroy_poisson_solver(poisson_solver fftps){
    free(fftps.kernel);
    free(fftps.workspace);
    fftw_destroy_plan(fftps.dctIn);
    fftw_destroy_plan(fftps.dctOut);
}




typedef struct{
    int *indices;
    int hullCount;
    
}convex_hull;


int sgn(double x){
    
    int truth=(x>0)-(x<0);
    return truth;
    
}


void init_hull(convex_hull *hull, int n){
    hull->indices=calloc(n,sizeof(double));
    hull->hullCount=0;
    
}

void destroy_hull(convex_hull *hull){
    free(hull->indices);
}

void transpose_doubles(double *transpose, double *data, int n1, int n2){
    
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
         
            transpose[j*n2+i]=data[i*n1+j];
        }
    }
}




double interpolate_function(double *function, double x, double y, int n1, int n2){
    
    int xIndex=fmin(fmax(x*n1-.5 ,0),n1-1);
    int yIndex=fmin(fmax(y*n2-.5 ,0),n2-1);
    
    double xfrac=x*n1-xIndex-.5;
    double yfrac=y*n2-yIndex-.5;
    
    int xOther=xIndex+sgn(xfrac);
    int yOther=yIndex+sgn(yfrac);
    
    xOther=fmax(fmin(xOther, n1-1),0);
    yOther=fmax(fmin(yOther, n2-1),0);
    
    double v1=(1-fabs(xfrac))*(1-fabs(yfrac))*function[yIndex*n1+xIndex];
    double v2=fabs(xfrac)*(1-fabs(yfrac))*function[yIndex*n1+xOther];
    double v3=(1-fabs(xfrac))*fabs(yfrac)*function[yOther*n1+xIndex];
    double v4=fabs(xfrac)*fabs(yfrac)*function[yOther*n1+xOther];
    
    double v=v1+v2+v3+v4;
    
    return v;
    
}


void add_point(double *u, convex_hull *hull, int i){
    
    
    if(hull->hullCount<2){
        hull->indices[1]=i;
        hull->hullCount++;
    }else{
        int hc=hull->hullCount;
        int ic1=hull->indices[hc-1];
        int ic2=hull->indices[hc-2];
        
        double oldSlope=(u[ic1]-u[ic2])/(ic1-ic2);
        double slope=(u[i]-u[ic1])/(i-ic1);
        
        if(slope>=oldSlope){
            int hc=hull->hullCount;
            hull->indices[hc]=i;
            hull->hullCount++;
        }else{
            hull->hullCount--;
            add_point(u, hull, i);
        }
    }
}


void get_convex_hull(double *u, convex_hull *hull, int n){
    
    hull->indices[0]=0;
    hull->indices[1]=1;
    hull->hullCount=2;
    
    for(int i=2;i<n;i++){
        
        add_point(u, hull, i);
        
    }
    
}


void compute_dual_indices(int *dualIndicies, double *u, convex_hull *hull, int n){
    
    int counter=1;
    int hc=hull->hullCount;
    
    for(int i=0;i<n;i++){
       
        double s=(i+.5)/(n*1.0);
        int ic1=hull->indices[counter];
        int ic2=hull->indices[counter-1];
        
        double slope=n*(u[ic1]-u[ic2])/(ic1-ic2);
        while(s>slope&&counter<hc-1){
            counter++;
            ic1=hull->indices[counter];
            ic2=hull->indices[counter-1];
            slope=n*(u[ic1]-u[ic2])/(ic1-ic2);
        }
        dualIndicies[i]=hull->indices[counter-1];
        
    }
}


void compute_dual(double *dual, double *u, int *dualIndicies, convex_hull *hull, int n){
    
    get_convex_hull(u, hull, n);
   
    
    compute_dual_indices(dualIndicies, u, hull, n);
    
    for(int i=0;i<n;i++){
        double s=(i+.5)/(n*1.0);
        int index=dualIndicies[i];
        double x=(index+.5)/(n*1.0);
        double v1=s*x-u[dualIndicies[i]];
        double v2=s*(n-.5)/(n*1.0)-u[n-1];
        if(v1>v2){
            dual[i]=v1;
        }else{
            dualIndicies[i]=n-1;
            dual[i]=v2;
        }
        
    }
    
}




void compute_2d_dual(double *dual, double *u, convex_hull *hull, int n1, int n2){
    
    int pcount=n1*n2;
    
    int n=fmax(n1,n2);
    
    int *argmin=calloc(n,sizeof(int));
    
    double *temp=calloc(pcount,sizeof(double));
    
    memcpy(temp, u, pcount*sizeof(double));
    
    
    for(int i=0;i<n2;i++){
        compute_dual(&dual[i*n1], &temp[i*n1], argmin, hull, n1);
        
    }
    transpose_doubles(temp, dual, n1, n2);
    for(int i=0;i<n1*n2;i++){
        dual[i]=-temp[i];
    }
    for(int j=0;j<n1;j++){
        compute_dual(&temp[j*n2], &dual[j*n2], argmin, hull, n2);
        
    }
    transpose_doubles(dual, temp, n2, n1);
    
    free(temp);
    free(argmin);
    
}



void convexify(double *phi, double *dual, convex_hull *hull, int n1, int n2){
    
    compute_2d_dual(dual, phi, hull, n1, n2);
    
    compute_2d_dual(phi, dual, hull, n1, n2);
    
}





void calc_pushforward_map(double *xMap, double *yMap, double *dual, int n1, int n2){
    
    
    double xStep=1.0/n1;
    double yStep=1.0/n2;
    
    
    for(int i=0;i<n2+1;i++){
        for(int j=0;j<n1+1;j++){
            double x=j/(n1*1.0);
            double y=i/(n2*1.0);
            
            double dualxp=interpolate_function(dual, x+xStep, y, n1, n2);
            double dualxm=interpolate_function(dual, x-xStep, y, n1, n2);
            
            double dualyp=interpolate_function(dual, x, y+yStep, n1, n2);
            double dualym=interpolate_function(dual, x, y-yStep, n1, n2);
            
            xMap[i*(n1+1)+j]=.5*n1*(dualxp-dualxm);
            yMap[i*(n1+1)+j]=.5*n2*(dualyp-dualym);
            
            
        }
    }
    
}





void sampling_pushforward(double *rho, double *mu, double totalMass, double *xMap, double *yMap, int n1, int n2){
    
    int pcount=n1*n2;
    
    memset(rho,0,pcount*sizeof(double));
    
    
    double xCut=pow(1.0/n1,1.0/3);
    double yCut=pow(1.0/n2,1.0/3);
    
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            
            double mass=mu[i*n1+j];
            
            if(mass>0){
                
                double xStretch0=fabs(xMap[i*(n1+1)+j+1]-xMap[i*(n1+1)+j]);
                double xStretch1=fabs(xMap[(i+1)*(n1+1)+j+1]-xMap[(i+1)*(n1+1)+j]);
                
                double yStretch0=fabs(yMap[(i+1)*(n1+1)+j]-yMap[i*(n1+1)+j]);
                double yStretch1=fabs(yMap[(i+1)*(n1+1)+j+1]-yMap[i*(n1+1)+j+1]);
                
                double xStretch=fmax(xStretch0, xStretch1);
                double yStretch=fmax(yStretch0, yStretch1);
                
                int xSamples=2*fmax(n1*xStretch,1);
                int ySamples=2*fmax(n2*yStretch,1);
                
                if(xStretch<xCut&&yStretch<yCut){
                    
                    double factor=1/(xSamples*ySamples*1.0);
                    
                    for(int l=0;l<ySamples;l++){
                        for(int k=0;k<xSamples;k++){
                            
                            double a=(k+.5)/(xSamples*1.0);
                            double b=(l+.5)/(ySamples*1.0);
                            
                            double xPoint=(1-b)*(1-a)*xMap[i*(n1+1)+j]+(1-b)*a*xMap[i*(n1+1)+j+1]+b*(1-a)*xMap[(i+1)*(n1+1)+j]+a*b*xMap[i*(n1+1)+j+1];
                            double yPoint=(1-b)*(1-a)*yMap[i*(n1+1)+j]+(1-b)*a*yMap[i*(n1+1)+j+1]+b*(1-a)*yMap[(i+1)*(n1+1)+j]+a*b*yMap[i*(n1+1)+j+1];
                            
                            double X=xPoint*n1-.5;
                            double Y=yPoint*n2-.5;
                            
                            int xIndex=X;
                            int yIndex=Y;
                            
                            double xFrac=X-xIndex;
                            double yFrac=Y-yIndex;
                            
                            int xOther=xIndex+1;
                            int yOther=yIndex+1;
                            
                            xIndex=fmin(fmax(xIndex,0),n1-1);
                            xOther=fmin(fmax(xOther,0),n1-1);
                            
                            yIndex=fmin(fmax(yIndex,0),n2-1);
                            yOther=fmin(fmax(yOther,0),n2-1);
                            
                            
                            rho[yIndex*n1+xIndex]+=(1-xFrac)*(1-yFrac)*mass*factor;
                            rho[yOther*n1+xIndex]+=(1-xFrac)*yFrac*mass*factor;
                            rho[yIndex*n1+xOther]+=xFrac*(1-yFrac)*mass*factor;
                            rho[yOther*n1+xOther]+=xFrac*yFrac*mass*factor;
                            
                        }
                    }
                }
                
            }
            
        }
    }
    
    double sum=0;
    for(int i=0;i<pcount;i++){
        sum+=rho[i]/pcount;
    }
    for(int i=0;i<pcount;i++){
        rho[i]*=totalMass/sum;
    }
    
}







double update_potential(poisson_solver fftps, double *phi, double *rho, double *nu, double sigma, int n1, int n2){
    
    int pcount=n1*n2;
    
    double h1=0;
    
    for(int i=0;i<pcount;i++){
        fftps.workspace[i]=rho[i]-nu[i];
    }
    
    fftw_execute(fftps.dctIn);
    
    fftps.workspace[0]=0;
    
    for(int i=1;i<pcount;i++){
        
        fftps.workspace[i]/=4*pcount*fftps.kernel[i];
        
    }
   
    
    fftw_execute(fftps.dctOut);
    
    for(int i=0;i<pcount;i++){
        phi[i]+=sigma*fftps.workspace[i];
        h1+=fftps.workspace[i]*(rho[i]-nu[i]);
    }
    
    h1/=pcount;
    
    
    
    return h1;
    
}





double compute_w2(double *phi, double *dual, double *mu, double *nu, int n1, int n2){
    
    int pcount=n1*n2;
    
    double value=0;
    
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);
            
            value+=.5*(x*x+y*y)*(mu[i*n1+j]+nu[i*n1+j])-nu[i*n1+j]*phi[i*n1+j]-mu[i*n1+j]*dual[i*n1+j];
        }
    }
    
    value/=pcount;
    
    return value;
    
}



double step_update(double sigma, double value, double oldValue, double gradSq, double scaleUp, double scaleDown, double upper, double lower){
    
    double diff=value-oldValue;
    
    if(diff>gradSq*sigma*upper){
        return sigma*scaleUp;
    }else if(diff<gradSq*sigma*lower){
        return sigma*scaleDown;
    }else{
        return sigma;
    }
    
}

double compute_l2_ot(double *mu, double *nu, double *phi, double *dual, double totalMass, double sigma, int maxIters, int n1, int n2){
    
    int pcount=n1*n2;
    poisson_solver fftps=create_poisson_solver_workspace(n1,n2);
    

    
    double *xMap=calloc((n1+1)*(n2+1),sizeof(double));
    double *yMap=calloc((n1+1)*(n2+1),sizeof(double));
    
   
    
    int n=fmax(n1,n2);
    
    convex_hull hull;
    
    init_hull(&hull, n);
    
    
    
   
    
    for(int i=0;i<n2+1;i++){
        for(int j=0;j<n1+1;j++){
            
            double x=j/(n1*1.0);
            double y=i/(n2*1.0);
            
            xMap[i*n1+j]=x;
            yMap[i*n1+j]=y;
            
        }
    }
    
    for(int i=0;i<n2;i++){
        for(int j=0;j<n1;j++){
            double x=(j+.5)/(n1*1.0);
            double y=(i+.5)/(n2*1.0);
            
            phi[i*n1+j]=.5*(x*x+y*y);
            dual[i*n1+j]=.5*(x*x+y*y);
        }
    }
       
    
    double *rho=calloc(pcount,sizeof(double));
    memcpy(rho,mu,pcount*sizeof(double));
    
    
    double oldValue=compute_w2(phi, dual, mu, nu, n1, n2);
    
    double scaleDown=.8;
    double scaleUp=1/scaleDown;
    double upper=.75;
    double lower=.25;

    int numDigitsIter=floor(log10(maxIters) + 1);
    
    for(int i=0;i<maxIters+1;i++){
        
        
        double gradSq=update_potential(fftps, phi, rho, nu, sigma, n1, n2);
        
       
        
        convexify(phi, dual, &hull, n1, n2);
        
        
        
        double value=compute_w2(phi, dual, mu, nu, n1, n2);
        
        sigma=step_update(sigma, value, oldValue, gradSq , scaleUp, scaleDown, upper, lower);
        
        oldValue=value;
        
        
        calc_pushforward_map(xMap, yMap, phi, n1, n2);
        
        sampling_pushforward(rho, nu, totalMass, xMap, yMap, n1, n2);
        
        
        
        
        gradSq=update_potential(fftps, dual, rho, mu, sigma, n1, n2);
        
        
        
        convexify(dual, phi, &hull, n1, n2);
        
        calc_pushforward_map(xMap, yMap, dual, n1, n2);
        
        sampling_pushforward(rho, mu, totalMass, xMap, yMap, n1, n2);
        
        
        
        value=compute_w2(phi, dual, mu, nu, n1, n2);
        
        sigma=step_update(sigma, value, oldValue, gradSq , scaleUp, scaleDown, upper, lower);
        
        oldValue=value;
        
        
        sigma=fmax(sigma,.05);
        
        if(i%5==0){
                                    
            mexPrintf("iter %*d, W2 value: %5e\n", numDigitsIter, i, value);
            mexEvalString("pause(.001);");
            
        }
        
       
        
       
        
        
    }
    
    
    destroy_hull(&hull);
    free(rho);
    free(xMap);
    free(yMap);
    destroy_poisson_solver(fftps);
    
    return oldValue;
}




