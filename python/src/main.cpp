#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;



class convex_hull{
public:
    int* indices;
    int  hullCount;

    convex_hull(int n){
        indices=new int[n];
        hullCount=0;
    }

    ~convex_hull(){
        delete[] indices;
    }
};



// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---



class BFM{
public:
int n1;
int n2;

double totalMass;

double *xMap;
double *yMap;

double *rho;

int *argmin;
double *temp;
    
convex_hull* hull;

    BFM(int n1, int n2, py::array_t<double> & mu_np){

        py::buffer_info mu_buf = mu_np.request();
        double *mu = static_cast<double *>(mu_buf.ptr);

        this->n1 = n1;
        this->n2 = n2;

        int n=fmax(n1,n2);
        hull   = new convex_hull(n);    
        argmin = new int[n];
        temp   = new double[n1*n2];

        xMap=new double[(n1+1)*(n2+1)];
        yMap=new double[(n1+1)*(n2+1)];

        for(int i=0;i<n2+1;i++){
            for(int j=0;j<n1+1;j++){
                
                double x=j/(n1*1.0);
                double y=i/(n2*1.0);
                
                xMap[i*(n1+1)+j]=x;
                yMap[i*(n1+1)+j]=y;
                
            }
        }

        rho=new double[n1*n2];
        memcpy(rho,mu,n1*n2*sizeof(double));

        totalMass = 0;
        for(int i=0;i<n1*n2;++i){
            totalMass += mu[i];
        }
        totalMass /= n1*n2;

    }

    ~BFM(){
        delete [] xMap;
        delete [] yMap;
        delete [] rho;
        delete hull;
    }

    void ctransform(py::array_t<double> & dual_np, py::array_t<double> & phi_np){

        py::buffer_info phi_buf  = phi_np.request();
        py::buffer_info dual_buf = dual_np.request();

        double *phi  = static_cast<double *> (phi_buf.ptr);
        double *dual = static_cast<double *> (dual_buf.ptr);

        compute_2d_dual_inside(dual, phi, hull, n1, n2);
    }

    void pushforward(py::array_t<double> & rho_np, py::array_t<double> & phi_np, py::array_t<double> & nu_np){

        py::buffer_info phi_buf  = phi_np.request();
        py::buffer_info nu_buf   = nu_np.request();

        double *phi = static_cast<double *> (phi_buf.ptr);
        double *nu  = static_cast<double *> (nu_buf.ptr);

        calc_pushforward_map(phi);
        sampling_pushforward(nu);

        py::buffer_info rho_buf  = rho_np.request();
        memcpy(static_cast<double *> (rho_buf.ptr), rho, n1*n2*sizeof(double));
    }


    double compute_w2(py::array_t<double> & phi_np, py::array_t<double> & dual_np, py::array_t<double> & mu_np, py::array_t<double> & nu_np){

        py::buffer_info phi_buf  = phi_np.request();
        py::buffer_info dual_buf = dual_np.request();
        py::buffer_info mu_buf   = mu_np.request();
        py::buffer_info nu_buf   = nu_np.request();

        double *phi  = static_cast<double *> (phi_buf.ptr);
        double *dual = static_cast<double *> (dual_buf.ptr);
        double *mu   = static_cast<double *> (mu_buf.ptr);
        double *nu   = static_cast<double *> (nu_buf.ptr);
        
        int pcount=n1*n2;
        
        double value=0;
        
        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){
                double x=(j+.5)/(n1*1.0);
                double y=(i+.5)/(n2*1.0);
                
                value+=.5*(x*x+y*y)*(mu[i*n1+j]+nu[i*n1+j])-nu[i*n1+j]*phi[i*n1+j]-mu[i*n1+j]*dual[i*n1+j];
            }
        }

        printf("value calculated");

        value/=pcount;
        
        return value;
    }

    void compute_2d_dual_inside(double *dual, double *u, convex_hull *hull, int n1, int n2){
        
        int pcount=n1*n2;
        
        int n=fmax(n1,n2);
        
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


    int sgn(double x){
        
        int truth=(x>0)-(x<0);
        return truth;
        
    }


    void transpose_doubles(double *transpose, double *data, int n1, int n2){
        
        for(int i=0;i<n2;i++){
            for(int j=0;j<n1;j++){
             
                transpose[j*n2+i]=data[i*n1+j];
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



    void calc_pushforward_map(double *dual){
        
        
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



    void sampling_pushforward(double *mu){


        int pcount=n1*n2;
        
        memset(rho,0,pcount*sizeof(double));
        
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
                    
                    int xSamples=fmax(n1*xStretch,1);
                    int ySamples=fmax(n2*yStretch,1);

                    double factor=1/(xSamples*ySamples*1.0);
                    
                    for(int l=0;l<ySamples;l++){
                        for(int k=0;k<xSamples;k++){
                            
                            double a=(k+.5)/(xSamples*1.0);
                            double b=(l+.5)/(ySamples*1.0);
                            
                            double xPoint=(1-b)*(1-a)*xMap[i*(n1+1)+j]+(1-b)*a*xMap[i*(n1+1)+j+1]+b*(1-a)*xMap[(i+1)*(n1+1)+j]+a*b*xMap[(i+1)*(n1+1)+(j+1)];
                            double yPoint=(1-b)*(1-a)*yMap[i*(n1+1)+j]+(1-b)*a*yMap[i*(n1+1)+j+1]+b*(1-a)*yMap[(i+1)*(n1+1)+j]+a*b*yMap[(i+1)*(n1+1)+(j+1)];
                            
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
        
        double sum=0;
        for(int i=0;i<pcount;i++){
            sum+=rho[i]/pcount;
        }
        for(int i=0;i<pcount;i++){
            rho[i]*=totalMass/sum;
        }
    }

};


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


PYBIND11_MODULE(w2, m) {
    // optional module docstring
    m.doc() = "pybind11 for w2 code";

    py::class_<BFM>(m, "BFM")
        .def(py::init<int, int, py::array_t<double> &>())
        .def("ctransform", &BFM::ctransform)
        .def("pushforward", &BFM::pushforward);
}
