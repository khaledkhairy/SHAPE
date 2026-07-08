#pragma once
#include "utility.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include <math.h>
using namespace std;
//#define TRUE 1
#define YES 1
//#define FALSE 0
#define NO 0
#define KKDEBUG 1



class spherical_mesh
{
public:
    // variables applicable to any triangular mesh
    int tri_n;	// index of the spherical triangulation to use
    int n_points, n_faces;
    int *f0, *f1, *f2;
    double* x1, *y1, *z1;
    double* x2, *y2, *z2;
    double* x3, *y3, *z3, *twoA;
    double **X;				// actual surface 3-vector columns come in here
    double **q, **r, **crossqpr, **n;
    double *H, *dA;
    double A, V, v, Eb, h;
    double *sf;		// store the scalar field n_points x 1 array
    std::vector< std::vector<double> > vN;	// normal at vertex
    
    //member variables specific for the mesh being spherical at initialization and being used in conjuction with shp_surface
    int L_max;
    MatrixXd p,t;	// column matrices
    MatrixXd YLK,Y_P, Y_T, Y_TT, Y_PP, Y_TP;	// basis matrix
    MatrixXd Xu, Xt, Xp, Xtt, Xpp, Xtp, E, F, G, L, M, N, nu,SSn, dXn, Hu, KGu;
    double *Han, *KGan;	// local mean curvature and local Gaussian curvature --- analytically evaluated
    std::string fn_str;	//filename with mesh information
    vertex_info* _vip;	//pointer to array of vertex_info structs;
    bool _curv_calc;
    // methods
    
    void mesh_init()		// reading edge information for curvature calculation from file is specific to this class
    {
        // initialize arrays necessary for analytical calculation of derivatives
        Xu = MatrixXd(this->n_points, 3);	// 3-vector of columns
        Xt = MatrixXd(this->n_points, 3);	// 3-vector of columns
        Xp = MatrixXd(this->n_points, 3);	// 3-vector of columns
        Xtt = MatrixXd(this->n_points, 3);	// 3-vector of columns
        Xpp = MatrixXd(this->n_points, 3);	// 3-vector of columns
        Xtp = MatrixXd(this->n_points, 3);	// 3-vector of columns
        
        E = MatrixXd(this->n_points, 1);	// 1-vector
        F = MatrixXd(this->n_points, 1);	// 1-vector
        G = MatrixXd(this->n_points, 1);	// 1-vector
        L = MatrixXd(this->n_points, 1);	// 1-vector
        M = MatrixXd(this->n_points, 1);	// 1-vector
        N = MatrixXd(this->n_points, 1);	// 1-vector
        
        SSn = MatrixXd(this->n_points, 1);	// 1-vector
        nu = MatrixXd(this->n_points, 3);	// 3-vector of columns resulting from the cross product of Xt and Xp
        
        Hu = MatrixXd(this->n_points, 1);	// 1-vector
        KGu = MatrixXd(this->n_points, 1);	// 1-vector
        this->Han = (double*)calloc(this->n_points,sizeof(double));
        this->KGan = (double*)calloc(this->n_points,sizeof(double));
        vN.resize(n_points);for(int i = 0;i<n_points;i++){vN[i].resize(3);vN[i][0] = 0.0;vN[i][1] = 0.0;vN[i][2] = 0.0;}
        // initialize tiangular mesh arrays
        //double A = 0.0;		// total surface area
        //double V = 0.0;		// total volume
        //double Vo= 0.0;		// volume of sphere of same surface area
        
        this->sf = (double*)calloc(this->n_points,sizeof(double));
        this->X = (double**)calloc(this->n_points,sizeof(double*));
        for(int i=0;i<this->n_points;i++){X[i] = (double*)calloc(3,sizeof(double));}
        ///////////////////////////////////////
        this->x1 = (double*)calloc(this->n_faces,sizeof(double));
        this->y1 = (double*)calloc(this->n_faces,sizeof(double));
        this->z1 = (double*)calloc(this->n_faces,sizeof(double));
        this->x2 = (double*)calloc(this->n_faces,sizeof(double));
        this->y2 = (double*)calloc(this->n_faces,sizeof(double));
        this->z2 = (double*)calloc(this->n_faces,sizeof(double));
        this->x3 = (double*)calloc(this->n_faces,sizeof(double));
        this->y3 = (double*)calloc(this->n_faces,sizeof(double));
        this->z3 = (double*)calloc(this->n_faces,sizeof(double));
        this->twoA = (double*)calloc(this->n_faces,sizeof(double));	//double the surface area of each triangle
        
        this->q = (double**)calloc(this->n_faces,sizeof(double*));	for(int i=0;i<this->n_faces;i++){this->q[i] = (double*)calloc(3,sizeof(double));}
        this->r = (double**)calloc(this->n_faces,sizeof(double*));	for(int i=0;i<this->n_faces;i++){this->r[i] = (double*)calloc(3,sizeof(double));}
        this->crossqpr = (double**)calloc(this->n_faces,sizeof(double*));	for(int i=0;i<this->n_faces;i++){this->crossqpr[i] = (double*)calloc(3,sizeof(double));}
        this->n = (double**)calloc(this->n_faces,sizeof(double*));	for(int i=0;i<this->n_faces;i++){this->n[i] = (double*)calloc(3,sizeof(double));}
        
        /// initialize arrays specific to curvature and bending energy calculation
        this->H = (double*)calloc(this->n_points,sizeof(double));
        this->dA = (double*)calloc(this->n_points,sizeof(double));
        
        /// prepare intermediate arrays of indices needed for mesh-based
        /// curvature calculation (read from disk -- this is messy)
        _vip = (struct vertex_info *)(calloc(this->n_points,sizeof(struct vertex_info)));	//allocate memory for the vertex infos
        int iReturn;
        int *vr;		// indices of faces a verted is member of
        int *indx;		// indices of other vertices a vertex is connected to
        int *r1;
        int *r2;
        int *r3;
        int vn;			//number of other vertices a vertex is connected to
        FILE* f;
        char fullFileName[500];
        int val2;
        sprintf(fullFileName,"./%s",fn_str.c_str());//
        f = fopen(fullFileName,"r");//
        if(f!=NULL)
        {
            for (int ix=0;ix<this->n_points;ix++)
            {	//loop over the number of vertices
                iReturn = fscanf(f,"%d",&vn);//read in the number if edge/triangle members of the vertex
                vr = (int*)(calloc(vn,sizeof(int)));
                r1 = (int*)(calloc(vn,sizeof(int)));
                r2 = (int*)(calloc(vn,sizeof(int)));
                r3 = (int*)(calloc(vn,sizeof(int)));
                indx = (int*)(calloc(vn,sizeof(int)));
                for(int rix = 0;rix<vn;rix++)
                {
                    iReturn = fscanf(f,"%d\t",&val2);//read in the values for r
                    vr[rix] = val2-1;
                }
                
                for(int rix = 0;rix<vn;rix++)
                {
                    iReturn = fscanf(f,"%d\t%d\t%d\t%d",&(indx[rix]),&(r1[rix]), &(r2[rix]), &(r3[rix]));//read in the values for Ve r1 r2
                    indx[rix]-=1;r1[rix]-=1;r2[rix]-=1;r3[rix]-=1;
                }
                _vip[ix] = vertex_info(vn,vr,indx,r1,r2, r3);
            }
            fclose(f);
            _curv_calc = true;
        }else
        {
            if(KKDEBUG){std::cout<<"Could not initialize curvature information for curvature calculation from triangular mesh"<<std::endl;}
            _curv_calc = false;
        }
    }
    
    void sfGen(const MatrixXd &xc) //// generate surface field based on the spherical harmonics coefficients in xc
    {
        //int nc = (this->L_max + 1) *(this->L_max + 1);
        MatrixXd tmp(this->n_points,1);
        tmp = this->YLK*xc;
        for (int i = 0; i< this->n_points;i++){	this->sf[i] = tmp(i,0);}
        tmp.resize(0,0);
    }
    void surfaceGenX(const MatrixXd &xc, const MatrixXd &yc, const MatrixXd &zc)
    //// Spherical harmonics-based generation of surface
    //// generate X position vector for all points on triangular mesh
    //// using a new set of coefficients (synthesis)
    {
        // calculate X
        //int nc = (this->L_max + 1) *(this->L_max + 1);
        Xu.col(0)= this->YLK*xc;
        Xu.col(1)= this->YLK*yc;
        Xu.col(2)= this->YLK*zc;
        for (int i = 0; i< this->n_points;i++)	// we copy to the double** X to read it from shape class which doesn't know about Eigen
        {
            this->X[i][0] = Xu(i,0);
            this->X[i][1] = Xu(i,1);
            this->X[i][2] = Xu(i,2);
        }
    }
    void surfaceGen(const MatrixXd &xc, const MatrixXd &yc, const MatrixXd &zc)
    //// Spherical harmonics-based generation of surface
    //// generate the X position vector for all points on triangular mesh
    //// using a new set of coefficients (synthesis)
    //// generates first and higher derivatives (first and second funamental forms)
    //// surface metric and normals
    {
        // calculate X
        //int nc = (this->L_max + 1) *(this->L_max + 1);
        Xu.col(0)= this->YLK*xc;
        Xu.col(1)= this->YLK*yc;
        Xu.col(2)= this->YLK*zc;
        for (int i = 0; i< this->n_points;i++)	// we copy to double** X to read it from shape class which doesn't know about Eigen
        {
            this->X[i][0] = Xu(i,0);
            this->X[i][1] = Xu(i,1);
            this->X[i][2] = Xu(i,2);
        }
        //
        this->Xt.col(0) = this->Y_T*xc;
        this->Xt.col(1) = this->Y_T*yc;
        this->Xt.col(2) = this->Y_T*zc;
        //if(KKDEBUG){std::cout<<"Xt"<<std::endl<<"length: "<<Xt.rows()<<std::endl<<Xt<<std::endl;}
        /// calculate Xp
        this->Xp.col(0) = this->Y_P*xc;
        this->Xp.col(1) = this->Y_P*yc;
        this->Xp.col(2) = this->Y_P*zc;
        //if(KKDEBUG){std::cout<<"Xp"<<std::endl<<"length: "<<Xp.rows()<<std::endl<<Xp<<std::endl;}
        /// calculate Xpp
        this->Xpp.col(0) = this->Y_PP*xc;
        this->Xpp.col(1) = this->Y_PP*yc;
        this->Xpp.col(2) = this->Y_PP*zc;
        //if(KKDEBUG){std::cout<<"Xpp"<<std::endl<<"length: "<<Xpp.rows()<<std::endl<<Xpp<<std::endl;}
        /// calculate Xtt
        this->Xtt.col(0) = this->Y_TT*xc;
        this->Xtt.col(1) = this->Y_TT*yc;
        this->Xtt.col(2) = this->Y_TT*zc;
        //if(KKDEBUG){std::cout<<"Xtt"<<std::endl<<"length: "<<Xtt.rows()<<std::endl<<Xtt<<std::endl;}
        /// calculate Xtp
        this->Xtp.col(0) = this->Y_TP*xc;
        this->Xtp.col(1) = this->Y_TP*yc;
        this->Xtp.col(2) = this->Y_TP*zc;
        //if(KKDEBUG){std::cout<<"Xtp"<<std::endl<<"length: "<<Xtp.rows()<<std::endl<<Xtp<<std::endl;}
        // Calculate first fundamental form
        this->E = (Xt.col(0).array())*(Xt.col(0).array()) +(Xt.col(1).array())*(Xt.col(1).array()) +(Xt.col(2).array())*(Xt.col(2).array());
        this->F = (Xt.col(0).array())*(Xp.col(0).array()) +(Xt.col(1).array())*(Xp.col(1).array()) +(Xt.col(2).array())*(Xp.col(2).array());
        this->G = (Xp.col(0).array())*(Xp.col(0).array()) +(Xp.col(1).array())*(Xp.col(1).array()) +(Xp.col(2).array())*(Xp.col(2).array());
        //if(KKDEBUG){std::cout<<"E"<<std::endl<<"length: "<<E.rows()<<std::endl<<E<<std::endl;}
        // Calculate  surface metric and normals
        this->SSn = (this->E.array()*this->G.array()-this->F.array()*this->F.array()).array().sqrt();	// the surface metric
        this->nu.col(0) = (Xt.col(1).array()*Xp.col(2).array()-Xt.col(2).array()*Xp.col(1).array()).array()/SSn.array();
        this->nu.col(1) = (Xt.col(2).array()*Xp.col(0).array()-Xt.col(0).array()*Xp.col(2).array()).array()/SSn.array();
        this->nu.col(2) = (Xt.col(0).array()*Xp.col(1).array()-Xt.col(1).array()*Xp.col(0).array()).array()/SSn.array();
        
        
        // Calculate second fundamental form
        this->L = (Xtt.col(0).array())*(nu.col(0).array()) +(Xtt.col(1).array())*(nu.col(1).array()) +(Xtt.col(2).array())*(nu.col(2).array());
        this->M = (Xtp.col(0).array())*(nu.col(0).array()) +(Xtp.col(1).array())*(nu.col(1).array()) +(Xtp.col(2).array())*(nu.col(2).array());
        this->N = (Xpp.col(0).array())*(nu.col(0).array()) +(Xpp.col(1).array())*(nu.col(1).array()) +(Xpp.col(2).array())*(nu.col(2).array());
        
        // Evaluate geometrical properties
        this->Hu = -(E.array()*N.array() + G.array()*L.array() - 2*F.array()*M.array()).array()/(2 * (E.array()*G.array()-F.array()*F.array()).array());
        this->KGu = (L.array()*N.array() - M.array()*M.array()).array()/(E.array()*G.array()-F.array()*F.array()).array();
        // copy Hu and KGu to double** arrays
        for (int i = 0; i< this->n_points;i++)	// we copy to the double** to read it from shape class which doesn't know about Eigen
        {
            this->Han[i] = Hu(i,0);
            this->KGan[i] = KGu(i,0);
        }
    }
    void surfaceGen(int coordix, int cix, double old_val, double new_val)
    //// generate the X position vector for generation of triangular mesh
    //// generates first and higher derivative fields
    {
        /// calculate X
        this->Xu.col(coordix) = this->Xu.col(coordix)-(this->YLK.col(cix) * old_val);	// subtract the effect of the basis vector
        this->Xu.col(coordix) = this->Xu.col(coordix)+(this->YLK.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
        for (int i = 0; i< this->n_points;i++)	// we copy to the double** X to read it from shape class which doesn't know about Eigen
        {
            this->X[i][0] = Xu(i,0);
            this->X[i][1] = Xu(i,1);
            this->X[i][2] = Xu(i,2);
        }
        //
        /// calculate Xt
        this->Xt.col(coordix) = this->Xt.col(coordix)-(this->Y_T.col(cix) * old_val);	// subtract the effect of the basis vector
        this->Xt.col(coordix) = this->Xt.col(coordix)+(this->Y_T.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
        /// calculate Xp
        this->Xp.col(coordix) = this->Xp.col(coordix)-(this->Y_P.col(cix) * old_val);	// subtract the effect of the basis vector
        this->Xp.col(coordix) = this->Xp.col(coordix)+(this->Y_P.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
        /// calculate Xpp
        this->Xpp.col(coordix) = this->Xpp.col(coordix)-(this->Y_PP.col(cix) * old_val);	// subtract the effect of the basis vector
        this->Xpp.col(coordix) = this->Xpp.col(coordix)+(this->Y_PP.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
        /// calculate Xtt
        this->Xtt.col(coordix) = this->Xtt.col(coordix)-(this->Y_TT.col(cix) * old_val);	// subtract the effect of the basis vector
        this->Xtt.col(coordix) = this->Xtt.col(coordix)+(this->Y_TT.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
        /// calculate Xtp
        this->Xtp.col(coordix) = this->Xtp.col(coordix)-(this->Y_TP.col(cix) * old_val);	// subtract the effect of the basis vector
        this->Xtp.col(coordix) = this->Xtp.col(coordix)+(this->Y_TP.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
        
        //if(KKDEBUG){std::cout<<"Xtp"<<std::endl<<"length: "<<Xtp.rows()<<std::endl<<Xtp<<std::endl;}
        // Calculate first fundamental form
        this->E = (Xt.col(0).array())*(Xt.col(0).array()) +(Xt.col(1).array())*(Xt.col(1).array()) +(Xt.col(2).array())*(Xt.col(2).array());
        this->F = (Xt.col(0).array())*(Xp.col(0).array()) +(Xt.col(1).array())*(Xp.col(1).array()) +(Xt.col(2).array())*(Xp.col(2).array());
        this->G = (Xp.col(0).array())*(Xp.col(0).array()) +(Xp.col(1).array())*(Xp.col(1).array()) +(Xp.col(2).array())*(Xp.col(2).array());
        //if(KKDEBUG){std::cout<<"E"<<std::endl<<"length: "<<E.rows()<<std::endl<<E<<std::endl;}
        // Calculate  the surface metric and normal
        this->SSn = (this->E.array()*this->G.array()-this->F.array()*this->F.array()).array().sqrt();	// the surface metric
        this->nu.col(0) = (Xt.col(1).array()*Xp.col(2).array()-Xt.col(2).array()*Xp.col(1).array()).array()/SSn.array();
        this->nu.col(1) = (Xt.col(2).array()*Xp.col(0).array()-Xt.col(0).array()*Xp.col(2).array()).array()/SSn.array();
        this->nu.col(2) = (Xt.col(0).array()*Xp.col(1).array()-Xt.col(1).array()*Xp.col(0).array()).array()/SSn.array();
        
        
        // Calculate the second fundamental form
        this->L = (Xtt.col(0).array())*(nu.col(0).array()) +(Xtt.col(1).array())*(nu.col(1).array()) +(Xtt.col(2).array())*(nu.col(2).array());
        this->M = (Xtp.col(0).array())*(nu.col(0).array()) +(Xtp.col(1).array())*(nu.col(1).array()) +(Xtp.col(2).array())*(nu.col(2).array());
        this->N = (Xpp.col(0).array())*(nu.col(0).array()) +(Xpp.col(1).array())*(nu.col(1).array()) +(Xpp.col(2).array())*(nu.col(2).array());
        
        // Evaluate geometrical properties
        this->Hu = -(E.array()*N.array() + G.array()*L.array() - 2*F.array()*M.array()).array()/(2 * (E.array()*G.array()-F.array()*F.array()).array());
        this->KGu = (L.array()*N.array() - M.array()*M.array()).array()/(E.array()*G.array()-F.array()*F.array()).array();
        // copy Hu and KGu to double** arrays
        for (int i = 0; i< this->n_points;i++)	// we copy to the double** to read it from shape class which doesn't know about Eigen
        {
            this->Han[i] = Hu(i,0);
            this->KGan[i] = KGu(i,0);
        }
    }
    void update_tri()
    {
        // calculates geometric properties based on Frank Julicher's PhD thesis appendix b
        this->A = 0.0;
        this->V = 0.0;
        // written out in full for readability
        for(int i = 0;i<this->n_faces;i++)
        {		//loop over the faces
            x1[i] = X[f0[i]-1][0];
            y1[i] = X[f0[i]-1][1];
            z1[i] = X[f0[i]-1][2];
            x2[i] = X[f1[i]-1][0];
            y2[i] = X[f1[i]-1][1];
            z2[i] = X[f1[i]-1][2];
            x3[i] = X[f2[i]-1][0];
            y3[i] = X[f2[i]-1][1];
            z3[i] = X[f2[i]-1][2];
            q[i][0] = x2[i]-x1[i];q[i][1]=y2[i]-y1[i];q[i][2]=z2[i]-z1[i];
            r[i][0] = x3[i]-x1[i];r[i][1]=y3[i]-y1[i];r[i][2]=z3[i]-z1[i];
            //take the cross product
            crossqpr[i][0]= (q[i][1]*r[i][2])-(q[i][2]*r[i][1]);
            crossqpr[i][1]= (q[i][2]*r[i][0])-(q[i][0]*r[i][2]);
            crossqpr[i][2]= (q[i][0]*r[i][1])-(q[i][1]*r[i][0]);
            twoA[i]= sqrt((crossqpr[i][0])*(crossqpr[i][0])+(crossqpr[i][1])*(crossqpr[i][1])+(crossqpr[i][2])*(crossqpr[i][2]));
            A+=twoA[i]/2;	// accumulate the areas of surface triangles
            //normalized surface normal vector
            n[i][0] =crossqpr[i][0]/twoA[i]; n[i][1] =crossqpr[i][1]/twoA[i];n[i][2] =crossqpr[i][2]/twoA[i];
            
            V+=-1.0/3.0*(n[i][0]*x1[i] + n[i][1]*y1[i]+ n[i][2]*z1[i])*twoA[i]/2.0;
        }
        double Vo = 4.0/3.0*PI*sqrt((A/4.0/PI)*(A/4.0/PI)*(A/4.0/PI));
        this->A = std::abs(this->A);
        this->V = std::abs(this->V);
        this->v = std::abs(this->V/Vo);
        
        Eb = 0.0;		//initialize bending energy
        h =  0.0;		//initialize total curvature
        if(_curv_calc)
        {
            // Let's calculate the curvatures at all vertices
            // this calculation heavily relies on precalculated quantities
            // that are generated with the matlabl file output triangulated props.m
            // this information is loaded into memory when reading mesh_info.txt.
            // variable meanings are similar to the matlab file.
            double theta,Lij = 0.0;
            int nm;		//number of member triangles
            double* n1;
            double* n2;    // working normals
            int *rr;
            int *indx;
            int *r1;
            int *r2;
            int *r3;
            double Va[3];
            double Ve[3];
            double Vfar[3];
            double P1;
            //double dA[_gdim];	//surface area associated with each vertex
            Eb = 0.0;		//initialize bending energy
            h =  0.0;		//initialize total curvature
            
            for(int ix=0;ix<this->n_points;ix++)// loop over the vertices
            {
                dA[ix] = 0.0;
                H[ix] = 0.0;
                vN[ix][0] = 0.0;
                vN[ix][1] = 0.0;
                vN[ix][2] = 0.0;
                Va[0] = X[ix][0];Va[1] = X[ix][1];Va[2] = X[ix][2];
                nm = _vip[ix].n;
                rr   = _vip[ix].r;
                indx = _vip[ix].indx;
                r1   = _vip[ix].r1;
                r2   = _vip[ix].r2;
                r3   = _vip[ix].r3;
                
                for(int rix=0;rix<nm;rix++)	//loop over edges
                {
                    // Calculate the surface normal at this vertex
                    this->vN[ix][0] += n[rr[rix]][0]/nm;
                    this->vN[ix][1] += n[rr[rix]][1]/nm;
                    this->vN[ix][2] += n[rr[rix]][2]/nm;
                    
                    // Calculate the curvature
                    Ve[0] = X[indx[rix]][0];
                    Ve[1] = X[indx[rix]][1];
                    Ve[2] = X[indx[rix]][2];
                    Vfar[0] = X[r3[rix]][0];Vfar[1] = X[r3[rix]][1];Vfar[2] = X[r3[rix]][2];
                    Lij = sqrt(  (Va[0]-Ve[0])*(Va[0]-Ve[0]) + (Va[1]-Ve[1])*(Va[1]-Ve[1]) + (Va[2]-Ve[2])*(Va[2]-Ve[2]) ) ;
                    n1 = n[r1[rix]];
                    n2 = n[r2[rix]];
                    theta = acos((n1[0]*n2[0] + n1[1]*n2[1]+ n1[2]*n2[2]));
                    P1 = -(n1[0]*Va[0]+n1[1]*Va[1]+n1[2]*Va[2]);
                    if ( ((n1[0]*Vfar[0]+n1[1]*Vfar[1]+n1[2]*Vfar[2])+P1)>=0 ){	H[ix]  += -Lij*theta/4;	}
                    else
                    {
                        H[ix]  += Lij*theta/4;	// in this case H is still actually H*dA
                    }
                    dA[ix] += twoA[rr[rix]]/2/3;
                }
                
                h+=H[ix]/A;
                Eb+=(2*H[ix]*H[ix]/dA[ix])/8/PI;
                H[ix] = H[ix]/dA[ix];	// now we have the real H value
                
            }
        }
    }
    
    void tri_props(vector< vector<double> > XX, double &a, double &v)
    // should be a static function. just calculates area and volume for some XX
    {
        // calculates geometric properties
        a = 0.0;
        v = 0.0;
        // written out in full for readability
        for(int i = 0;i<this->n_faces;i++)
        {		//loop over the faces
            x1[i] = XX[f0[i]-1][0];
            y1[i] = XX[f0[i]-1][1];
            z1[i] = XX[f0[i]-1][2];
            x2[i] = XX[f1[i]-1][0];
            y2[i] = XX[f1[i]-1][1];
            z2[i] = XX[f1[i]-1][2];
            x3[i] = XX[f2[i]-1][0];
            y3[i] = XX[f2[i]-1][1];
            z3[i] = XX[f2[i]-1][2];
            q[i][0] = x2[i]-x1[i];q[i][1]=y2[i]-y1[i];q[i][2]=z2[i]-z1[i];
            r[i][0] = x3[i]-x1[i];r[i][1]=y3[i]-y1[i];r[i][2]=z3[i]-z1[i];
            //take the cross product
            crossqpr[i][0]= (q[i][1]*r[i][2])-(q[i][2]*r[i][1]);
            crossqpr[i][1]= (q[i][2]*r[i][0])-(q[i][0]*r[i][2]);
            crossqpr[i][2]= (q[i][0]*r[i][1])-(q[i][1]*r[i][0]);
            twoA[i]= sqrt((crossqpr[i][0])*(crossqpr[i][0])+(crossqpr[i][1])*(crossqpr[i][1])+(crossqpr[i][2])*(crossqpr[i][2]));
            a+=twoA[i]/2;	// accumulate the areas of surface triangles
            
            //normalized surface normal vector
            n[i][0] =crossqpr[i][0]/twoA[i]; n[i][1] =crossqpr[i][1]/twoA[i];n[i][2] =crossqpr[i][2]/twoA[i];
            
            v+=-1.0/3.0*(n[i][0]*x1[i] + n[i][1]*y1[i]+ n[i][2]*z1[i])*twoA[i]/2.0;
        }
        a = std::abs(a);
        v = std::abs(v);
    }
    
    int self_intersect()
    {
        double* V0;
        double* V1;
        double* V2;
        double* U0;
        double* U1;
        double* U2;
        for(int i = 0;i<this->n_faces-1;i++)
        {
            V0 = (X[int(this->f0[i]-1)]);
            V1 = (X[int(this->f1[i]-1)]);
            V2 = (X[int(this->f2[i]-1)]);
            for(int j = i+1;j<this->n_faces;j++)
            {
                U0 = (X[int(this->f0[j]-1)]);
                U1 = (X[int(this->f1[j]-1)]);
                U2 = (X[int(this->f2[j]-1)]);
                if(tri_tri_intersect(V0,V1,V2,U0,U1,U2))
                {
                    int sv = share_vert((this->f0[i]-1),(this->f1[i]-1),(this->f2[i]-1),(this->f0[j]-1),(this->f1[j]-1),(this->f2[j]-1));
                    if(sv==0)
                    {
                        return 1;
                    }
                }
            }
        }
        return 0;
    }
    // constructor
    spherical_mesh(int L = 6, int ix = 2)
    {
        L_max = L;
        _curv_calc = true;
        tri_n = ix;
        if (tri_n == 1)
        {
            fn_str = "ico2.tri";
            n_points = 162;
            n_faces  = 320;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_01,p);
            t.resize(n_points,1);d2eig(THETA_01,t);
            f0 = F0_01;
            f1 = F1_01;
            f2 = F2_01;
        }
        if (tri_n == 2)
        {
            fn_str = "ico3.tri";
            n_points = 642;
            n_faces  = 1280;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_02,p);
            t.resize(n_points,1);d2eig(THETA_02,t);
            f0 = F0_02;
            f1 = F1_02;
            f2 = F2_02;
        }
        
        if (tri_n == 3)
        {
            fn_str = "uni900.tri";
            n_points = 900;
            n_faces  = 1796;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_03,p);
            t.resize(n_points,1);d2eig(THETA_03,t);
            f0 = F0_03;
            f1 = F1_03;
            f2 = F2_03;
        }
        if (tri_n == 4)
        {
            fn_str = "ico4.tri";
            n_points = 2562;
            n_faces  = 5120;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_04,p);
            t.resize(n_points,1);d2eig(THETA_04,t);
            f0 = F0_04;
            f1 = F1_04;
            f2 = F2_04;
        }
        if (tri_n == 5)
        {
            fn_str = "uni10k.tri";
            n_points = 10000;
            n_faces  = 19996;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_05,p);
            t.resize(n_points,1);d2eig(THETA_05,t);
            f0 = F0_05;
            f1 = F1_05;
            f2 = F2_05;
        }
        if (tri_n == 6)
        {
            fn_str = "ico5.tri";
            n_points = 10242;
            n_faces  = 20480;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_06,p);
            t.resize(n_points,1);d2eig(THETA_06,t);
            f0 = F0_06;
            f1 = F1_06;
            f2 = F2_06;
        }
        MatrixXd PLK, P_T;
        ylk_cos_sin_bosh(L_max, p, t,YLK, PLK);
        ylk_cos_sin_dphi_bosh(L_max, p, t, PLK, Y_P);
        ylk_cos_sin_dphiphi_bosh(L_max, p, t, PLK, Y_PP);
        ylk_cos_sin_dtheta_bosh(L_max, p, t, PLK, Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(L_max, p, t, PLK, P_T, Y_TP);
        ylk_cos_sin_dthetatheta_bosh(L_max, p, t, P_T, Y_TT);
        
        this->mesh_init();		// initialize arrays used to calculate surface properties
        
        // cleanup
        PLK.resize(0,0);
        P_T.resize(0,0);
    }
    
    spherical_mesh(const spherical_mesh& that)	// copy constructor
    {
        L_max = that.L_max;
        _curv_calc = true;
        tri_n = that.tri_n;
        if (tri_n == 1)
        {
            fn_str = "ico2.tri";
            n_points = 162;
            n_faces  = 320;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_01,p);
            t.resize(n_points,1);d2eig(THETA_01,t);
            f0 = F0_01;
            f1 = F1_01;
            f2 = F2_01;
        }
        if (tri_n == 2)
        {
            fn_str = "ico3.tri";
            n_points = 642;
            n_faces  = 1280;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_02,p);
            t.resize(n_points,1);d2eig(THETA_02,t);
            f0 = F0_02;
            f1 = F1_02;
            f2 = F2_02;
        }
        
        if (tri_n == 3)
        {
            fn_str = "uni900.tri";
            n_points = 900;
            n_faces  = 1796;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_03,p);
            t.resize(n_points,1);d2eig(THETA_03,t);
            f0 = F0_03;
            f1 = F1_03;
            f2 = F2_03;
        }
        if (tri_n == 4)
        {
            fn_str = "ico4.tri";
            n_points = 2562;
            n_faces  = 5120;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_04,p);
            t.resize(n_points,1);d2eig(THETA_04,t);
            f0 = F0_04;
            f1 = F1_04;
            f2 = F2_04;
        }
        if (tri_n == 5)
        {
            fn_str = "uni10k.tri";
            n_points = 10000;
            n_faces  = 19996;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_05,p);
            t.resize(n_points,1);d2eig(THETA_05,t);
            f0 = F0_05;
            f1 = F1_05;
            f2 = F2_05;
        }
        if (tri_n == 6)
        {
            fn_str = "ico5.tri";
            n_points = 10242;
            n_faces  = 20480;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_06,p);
            t.resize(n_points,1);d2eig(THETA_06,t);
            f0 = F0_06;
            f1 = F1_06;
            f2 = F2_06;
        }
        MatrixXd PLK, P_T;
        ylk_cos_sin_bosh(L_max, p, t,YLK, PLK);
        ylk_cos_sin_dphi_bosh(L_max, p, t, PLK, Y_P);
        ylk_cos_sin_dphiphi_bosh(L_max, p, t, PLK, Y_PP);
        ylk_cos_sin_dtheta_bosh(L_max, p, t, PLK, Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(L_max, p, t, PLK, P_T, Y_TP);
        ylk_cos_sin_dthetatheta_bosh(L_max, p, t, P_T, Y_TT);
        this->mesh_init();		// initialize arrays used to calculate surface properties
        // cleanup
        PLK.resize(0,0);
        P_T.resize(0,0);
        
    }
    spherical_mesh& operator=(const spherical_mesh& that)	//copy assignment operator -- copies member variables and then initializes
    {
        L_max = that.L_max;
        _curv_calc = true;
        tri_n = that.tri_n;
        if (tri_n == 1)
        {
            fn_str = "ico2.tri";
            n_points = 162;
            n_faces  = 320;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_01,p);
            t.resize(n_points,1);d2eig(THETA_01,t);
            f0 = F0_01;
            f1 = F1_01;
            f2 = F2_01;
        }
        if (tri_n == 2)
        {
            fn_str = "ico3.tri";
            n_points = 642;
            n_faces  = 1280;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_02,p);
            t.resize(n_points,1);d2eig(THETA_02,t);
            f0 = F0_02;
            f1 = F1_02;
            f2 = F2_02;
        }
        
        if (tri_n == 3)
        {
            fn_str = "uni900.tri";
            n_points = 900;
            n_faces  = 1796;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_03,p);
            t.resize(n_points,1);d2eig(THETA_03,t);
            f0 = F0_03;
            f1 = F1_03;
            f2 = F2_03;
        }
        if (tri_n == 4)
        {
            fn_str = "ico4.tri";
            n_points = 2562;
            n_faces  = 5120;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_04,p);
            t.resize(n_points,1);d2eig(THETA_04,t);
            f0 = F0_04;
            f1 = F1_04;
            f2 = F2_04;
        }
        if (tri_n == 5)
        {
            fn_str = "uni10k.tri";
            n_points = 10000;
            n_faces  = 19996;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_05,p);
            t.resize(n_points,1);d2eig(THETA_05,t);
            f0 = F0_05;
            f1 = F1_05;
            f2 = F2_05;
        }
        if (tri_n == 6)
        {
            fn_str = "ico5.tri";
            n_points = 10242;
            n_faces  = 20480;
            // calculate the basis matrix and import the faces for the visualization
            p.resize(n_points,1);d2eig(PHI_06,p);
            t.resize(n_points,1);d2eig(THETA_06,t);
            f0 = F0_06;
            f1 = F1_06;
            f2 = F2_06;
        }
        MatrixXd PLK, P_T;
        ylk_cos_sin_bosh(L_max, p, t,YLK, PLK);
        ylk_cos_sin_dphi_bosh(L_max, p, t, PLK, Y_P);
        ylk_cos_sin_dphiphi_bosh(L_max, p, t, PLK, Y_PP);
        ylk_cos_sin_dtheta_bosh(L_max, p, t, PLK, Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(L_max, p, t, PLK, P_T, Y_TP);
        ylk_cos_sin_dthetatheta_bosh(L_max, p, t, P_T, Y_TT);
        
        this->mesh_init();		// initialize arrays used to calculate surface properties
        // cleanup
        PLK.resize(0,0);
        P_T.resize(0,0);
        
        return *this;
        
    }
    ~spherical_mesh()		// destructor
    {
        this->p.resize(0,0);
        this->t.resize(0,0);
        this->YLK.resize(0,0);
        this->Y_P.resize(0,0);
        this->Y_T.resize(0,0);
        this->Y_TT.resize(0,0);
        this->Y_PP.resize(0,0);
        this->Y_TP.resize(0,0);
        this->Xu.resize(0,0);
        this->Xt.resize(0,0);
        this->Xp.resize(0,0);
        this->Xtt.resize(0,0);
        this->Xtp.resize(0,0);
        this->Xpp.resize(0,0);
        this->E.resize(0,0);
        this->F.resize(0,0);
        this->G.resize(0,0);
        this->L.resize(0,0);
        this->M.resize(0,0);
        this->N.resize(0,0);
        this->nu.resize(0,0);
        this->SSn.resize(0,0);
        this->dXn.resize(0,0);
        this->Hu.resize(0,0);
        this->KGu.resize(0,0);
        
        free(this->sf);
        free(this->X);
        free(this->H);
        free(this->Han);
        free(this->KGan);
        free(this->dA);
        free(this->_vip);
        free(this->q);
        free(this->r);
        free(this->crossqpr);
        free(this->n);
        free(this->twoA);
        free(this->x1);
        free(this->y1);
        free(this->z1);
        free(this->x2);
        free(this->y2);
        free(this->z2);
        free(this->x3);
        free(this->y3);
        free(this->z3);
    }
    
public:
    ostream& disp(ostream& os, const spherical_mesh& b)
    {
        /**/
        os<<"L_max: "<<b.L_max<<std::endl;
        os<<"Dimensions: "<<std::endl;
        os<<"YLK: "<<b.YLK.rows()<<" x "<<b.YLK.cols()<<std::endl;
        os<<"n_points: "<<b.n_points<<std::endl;
        os<<"n_faces: "<<b.n_faces<<std::endl;
        return os;
    };
    
    friend ostream& operator<<(ostream& os, const spherical_mesh& s)
    {
        
        os<<s.tri_n<<std::endl;
        os<<s.L_max<<std::endl;
        //
        return os;
    };
    
    
    friend istream& operator>>(istream& in, spherical_mesh& s)
    {
        in>>s.tri_n;
        in>>s.L_max;
        // initialize accordingly
        s._curv_calc = true;
        if (s.tri_n == 1)
        {
            s.fn_str = "ico2.tri";
            s.n_points = 162;
            s.n_faces  = 320;
            // calculate the basis matrix and import the faces for the visualization
            s.p.resize(s.n_points,1);d2eig(PHI_01,s.p);
            s.t.resize(s.n_points,1);d2eig(THETA_01,s.t);
            s.f0 = F0_01;
            s.f1 = F1_01;
            s.f2 = F2_01;
        }
        if (s.tri_n == 2)
        {
            s.fn_str = "ico3.tri";
            s.n_points = 642;
            s.n_faces  = 1280;
            // calculate the basis matrix and import the faces for the visualization
            s.p.resize(s.n_points,1);d2eig(PHI_02,s.p);
            s.t.resize(s.n_points,1);d2eig(THETA_02,s.t);
            s.f0 = F0_02;
            s.f1 = F1_02;
            s.f2 = F2_02;
        }
        
        if (s.tri_n == 3)
        {
            s.fn_str = "uni900.tri";
            s.n_points = 900;
            s.n_faces  = 1796;
            // calculate the basis matrix and import the faces for the visualization
            s.p.resize(s.n_points,1);d2eig(PHI_03,s.p);
            s.t.resize(s.n_points,1);d2eig(THETA_03,s.t);
            s.f0 = F0_03;
            s.f1 = F1_03;
            s.f2 = F2_03;
        }
        if (s.tri_n == 4)
        {
            s.fn_str = "ico4.tri";
            s.n_points = 2562;
            s.n_faces  = 5120;
            // calculate the basis matrix and import the faces for the visualization
            s.p.resize(s.n_points,1);d2eig(PHI_04,s.p);
            s.t.resize(s.n_points,1);d2eig(THETA_04,s.t);
            s.f0 = F0_04;
            s.f1 = F1_04;
            s.f2 = F2_04;
        }
        if (s.tri_n == 5)
        {
            s.fn_str = "uni10k.tri";
            s.n_points = 10000;
            s.n_faces  = 19996;
            // calculate the basis matrix and import the faces for the visualization
            s.p.resize(s.n_points,1);d2eig(PHI_05,s.p);
            s.t.resize(s.n_points,1);d2eig(THETA_05,s.t);
            s.f0 = F0_05;
            s.f1 = F1_05;
            s.f2 = F2_05;
        }
        if (s.tri_n == 6)
        {
            s.fn_str = "ico5.tri";
            s.n_points = 10242;
            s.n_faces  = 20480;
            // calculate the basis matrix and import the faces for the visualization
            s.p.resize(s.n_points,1);d2eig(PHI_06,s.p);
            s.t.resize(s.n_points,1);d2eig(THETA_06,s.t);
            s.f0 = F0_06;
            s.f1 = F1_06;
            s.f2 = F2_06;
        }
        MatrixXd PLK, P_T;
        ylk_cos_sin_bosh(s.L_max, s.p, s.t,s.YLK, PLK);
        ylk_cos_sin_dphi_bosh(s.L_max, s.p, s.t, PLK, s.Y_P);
        ylk_cos_sin_dphiphi_bosh(s.L_max, s.p, s.t, PLK, s.Y_PP);
        ylk_cos_sin_dtheta_bosh(s.L_max, s.p, s.t, PLK, s.Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(s.L_max, s.p, s.t, PLK, P_T, s.Y_TP);
        ylk_cos_sin_dthetatheta_bosh(s.L_max, s.p, s.t, P_T, s.Y_TT);
        
        s.mesh_init();		// initialize arrays used to calculate surface properties
        s.update_tri();
        // cleanup
        PLK.resize(0,0);
        P_T.resize(0,0);
        
        //
        return in;
    };
};



class sh_basis	// main basis class initializing to Gaussian quadrature weights and visualization basis
{
public:
    int dim;
    int L_max;
    //member variables for Gaussian quadrature
    MatrixXd p,t,wt,wp;	// column matrices holding the GA
    MatrixXd phi, theta,w;
    MatrixXd YLK, Y_P, Y_T, Y_TT, Y_PP, Y_TP;
    
    // member variables for visualization only
    //MatrixXd YLK_tri;
    
    // constructor
    sh_basis(int arg1 = 6, int arg2 = 30)
    {
        //replicate the basis generation code in matlab
        dim = arg2;
        L_max = arg1;
        t.resize(dim,1);
        wt.resize(dim,1);
        p.resize(dim,1);
        wp.resize(dim,1);
        Eigen::MatrixXd wphi,wtheta, PLK,P_T;
        
        // prepare the Gaussian base points and weights
        gaussquad(0, PI, t, wt);	// generate the base points and weights for Gaussian quadrature
        gaussquad(0, 2*PI, p, wp);	// generate the base points and weights for Gaussian quadrature
        meshgrid(p, t, phi, theta);
        meshgrid(wp, wt, wphi, wtheta);
        reshape(wphi,wphi.rows()*wphi.cols(), 1);
        reshape(wtheta,wtheta.rows()*wtheta.cols(), 1);
        w = wphi.array() * wtheta.array();
        // calculate the basis vectors and the derivatives and store them
        ylk_cos_sin_bosh(L_max, phi, theta,YLK, PLK);
        ylk_cos_sin_dphi_bosh(L_max, phi, theta, PLK, Y_P);
        ylk_cos_sin_dphiphi_bosh(L_max, phi, theta, PLK, Y_PP);
        ylk_cos_sin_dtheta_bosh(L_max, phi, theta, PLK, Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(L_max, phi, theta, PLK, P_T, Y_TP);
        ylk_cos_sin_dthetatheta_bosh(L_max, phi, theta, P_T, Y_TT);
        //cleanup
        wphi.resize(0,0);
        wtheta.resize(0,0);
        PLK.resize(0,0);
        P_T.resize(0,0);
        
    }
    sh_basis(const sh_basis& that)		// copy constructor
    {
        //replicate the basis generation code in matlab
        dim = that.dim;
        L_max = that.L_max;
        t.resize(dim,1);
        wt.resize(dim,1);
        p.resize(dim,1);
        wp.resize(dim,1);
        Eigen::MatrixXd wphi,wtheta, PLK,P_T;
        
        // prepare the Gaussian base points and weights
        gaussquad(0, PI, t, wt);	// generate the base points and weights for Gaussian quadrature
        gaussquad(0, 2*PI, p, wp);	// generate the base points and weights for Gaussian quadrature
        meshgrid(p, t, phi, theta);
        meshgrid(wp, wt, wphi, wtheta);
        reshape(wphi,wphi.rows()*wphi.cols(), 1);
        reshape(wtheta,wtheta.rows()*wtheta.cols(), 1);
        w = wphi.array() * wtheta.array();
        // calculate the basis vectors and the derivatives and store them
        ylk_cos_sin_bosh(L_max, phi, theta,YLK, PLK);
        ylk_cos_sin_dphi_bosh(L_max, phi, theta, PLK, Y_P);
        ylk_cos_sin_dphiphi_bosh(L_max, phi, theta, PLK, Y_PP);
        ylk_cos_sin_dtheta_bosh(L_max, phi, theta, PLK, Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(L_max, phi, theta, PLK, P_T, Y_TP);
        ylk_cos_sin_dthetatheta_bosh(L_max, phi, theta, P_T, Y_TT);
        
        //cleanup
        wphi.resize(0,0);
        wtheta.resize(0,0);
        PLK.resize(0,0);
        P_T.resize(0,0);
        
    }
    sh_basis& operator=(const sh_basis& that)	// copy assignment
    {
        this->dim = that.dim;
        this->L_max = that.L_max;
        t.resize(dim,1);
        wt.resize(dim,1);
        p.resize(dim,1);
        wp.resize(dim,1);
        Eigen::MatrixXd wphi,wtheta, PLK,P_T;
        
        // prepare the Gaussian base points and weights
        gaussquad(0, PI, t, wt);	// generate the base points and weights for Gaussian quadrature
        gaussquad(0, 2*PI, p, wp);	// generate the base points and weights for Gaussian quadrature
        meshgrid(p, t, phi, theta);
        meshgrid(wp, wt, wphi, wtheta);
        reshape(wphi,wphi.rows()*wphi.cols(), 1);
        reshape(wtheta,wtheta.rows()*wtheta.cols(), 1);
        w = wphi.array() * wtheta.array();
        // calculate the basis vectors and the derivatives, and store them
        ylk_cos_sin_bosh(L_max, phi, theta,YLK, PLK);
        ylk_cos_sin_dphi_bosh(L_max, phi, theta, PLK, Y_P);
        ylk_cos_sin_dphiphi_bosh(L_max, phi, theta, PLK, Y_PP);
        ylk_cos_sin_dtheta_bosh(L_max, phi, theta, PLK, Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(L_max, phi, theta, PLK, P_T, Y_TP);
        ylk_cos_sin_dthetatheta_bosh(L_max, phi, theta, P_T, Y_TT);
        
        //cleanup
        wphi.resize(0,0);
        wtheta.resize(0,0);
        PLK.resize(0,0);
        P_T.resize(0,0);
        return *this;
    }
    ~sh_basis()
    {
        p.resize(0,0);
        t.resize(0,0);
        wt.resize(0,0);
        wp.resize(0,0);
        phi.resize(0,0);
        theta.resize(0,0);
        w.resize(0,0);
        YLK.resize(0,0);
        Y_P.resize(0,0);
        Y_T.resize(0,0);
        Y_TT.resize(0,0);
        Y_PP.resize(0,0);
        Y_TP.resize(0,0);
    };
    
public:
    ostream& disp(ostream& os, const sh_basis& b)
    {
        /**/
        os<<"L_max: "<<b.L_max<<std::endl;
        os<<"dim  : "<<b.dim<<std::endl;
        os<<"Dimensions: "<<std::endl;
        os<<"			t: "<<b.t.rows()<<" x "<<b.t.cols()<<std::endl;
        os<<"		theta: "<<b.theta.rows()<<" x "<<b.theta.cols()<<std::endl;
        os<<"			p: "<<b.p.rows()<<" x "<<b.p.cols()<<std::endl;
        os<<"		  phi: "<<b.phi.rows()<<" x "<<b.phi.cols()<<std::endl;
        os<<"		   wt: "<<b.wt.rows()<<" x "<<b.wt.cols()<<std::endl;
        os<<"		   wp: "<<b.wp.rows()<<" x "<<b.wp.cols()<<std::endl;
        os<<"		    w: "<<b.w.rows()<<" x "<<b.w.cols()<<std::endl;
        os<<"		   YLK: "<<b.YLK.rows()<<" x "<<b.YLK.cols()<<std::endl;
        os<<"		   Y_P: "<<b.Y_P.rows()<<" x "<<b.Y_P.cols()<<std::endl;
        os<<"		   Y_T: "<<b.Y_T.rows()<<" x "<<b.Y_T.cols()<<std::endl;
        //os<<"	   YLK_tri: "<<b.YLK_tri.rows()<<" x "<<b.YLK_tri.cols()<<std::endl;
        
        return os;
    };
    friend ostream& operator<<(ostream& os, const sh_basis& s)
    {
        
        os<<s.dim<<std::endl;
        os<<s.L_max<<std::endl;
        //
        return os;
    };
    
    
    friend istream& operator>>(istream& in, sh_basis& s)
    {
        in>>s.dim;
        in>>s.L_max;
        // initialize accordingly
        s.t.resize(s.dim,1);
        s.wt.resize(s.dim,1);
        s.p.resize(s.dim,1);
        s.wp.resize(s.dim,1);
        Eigen::MatrixXd wphi,wtheta, PLK,P_T;
        
        // prepare the Gaussian base points and weights
        gaussquad(0, PI, s.t, s.wt);	// generate the base points and weights for Gaussian quadrature
        gaussquad(0, 2*PI, s.p, s.wp);	// generate the base points and weights for Gaussian quadrature
        meshgrid(s.p, s.t, s.phi, s.theta);
        meshgrid(s.wp, s.wt, wphi, wtheta);
        reshape(wphi,wphi.rows()*wphi.cols(), 1);
        reshape(wtheta,wtheta.rows()*wtheta.cols(), 1);
        s.w = wphi.array() * wtheta.array();
        // calculate the basis vectors and the derivatives and store them
        ylk_cos_sin_bosh(s.L_max, s.phi, s.theta,s.YLK, PLK);
        ylk_cos_sin_dphi_bosh(s.L_max, s.phi, s.theta, PLK, s.Y_P);
        ylk_cos_sin_dphiphi_bosh(s.L_max, s.phi, s.theta, PLK, s.Y_PP);
        ylk_cos_sin_dtheta_bosh(s.L_max, s.phi, s.theta, PLK, s.Y_T, P_T);
        ylk_cos_sin_dthetaphi_bosh(s.L_max, s.phi, s.theta, PLK, P_T, s.Y_TP);
        ylk_cos_sin_dthetatheta_bosh(s.L_max, s.phi, s.theta, P_T, s.Y_TT);
        
        //cleanup
        wphi.resize(0,0);
        wtheta.resize(0,0);
        PLK.resize(0,0);
        P_T.resize(0,0);
        
        
        //
        return in;
    };
    
};

////////////////////////////////////////////////////////////////////////////////
class shp_surface
// class depends on spherical_mesh and sh_basis. uses Bosh
// 2000 normalized spherical harmonics
// Defines the actual surface with all scalar fields. It maintains an updated
// triangular mesh (through spherical_mesh) only for self-intersection tests
// and visualization (which usually --eg. vtk-- requires a surface tesselation). Also, since
// surface properties can be calculated approximately on the triangular mesh
// it comes in handy as a self-check and comparison with the analytical formula(for academic purposes).
// shp_surface is able to: enforce symmetry, change L_max, rotate surface around itself,
// calculate total mean and Gaussian curvature, surface area, volume, bending
// energy and surface normals.
//
{
public:
	//member variables
	string name;
	bool needsUpdating;
	double A, V, v, Eb,wb, h, T;
	MatrixXd xc, yc, zc, sfmx;	// these spherical harmonics coefficients represent the Cartesian coordinates as functions on the sphere
								// sf represents all other scalar fields
	vector<string> sf_tags;
	// for efficiency we will initialize (see mx_init()) the intermediate matrices below (all used in the function update()), and keep them in memory
	MatrixXd Xu, Xt, Xp, Xtt, Xpp, Xtp, E, F, G, L, M, N, n,SSn, dXn, H, KG;
	//
	sh_basis* b;	//
	spherical_mesh * sm;
	//debug matrices
	
	//methods
	void enforce_mirror_yz(std::vector<bool> &ixc, std::vector<bool> &iyc, std::vector<bool> &izc)
	{	// use the symmetry rules to enforce the mirror plane
		// Assuming xc1m1 is leading in x, and yc00 is leading in y, and zc11 in z
		// values set to zero are labeled false in their corresponding vectors ixc, iyc and izc
		int nc = this->xc.rows();
		ixc.resize(nc);std::fill(ixc.begin(), ixc.end(), true);
		iyc.resize(nc);std::fill(iyc.begin(), iyc.end(), true);
		izc.resize(nc);std::fill(izc.begin(), izc.end(), true);
		int counter = 0;
		for (int l = 0; l<this->b->L_max;l++)
		{
			for(int m = -l;m<=l;m++)
			{
				if (m >=0) {this->xc(counter) = 0;ixc[counter] = false;}
				if (m < 0) 
				{
					this->yc(counter) = 0;iyc[counter] = false;
					this->zc(counter) = 0;izc[counter] = false;
				}
				counter++;
			}
		}
		this->update();
	}
	void enforce_mirror_yz()
	{	// use the symmetry rules to enforce the mirror plane
		// Assuming xc1m1 is leading in x, and yc00 is leading in y, and zc11 in z
		// values set to zero are labeled false in their corresponding vectors ixc, iyc and izc
		int counter = 0;
		for (int l = 0; l<this->b->L_max;l++)
		{
			for(int m = -l;m<=l;m++)
			{
				if (m >=0) {this->xc(counter) = 0;}
				if (m < 0) 
				{
					this->yc(counter) = 0;
					this->zc(counter) = 0;
				}
				counter++;
			}
		}
		this->update();
	}
	void rotate_around_self(double*ang)
	{
		MatrixXd R, Rg, Rb, Ra, C, cr;
		rot_mx(ang[0], 2, Ra);	//rotate around y
		rot_mx(ang[1], 1, Rb);  //rotate around z
		rot_mx(ang[2], 2, Rg);  //rotate around y again
		double tmpx, tmpy, tmpz;
		tmpx = this->xc(0);
		tmpy = this->yc(0);
		tmpz = this->zc(0);
		this->xc(0) = 0.0;
		this->yc(0) = 0.0;
		this->zc(0) = 0.0;
		C.resize(this->xc.rows(), 3);
		C.col(0) = xc;
		C.col(1) = yc;
		C.col(2) = zc;
		C.transpose();
		R = (Ra*Rb*Rg);
		cr = C*R;
		cr.transpose();
		this->xc = cr.col(0);
		this->yc = cr.col(1);
		this->zc = cr.col(2);

		this->xc(0) = tmpx;
		this->yc(0) = tmpy;
		this->zc(0) = tmpz;
	}
	void initializeToSphere()
	{
		int nc = (this->b->L_max + 1)*(this->b->L_max + 1);
		double R = 1;
		this->xc.resize(nc,1);this->xc.fill(0.0);
		this->yc.resize(nc,1);this->yc.fill(0.0);
		this->zc.resize(nc,1);this->zc.fill(0.0);
		this->sfmx.resize(0,0);
		// initialize to a sphere
		this->xc(1) = R/N_LK_bosh(1,1);
		this->yc(2) = -R/N_LK_bosh(1,-1);
		this->zc(3) = -R/N_LK_bosh(1,0);
		this->sf_tags.clear();
		this->sf_tags.push_back("x");
		this->sf_tags.push_back("y");
		this->sf_tags.push_back("z");
		this->needsUpdating = YES;
	}
	void mx_init()
	{

		Xu = MatrixXd(this->b->dim * this->b->dim, 3);	// 3-vector of columns
		Xt = MatrixXd(this->b->dim * this->b->dim, 3);	// 3-vector of columns
		Xp = MatrixXd(this->b->dim * this->b->dim, 3);	// 3-vector of columns
		Xtt = MatrixXd(this->b->dim * this->b->dim, 3);	// 3-vector of columns
		Xpp = MatrixXd(this->b->dim * this->b->dim, 3);	// 3-vector of columns
		Xtp = MatrixXd(this->b->dim * this->b->dim, 3);	// 3-vector of columns

		E = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		F = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		G = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		L = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		M = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		N = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 

		SSn = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		n = MatrixXd(this->b->dim * this->b->dim, 3);	// 3-vector of columns resulting from the cross product of Xt and Xp

		dXn = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		H = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
		KG = MatrixXd(this->b->dim * this->b->dim, 1);	// 1-vector 
	}
	void set_clks(MatrixXd xc, MatrixXd yc, MatrixXd zc)
	{
		if(xc.rows()!= this->xc.rows())
		{
			int new_L = int(sqrt(double(xc.rows())))-1;
			this->set_L_max(new_L);
		}
		this->xc.resize(0,0);
		this->yc.resize(0,0);
		this->zc.resize(0,0);
		
		this->xc = xc;
		this->yc = yc;
		this->zc = zc;
		this->needsUpdating = 1;
		this->update();
	}
	void update()
	{
		/// calculate X
		this->Xu.col(0) = this->b->YLK*this->xc;
		this->Xu.col(1) = this->b->YLK*this->yc;
		this->Xu.col(2) = this->b->YLK*this->zc;
//if(KKDEBUG){std::cout<<"X"<<std::endl<<"length: "<<X.rows()<<std::endl<<X<<std::endl;}
		/// calculate Xt
		this->Xt.col(0) = this->b->Y_T*this->xc;
		this->Xt.col(1) = this->b->Y_T*this->yc;
		this->Xt.col(2) = this->b->Y_T*this->zc;
//if(KKDEBUG){std::cout<<"Xt"<<std::endl<<"length: "<<Xt.rows()<<std::endl<<Xt<<std::endl;}
		/// calculate Xp
		this->Xp.col(0) = this->b->Y_P*this->xc;
		this->Xp.col(1) = this->b->Y_P*this->yc;
		this->Xp.col(2) = this->b->Y_P*this->zc;
//if(KKDEBUG){std::cout<<"Xp"<<std::endl<<"length: "<<Xp.rows()<<std::endl<<Xp<<std::endl;}
		/// calculate Xpp
		this->Xpp.col(0) = this->b->Y_PP*this->xc;
		this->Xpp.col(1) = this->b->Y_PP*this->yc;
		this->Xpp.col(2) = this->b->Y_PP*this->zc;
//if(KKDEBUG){std::cout<<"Xpp"<<std::endl<<"length: "<<Xpp.rows()<<std::endl<<Xpp<<std::endl;}
		/// calculate Xtt
		this->Xtt.col(0) = this->b->Y_TT*this->xc;
		this->Xtt.col(1) = this->b->Y_TT*this->yc;
		this->Xtt.col(2) = this->b->Y_TT*this->zc;
//if(KKDEBUG){std::cout<<"Xtt"<<std::endl<<"length: "<<Xtt.rows()<<std::endl<<Xtt<<std::endl;}
		/// calculate Xtp
		this->Xtp.col(0) = this->b->Y_TP*this->xc;
		this->Xtp.col(1) = this->b->Y_TP*this->yc;
		this->Xtp.col(2) = this->b->Y_TP*this->zc;
//if(KKDEBUG){std::cout<<"Xtp"<<std::endl<<"length: "<<Xtp.rows()<<std::endl<<Xtp<<std::endl;}
		// Calculate first fundamental form
		 this->E = (Xt.col(0).array())*(Xt.col(0).array()) +(Xt.col(1).array())*(Xt.col(1).array()) +(Xt.col(2).array())*(Xt.col(2).array()); 
		 this->F = (Xt.col(0).array())*(Xp.col(0).array()) +(Xt.col(1).array())*(Xp.col(1).array()) +(Xt.col(2).array())*(Xp.col(2).array()); 
		 this->G = (Xp.col(0).array())*(Xp.col(0).array()) +(Xp.col(1).array())*(Xp.col(1).array()) +(Xp.col(2).array())*(Xp.col(2).array()); 
//if(KKDEBUG){std::cout<<"E"<<std::endl<<"length: "<<E.rows()<<std::endl<<E<<std::endl;}
		// Calculate  the surface metric and normal
		this->SSn = (this->E.array()*this->G.array()-this->F.array()*this->F.array()).array().sqrt();	// the surface metric
		this->n.col(0) = (Xt.col(1).array()*Xp.col(2).array()-Xt.col(2).array()*Xp.col(1).array()).array()/SSn.array();
		this->n.col(1) = (Xt.col(2).array()*Xp.col(0).array()-Xt.col(0).array()*Xp.col(2).array()).array()/SSn.array();
		this->n.col(2) = (Xt.col(0).array()*Xp.col(1).array()-Xt.col(1).array()*Xp.col(0).array()).array()/SSn.array();

		
		// Calculate the second fundamental form
		 this->L = (Xtt.col(0).array())*(n.col(0).array()) +(Xtt.col(1).array())*(n.col(1).array()) +(Xtt.col(2).array())*(n.col(2).array());
		 this->M = (Xtp.col(0).array())*(n.col(0).array()) +(Xtp.col(1).array())*(n.col(1).array()) +(Xtp.col(2).array())*(n.col(2).array());
		 this->N = (Xpp.col(0).array())*(n.col(0).array()) +(Xpp.col(1).array())*(n.col(1).array()) +(Xpp.col(2).array())*(n.col(2).array());

		// Evaluate geometrical properties
		this->dXn = (Xu.col(0).array())*(n.col(0).array()) +(Xu.col(1).array())*(n.col(1).array()) +(Xu.col(2).array())*(n.col(2).array());
		
		this->V = abs(double(((this->b->w.array()*this->dXn.array()*this->SSn.array()).sum()))/3);
		this->A = double((this->b->w.array()*this->SSn.array()).sum());

		this->H = -(E.array()*N.array() + G.array()*L.array() - 2*F.array()*M.array()).array()/(2 * (E.array()*G.array()-F.array()*F.array()).array());
		this->KG = (L.array()*N.array() - M.array()*M.array()).array()/(E.array()*G.array()-F.array()*F.array()).array();
		this->h = ((H.array()*this->b->w.array()*SSn.array()).sum())/this->A;       // total mean curvature
		this->T = (KG.array()*this->b->w.array()*SSn.array()).sum()/4/PI;           // total Gaussian curvature (constant for topology)
		this->Eb = (4*H.array()*H.array()*this->b->w.array()*SSn.array()).sum()/16/PI;
		//this->wb = this->Eb/8/PI;
		this->v = this->V/(4.0/3.0*PI*sqrt((this->A/4.0/PI)*(this->A/4.0/PI)*(this->A/4.0/PI)));
		this->needsUpdating = false;

	}
	void update_fast()
	{
		/// calculate X
		this->Xu.col(0) = this->b->YLK*this->xc;
		this->Xu.col(1) = this->b->YLK*this->yc;
		this->Xu.col(2) = this->b->YLK*this->zc;
//if(KKDEBUG){std::cout<<"X"<<std::endl<<"length: "<<X.rows()<<std::endl<<X<<std::endl;}
		/// calculate Xt
		this->Xt.col(0) = this->b->Y_T*this->xc;
		this->Xt.col(1) = this->b->Y_T*this->yc;
		this->Xt.col(2) = this->b->Y_T*this->zc;
//if(KKDEBUG){std::cout<<"Xt"<<std::endl<<"length: "<<Xt.rows()<<std::endl<<Xt<<std::endl;}
		/// calculate Xp
		this->Xp.col(0) = this->b->Y_P*this->xc;
		this->Xp.col(1) = this->b->Y_P*this->yc;
		this->Xp.col(2) = this->b->Y_P*this->zc;
//if(KKDEBUG){std::cout<<"Xp"<<std::endl<<"length: "<<Xp.rows()<<std::endl<<Xp<<std::endl;}
		/// calculate Xpp
		this->Xpp.col(0) = this->b->Y_PP*this->xc;
		this->Xpp.col(1) = this->b->Y_PP*this->yc;
		this->Xpp.col(2) = this->b->Y_PP*this->zc;
//if(KKDEBUG){std::cout<<"Xpp"<<std::endl<<"length: "<<Xpp.rows()<<std::endl<<Xpp<<std::endl;}
		/// calculate Xtt
		this->Xtt.col(0) = this->b->Y_TT*this->xc;
		this->Xtt.col(1) = this->b->Y_TT*this->yc;
		this->Xtt.col(2) = this->b->Y_TT*this->zc;
//if(KKDEBUG){std::cout<<"Xtt"<<std::endl<<"length: "<<Xtt.rows()<<std::endl<<Xtt<<std::endl;}
		/// calculate Xtp
		this->Xtp.col(0) = this->b->Y_TP*this->xc;
		this->Xtp.col(1) = this->b->Y_TP*this->yc;
		this->Xtp.col(2) = this->b->Y_TP*this->zc;
//if(KKDEBUG){std::cout<<"Xtp"<<std::endl<<"length: "<<Xtp.rows()<<std::endl<<Xtp<<std::endl;}
		// Calculate first fundamental form
		 this->E = (Xt.col(0).array())*(Xt.col(0).array()) +(Xt.col(1).array())*(Xt.col(1).array()) +(Xt.col(2).array())*(Xt.col(2).array()); 
		 this->F = (Xt.col(0).array())*(Xp.col(0).array()) +(Xt.col(1).array())*(Xp.col(1).array()) +(Xt.col(2).array())*(Xp.col(2).array()); 
		 this->G = (Xp.col(0).array())*(Xp.col(0).array()) +(Xp.col(1).array())*(Xp.col(1).array()) +(Xp.col(2).array())*(Xp.col(2).array()); 
//if(KKDEBUG){std::cout<<"E"<<std::endl<<"length: "<<E.rows()<<std::endl<<E<<std::endl;}
		// Calculate  the surface metric and normal
		this->SSn = (this->E.array()*this->G.array()-this->F.array()*this->F.array()).array().sqrt();	// the surface metric
		this->n.col(0) = (Xt.col(1).array()*Xp.col(2).array()-Xt.col(2).array()*Xp.col(1).array()).array()/SSn.array();
		this->n.col(1) = (Xt.col(2).array()*Xp.col(0).array()-Xt.col(0).array()*Xp.col(2).array()).array()/SSn.array();
		this->n.col(2) = (Xt.col(0).array()*Xp.col(1).array()-Xt.col(1).array()*Xp.col(0).array()).array()/SSn.array();

		
		// Calculate the second fundamental form
		 this->L = (Xtt.col(0).array())*(n.col(0).array()) +(Xtt.col(1).array())*(n.col(1).array()) +(Xtt.col(2).array())*(n.col(2).array());
		 this->M = (Xtp.col(0).array())*(n.col(0).array()) +(Xtp.col(1).array())*(n.col(1).array()) +(Xtp.col(2).array())*(n.col(2).array());
		 this->N = (Xpp.col(0).array())*(n.col(0).array()) +(Xpp.col(1).array())*(n.col(1).array()) +(Xpp.col(2).array())*(n.col(2).array());

		// Evaluate local mean curvature
		this->H = -(E.array()*N.array() + G.array()*L.array() - 2*F.array()*M.array()).array()/(2 * (E.array()*G.array()-F.array()*F.array()).array());
		

	}
	void update_fast_volume()
	{
		/// calculate X
		this->Xu.col(0) = this->b->YLK*this->xc;
		this->Xu.col(1) = this->b->YLK*this->yc;
		this->Xu.col(2) = this->b->YLK*this->zc;
//if(KKDEBUG){std::cout<<"X"<<std::endl<<"length: "<<X.rows()<<std::endl<<X<<std::endl;}
		/// calculate Xt
		this->Xt.col(0) = this->b->Y_T*this->xc;
		this->Xt.col(1) = this->b->Y_T*this->yc;
		this->Xt.col(2) = this->b->Y_T*this->zc;
//if(KKDEBUG){std::cout<<"Xt"<<std::endl<<"length: "<<Xt.rows()<<std::endl<<Xt<<std::endl;}
		/// calculate Xp
		this->Xp.col(0) = this->b->Y_P*this->xc;
		this->Xp.col(1) = this->b->Y_P*this->yc;
		this->Xp.col(2) = this->b->Y_P*this->zc;
		// Calculate first fundamental form
		 this->E = (Xt.col(0).array())*(Xt.col(0).array()) +(Xt.col(1).array())*(Xt.col(1).array()) +(Xt.col(2).array())*(Xt.col(2).array()); 
		 this->F = (Xt.col(0).array())*(Xp.col(0).array()) +(Xt.col(1).array())*(Xp.col(1).array()) +(Xt.col(2).array())*(Xp.col(2).array()); 
		 this->G = (Xp.col(0).array())*(Xp.col(0).array()) +(Xp.col(1).array())*(Xp.col(1).array()) +(Xp.col(2).array())*(Xp.col(2).array()); 
		// Calculate  the surface metric and normal
		this->SSn = (this->E.array()*this->G.array()-this->F.array()*this->F.array()).array().sqrt();	// the surface metric
		this->n.col(0) = (Xt.col(1).array()*Xp.col(2).array()-Xt.col(2).array()*Xp.col(1).array()).array()/SSn.array();
		this->n.col(1) = (Xt.col(2).array()*Xp.col(0).array()-Xt.col(0).array()*Xp.col(2).array()).array()/SSn.array();
		this->n.col(2) = (Xt.col(0).array()*Xp.col(1).array()-Xt.col(1).array()*Xp.col(0).array()).array()/SSn.array();
		// Evaluate geometrical properties
		this->dXn = (Xu.col(0).array())*(n.col(0).array()) +(Xu.col(1).array())*(n.col(1).array()) +(Xu.col(2).array())*(n.col(2).array());
		this->V = abs(double(((this->b->w.array()*this->dXn.array()*this->SSn.array()).sum()))/3);

	}
	void update(int coordix, int cix, double old_val, double new_val)
	{
		/// calculate X
		this->Xu.col(coordix) = this->Xu.col(coordix)-(this->b->YLK.col(cix) * old_val);	// subtract the effect of the basis vector
		this->Xu.col(coordix) = this->Xu.col(coordix)+(this->b->YLK.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
		/// calculate Xt
		this->Xt.col(coordix) = this->Xt.col(coordix)-(this->b->Y_T.col(cix) * old_val);	// subtract the effect of the basis vector
		this->Xt.col(coordix) = this->Xt.col(coordix)+(this->b->Y_T.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
		/// calculate Xp
		this->Xp.col(coordix) = this->Xp.col(coordix)-(this->b->Y_P.col(cix) * old_val);	// subtract the effect of the basis vector
		this->Xp.col(coordix) = this->Xp.col(coordix)+(this->b->Y_P.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
		/// calculate Xpp
		this->Xpp.col(coordix) = this->Xpp.col(coordix)-(this->b->Y_PP.col(cix) * old_val);	// subtract the effect of the basis vector
		this->Xpp.col(coordix) = this->Xpp.col(coordix)+(this->b->Y_PP.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
		/// calculate Xtt
		this->Xtt.col(coordix) = this->Xtt.col(coordix)-(this->b->Y_TT.col(cix) * old_val);	// subtract the effect of the basis vector
		this->Xtt.col(coordix) = this->Xtt.col(coordix)+(this->b->Y_TT.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
		/// calculate Xtp
		this->Xtp.col(coordix) = this->Xtp.col(coordix)-(this->b->Y_TP.col(cix) * old_val);	// subtract the effect of the basis vector
		this->Xtp.col(coordix) = this->Xtp.col(coordix)+(this->b->Y_TP.col(cix) * new_val);	// add the effect with the new value of the clk coefficient
		// Calculate first fundamental form
		 this->E = (Xt.col(0).array())*(Xt.col(0).array()) +(Xt.col(1).array())*(Xt.col(1).array()) +(Xt.col(2).array())*(Xt.col(2).array()); 
		 this->F = (Xt.col(0).array())*(Xp.col(0).array()) +(Xt.col(1).array())*(Xp.col(1).array()) +(Xt.col(2).array())*(Xp.col(2).array()); 
		 this->G = (Xp.col(0).array())*(Xp.col(0).array()) +(Xp.col(1).array())*(Xp.col(1).array()) +(Xp.col(2).array())*(Xp.col(2).array()); 
		// Calculate  the surface metric and normal
		this->SSn = (this->E.array()*this->G.array()-this->F.array()*this->F.array()).array().sqrt();	// the surface metric
		this->n.col(0) = (Xt.col(1).array()*Xp.col(2).array()-Xt.col(2).array()*Xp.col(1).array()).array()/SSn.array();
		this->n.col(1) = (Xt.col(2).array()*Xp.col(0).array()-Xt.col(0).array()*Xp.col(2).array()).array()/SSn.array();
		this->n.col(2) = (Xt.col(0).array()*Xp.col(1).array()-Xt.col(1).array()*Xp.col(0).array()).array()/SSn.array();
		// Calculate the second fundamental form
		 this->L = (Xtt.col(0).array())*(n.col(0).array()) +(Xtt.col(1).array())*(n.col(1).array()) +(Xtt.col(2).array())*(n.col(2).array());
		 this->M = (Xtp.col(0).array())*(n.col(0).array()) +(Xtp.col(1).array())*(n.col(1).array()) +(Xtp.col(2).array())*(n.col(2).array());
		 this->N = (Xpp.col(0).array())*(n.col(0).array()) +(Xpp.col(1).array())*(n.col(1).array()) +(Xpp.col(2).array())*(n.col(2).array());
		// Evaluate geometrical properties
		this->dXn = (Xu.col(0).array())*(n.col(0).array()) +(Xu.col(1).array())*(n.col(1).array()) +(Xu.col(2).array())*(n.col(2).array());
		this->V = abs(double(((this->b->w.array()*this->dXn.array()*this->SSn.array()).sum()))/3);
		this->A = double((this->b->w.array()*this->SSn.array()).sum());

		this->H = -(E.array()*N.array() + G.array()*L.array() - 2*F.array()*M.array()).array()/(2 * (E.array()*G.array()-F.array()*F.array()).array());
		this->KG = (L.array()*N.array() - M.array()*M.array()).array()/(E.array()*G.array()-F.array()*F.array()).array();
		this->h = ((H.array()*this->b->w.array()*SSn.array()).sum())/this->A;       // total mean curvature
		this->T = (KG.array()*this->b->w.array()*SSn.array()).sum()/4/PI;           // total Gaussian curvature (constant for topology)
		this->Eb = (4*H.array()*H.array()*this->b->w.array()*SSn.array()).sum()/16/PI;
		//this->wb = this->Eb/8/PI;
		this->v = this->V/(4.0/3.0*PI*sqrt((this->A/4.0/PI)*(this->A/4.0/PI)*(this->A/4.0/PI)));
		this->needsUpdating = false;

	}
	void update_tri(){this->sm->update_tri();}
	double area()  {if(this->needsUpdating){this->update();}return this->A;}
	double volume(){if(this->needsUpdating){this->update();}return this->V;}
	double energy(){if(this->needsUpdating){this->update();}return this->Eb;}
	int self_intersect(){ return sm->self_intersect();}
	void surfaceGen(){this->sm->surfaceGen(xc, yc, zc);}
	void surfaceGen(int coordix, int cix, double old_val, double new_val){this->sm->surfaceGen( coordix,  cix,  old_val,  new_val);}
	void sfGen(int n){this->sm->sfGen(sfmx.col(n));}
	void set_L_max(const int L)
	{
		int old_L_max = this->b->L_max;
		int dim = this->b->dim;
		int old_nc = this->xc.rows();
		if(old_L_max!=L)	// 
		{
			int new_nc = (L+1)*(L+1);  // new number of coefficients
			// generate a new basis and delete the old one
			sh_basis* tmp_b = this->b;
			this->b = new sh_basis(L,dim);	// make a new basis
			delete tmp_b;		// free memory of the old basis;
			// resize the scalar field matrix
			MatrixXd tmp_sfmx = this->sfmx;
			this->sfmx.resize(new_nc,sfmx.cols());
			
			MatrixXd tmp;

			// 
			if(L>old_L_max)
			{
				tmp = this->xc;xc.conservativeResize(new_nc,1);xc.fill(0.0);for(int i=0;i<old_nc;i++){xc(i) = tmp(i);}
				tmp = this->yc;yc.conservativeResize(new_nc,1);yc.fill(0.0);for(int i=0;i<old_nc;i++){yc(i) = tmp(i);}
				tmp = this->zc;zc.conservativeResize(new_nc,1);zc.fill(0.0);for(int i=0;i<old_nc;i++){zc(i) = tmp(i);}
				this->sfmx.fill(0.0);
				for(int r = 0;r<old_nc;r++)
					for(int c = 0;c<tmp_sfmx.cols();c++)
					{this->sfmx(r,c) = tmp_sfmx(r,c);}
			}
			if(L<old_L_max)
			{
				tmp = this->xc;xc.conservativeResize(new_nc,1);xc.fill(0.0);for(int i=0;i<new_nc;i++){xc(i) = tmp(i);}
				tmp = this->yc;yc.conservativeResize(new_nc,1);yc.fill(0.0);for(int i=0;i<new_nc;i++){yc(i) = tmp(i);}
				tmp = this->zc;zc.conservativeResize(new_nc,1);zc.fill(0.0);for(int i=0;i<new_nc;i++){zc(i) = tmp(i);}
				for(int r = 0;r<new_nc;r++)
					for(int c = 0;c<tmp_sfmx.cols();c++)
					{this->sfmx(r,c) = tmp_sfmx(r,c);}
			}
			this->mx_init();
			this->needsUpdating = YES;
			this->update();
			
			// generate a new spherical mesh object and delete the old one
			spherical_mesh*tmp_sm = this->sm;
			sm = new spherical_mesh(L,this->sm->tri_n);
			delete tmp_sm;
		}
	}
	void set_new_basis(const int L, const int dim)
	{
		int old_L_max = this->b->L_max;
		int old_nc = this->xc.rows();
		int new_nc = (L+1)*(L+1);
		sh_basis* tmp_b = this->b;
		this->b = new sh_basis(L,dim);	// make a new basis
		delete tmp_b;		// free memory of the old basis;
		if(L>old_L_max)
		{
			MatrixXd tmp = this->xc;xc.conservativeResize(new_nc,1);xc.fill(0.0);for(int i=0;i<old_nc;i++){xc(i) = tmp(i);}
					 tmp = this->yc;yc.conservativeResize(new_nc,1);yc.fill(0.0);for(int i=0;i<old_nc;i++){yc(i) = tmp(i);}
					 tmp = this->zc;zc.conservativeResize(new_nc,1);zc.fill(0.0);for(int i=0;i<old_nc;i++){zc(i) = tmp(i);}
		}
		if(L<old_L_max)
		{
			 MatrixXd tmp = this->xc;xc.conservativeResize(new_nc,1);xc.fill(0.0);for(int i=0;i<new_nc;i++){xc(i) = tmp(i);}
					  tmp = this->yc;yc.conservativeResize(new_nc,1);yc.fill(0.0);for(int i=0;i<new_nc;i++){yc(i) = tmp(i);}
					  tmp = this->zc;zc.conservativeResize(new_nc,1);zc.fill(0.0);for(int i=0;i<new_nc;i++){zc(i) = tmp(i);}
		}
		this->mx_init();
		this->needsUpdating = YES;
		this->update();
		if(old_L_max!=L)	// then we need to also update this->sm
		{
			spherical_mesh*tmp = this->sm;
			sm = new spherical_mesh(L,this->sm->tri_n);
			delete tmp;
		}

	}
	void set_new_spherical_mesh(int n)
	{
		spherical_mesh *tmp = this->sm;
		sm = new spherical_mesh(this->b->L_max,n);
		delete tmp;
	}
	void set_new_spherical_mesh(int L, int n)
	{
		spherical_mesh*tmp = this->sm;
		sm = new spherical_mesh(L,n);
		delete tmp;
	}

	void set_new_gdim(int gdim)
	{// recalculate basis to only change the Gaussian base point mesh
		if(gdim!=this->b->dim)
		{
			sh_basis* tmp_basis = this->b;
			this->b = new sh_basis(tmp_basis->L_max,gdim);
			delete tmp_basis;
			this->mx_init();
			this->update();
		}
	}
	void scale(double s)
	{
		double xtmp = this->xc(0);
		double ytmp = this->yc(0);
		double ztmp = this->zc(0);
		this->xc = this->xc * s;
		this->yc = this->yc * s;
		this->zc = this->zc * s;
		this->xc(0) = xtmp;
		this->yc(0) = ytmp;
		this->zc(0) = ztmp;
	}
	void center_to_zero()
	{
		this->xc(0) = 0.0;
		this->yc(0) = 0.0;
		this->zc(0) = 0.0;
	}
	void flush_fields_after_L(int Lin)
	{	// sets all field coefficients with L larger than Lin to zero
		if (Lin<this->b->L_max)
		{
			int n = (Lin + 1) * (Lin + 1);
			for(int i = n+1;i<this->sfmx.rows();i++)
			{
				for(int c = 0;c<this->sfmx.cols();c++)
				{
					this->sfmx(i,c) = 0.0;
				}
			}
		}

	}
	void flush_after_L(int Lin)
	{
		if (Lin<this->b->L_max)
		{
			int n = (Lin + 1) * (Lin + 1);
			for(int i = n+1;i<this->xc.rows();i++)
			{
				this->xc(i,0) = 0.0;
				this->yc(i,0) = 0.0;
				this->zc(i,0) = 0.0;
			}
		}

	}
	// file read write
	/*void read(const char* filename)	// hasn't been tested yet
	{	// reads the coefficients from a file (same format as clks.txt)
		// and generates the surface -- needs to handle the case when L is not equal to current L_max
		int iReturn;
		int n_shapes, n_components;
		int L = 0;
		int nc;
		float xclk, yclk, zclk, mxel;
		char cstr[50000];
		string str;
		FILE *f;
		errno_t err = fopen_s(&f,filename,"r");
		fscanf_s(f,"n_shapes = %d\n",&n_shapes);
		fscanf_s(f,"L_max = %d\n",&L);
		if(KKDEBUG){std::cout<<"Found shape with L_max = "<<L<<std::endl;}
		fscanf_s(f,"n_components = %d\n",&n_components);
		nc = (L+1)*(L+1);
		
		//// the next line contains tab separated strings for the tags associated with the coefficient lists
		this->sf_tags.clear();
		fgets(cstr , 50000 , f); 		// read contents of the tags line into str
		str = string(cstr);
		std::string delim = "\t";
		Tokenize(str, sf_tags, delim);
		// handle case when L is not equal to L_max
		if(L>(this->b->L_max))
		{
			if(KKDEBUG){std::cout<<"Please be patient: New L_max ("<<L<<") exceeds current L_max. Calculating new basis...";}
			set_new_basis(L, this->b->dim);
			if(KKDEBUG){std::cout<<"Done"<<std::endl;}
		}
		if(L<(this->b->L_max))
		{
			if(KKDEBUG){std::cout<<"Please be patient: New L_max ("<<L<<") below current L_max. Resizing basis...";}
			set_new_basis(L, this->b->dim);
			if(KKDEBUG){std::cout<<"Done"<<std::endl;}
		}
		if(L<6)
		{
			if(KKDEBUG){std::cout<<"Please be patient: New L_max ("<<L<<") below current L_max. \nResizing basis to minimum (L = 6) ...";}
			set_new_basis(6, this->b->dim);
			if(KKDEBUG){std::cout<<"Done"<<std::endl;}
		}

		if(n_components>3)
		{
			sfmx.resize(nc,n_components-3);sfmx.fill(0.0);
		}else
		{
			sfmx.resize(0,0);
		}
		this->xc.fill(0.0);
		this->yc.fill(0.0);
		this->zc.fill(0.0);
		for(int i=0;i<nc;i++)
		{
			iReturn = fscanf_s(f,"%e\t%e\t%e",&xclk,&yclk,&zclk);
			xc(i,0) = xclk;yc(i,0) = yclk;zc(i,0) = zclk;
			for(int j=0;j<(n_components-3);j++)
			{
				iReturn = fscanf_s(f,"\t%e",&mxel);
				sfmx(i,j) = mxel;
			}
		}

		fclose(f);
		this->needsUpdating = YES;
	}
	*/
	void read(const char* filename)
	{	// reads the coefficients from a file (same format as clks.txt)
		// and generates the surface -- needs to handle the case when L is not equal to current L_max
		int iReturn;
		int n_shapes, n_components;
		int L = 0;
		float xclk, yclk, zclk, mxel;
		char cstr[50000];
		string str;
		FILE *f = fopen(filename,"r");
		fscanf(f,"n_shapes = %d\n",&n_shapes);
		fscanf(f,"L_max = %d\n",&L);
		if(KKDEBUG){std::cout<<"Found shape with L_max = "<<L<<std::endl;}
		fscanf(f,"n_components = %d\n",&n_components);
		//// the next line contains tab separated strings for the tags associated with the coefficient lists
		this->sf_tags.clear();
		fgets(cstr , 50000 , f); 		// read contents of the tags line into str
		str = string(cstr);
		std::string delim = "\t";
		Tokenize(str, sf_tags, delim);
		//int nco = (this->b->L_max +1) * (this->b->L_max +1);		// this is the nc that we  need
		int nc = (L+1) * (L+1);
		this->xc.resize(nc,1);
		this->yc.resize(nc,1);
		this->zc.resize(nc,1);
		this->xc.fill(0.0);
		this->yc.fill(0.0);
		this->zc.fill(0.0);
		if(n_components>3){sfmx.resize(nc,n_components-3);sfmx.fill(0.0);}
		else{sfmx.resize(0,0);}

		// handle case when L is not equal to L_max
		if(L>(this->b->L_max))
		{
			if(KKDEBUG){std::cout<<"L_max in file larger than current L_max -- expanding basis set)"<<std::endl;}
		}
		if(L<(this->b->L_max))
		{
			if(KKDEBUG){std::cout<<"L_max in file smaller than current L_max -- decreasing basis size"<<std::endl;}	//nc = (L+1) * (L+1);
		}
		for(int i=0;i<nc;i++)
		{
			iReturn = fscanf(f,"%e\t%e\t%e",&xclk,&yclk,&zclk);
			xc(i,0) = xclk;yc(i,0) = yclk;zc(i,0) = zclk;
			for(int j=0;j<(n_components-3);j++)
			{
				iReturn = fscanf(f,"\t%e",&mxel);
				sfmx(i,j) = mxel;
			}
		}
		fclose(f);
		this->needsUpdating = YES;
		this->set_new_basis(L, this->b->dim);

	}
	void read_trunc(const char* filename)
	{	// reads the coefficients from a file (same format as clks.txt)
		// and generates the surface while keeping the same L_max (truncating or padding if needed)
		int iReturn;
		int n_shapes, n_components;
		int L = 0;
		float xclk, yclk, zclk, mxel;
		char cstr[50000];
		string str;
		FILE *f = fopen(filename,"r");
		fscanf(f,"n_shapes = %d\n",&n_shapes);
		fscanf(f,"L_max = %d\n",&L);
		if(KKDEBUG){std::cout<<"Found shape with L_max = "<<L<<std::endl;}
		fscanf(f,"n_components = %d\n",&n_components);
		//// the next line contains tab separated strings for the tags associated with the coefficient lists
		this->sf_tags.clear();
		fgets(cstr , 50000 , f); 		// read contents of the tags line into str
		str = string(cstr);
		std::string delim = "\t";
		Tokenize(str, sf_tags, delim);
		int nco = (this->b->L_max +1) * (this->b->L_max +1);		// this is the nc that we  need
		int nc = (L+1) * (L+1);
		this->xc.resize(nco,1);
		this->yc.resize(nco,1);
		this->zc.resize(nco,1);
		this->xc.fill(0.0);
		this->yc.fill(0.0);
		this->zc.fill(0.0);
		if(n_components>3){sfmx.resize(nco,n_components-3);sfmx.fill(0.0);}
		else{sfmx.resize(0,0);}

		// handle case when L is not equal to L_max
		if(L>(this->b->L_max))
		{
			if(KKDEBUG){std::cout<<"L_max in file larger than current L_max -- truncating (possible loss of shape data)"<<std::endl;}
			for(int i=0;i<nco;i++)
			{
				iReturn = fscanf(f,"%e\t%e\t%e",&xclk,&yclk,&zclk);
				xc(i,0) = xclk;yc(i,0) = yclk;zc(i,0) = zclk;
				for(int j=0;j<(n_components-3);j++)
				{
					iReturn = fscanf(f,"\t%e",&mxel);
					sfmx(i,j) = mxel;
				}
			}
		}
		if(L<=(this->b->L_max))
		{
			if(KKDEBUG){std::cout<<"L_max in file smaller than current L_max -- padding with zero amplitudes"<<std::endl;}	//nc = (L+1) * (L+1);
			for(int i=0;i<nc;i++)
			{
				iReturn = fscanf(f,"%e\t%e\t%e",&xclk,&yclk,&zclk);
				xc(i,0) = xclk;yc(i,0) = yclk;zc(i,0) = zclk;
				for(int j=0;j<(n_components-3);j++)
				{
					iReturn = fscanf(f,"\t%e",&mxel);
					sfmx(i,j) = mxel;
				}
			}
		}

		fclose(f);
		this->needsUpdating = YES;
		//this->set_new_basis(L, this->b->dim);

	}
	void write(const char* filename)
	{
		// writes the coefficients to a file (same format as clks.txt)
		int iReturn, nc;
		double xclk, yclk, zclk;
		int n_components = this->sf_tags.size();
		FILE *f= fopen(filename,"w");
		nc = (this->b->L_max + 1) *(this->b->L_max + 1); 
		fprintf(f,"n_shapes = %d\n",1);
		fprintf(f,"L_max = %d\n",this->b->L_max);
		fprintf(f,"n_components = %d\n",n_components);
		fprintf(f,"x\ty\tz");
		for(int i = 3;i<(this->sf_tags.size());i++){fprintf(f,"\t%s",this->sf_tags[i].c_str());}
		fprintf(f,"\n");
		for(int i=0;i<nc;i++)
		{
			xclk = xc(i,0);yclk = yc(i,0);zclk = zc(i,0);
			iReturn = fprintf(f,"%e\t%e\t%e",xclk,yclk,zclk);
			for(int j=0;j<this->sf_tags.size()-3;j++)
			{
				iReturn = fprintf(f,"\t%e", sfmx(i,j));
			}
			iReturn = fprintf(f,"\n");
		}
		fclose(f);
	}
	//Constructors
	shp_surface(const int L_max, const int dim, const int n = 3)
	{
		this->b = new sh_basis(L_max,dim);
		this->sm = new spherical_mesh(L_max, n); // initialize the triangular mesh
		this->name = "untitled";
		this->needsUpdating = true;
		this->mx_init();
		this->initializeToSphere();
		this->surfaceGen();
	}
	shp_surface(sh_basis *basis)
	{
		this->b = basis;
		this->sm = new spherical_mesh(b->L_max, 2); // initialize the triangular mesh using the default (3) mesh
		this->name = "untitled";
		this->needsUpdating = true;
		this->mx_init();
		this->initializeToSphere();
		this->surfaceGen();
	}
		//shp_surface& shp_surface::operator=( const shp_surface& other ) 
	//{
	//    this->b = new sh_basis(other.b->L_max,other.b->dim);
	//	this->sm = new spherical_mesh(other.b->L_max, other.sm->tri_n); // initialize the triangular mesh
	//	this->name = "untitled";
	//	this->needsUpdating = true;
	//	this->mx_init();
	//	this->initializeToSphere();
	//	this->surfaceGen();
	//	return*this;
   //}
	/*shp_surface()
	{
		b = sh_basis();
		this->name = "untitled";
		this->needsUpdating = TRUE;
		this->initializeToSphere();
		// for visualization and triangular mesh we need the coordinates in a double array
		this->X = (double**)calloc(N_POINTS,sizeof(double*));	for(int i=0;i<N_POINTS;i++){X[i] = (double*)calloc(3,sizeof(double));}
		this->surfaceGen();
	}*/
	// destructor
	~shp_surface()
	{
		xc.resize(0,0);
		yc.resize(0,0);
		zc.resize(0,0);
		sfmx.resize(0,0);
		Xu.resize(0,0);
		Xt.resize(0,0);
		Xp.resize(0,0);
		Xtt.resize(0,0);
		Xpp.resize(0,0);
		Xtp.resize(0,0);
		E.resize(0,0);
		F.resize(0,0);
		G.resize(0,0);
		L.resize(0,0);
		M.resize(0,0);
		N.resize(0,0);
		n.resize(0,0);
		SSn.resize(0,0);
		dXn.resize(0,0);
		H.resize(0,0);
		KG.resize(0,0);
		delete sm;
		//delete b;  // b should not be deleted here (in future: decrement count -- use smart pointer) 
	}

// display
public:
	ostream& disp(ostream& os, const shp_surface& s)
	{
		MatrixXd clks;
		clks.resize(s.xc.rows(),3);
		clks.col(0) = s.xc;
		clks.col(1) = s.yc;
		clks.col(2) = s.zc;
		os<<" --------------------------------------------"<<std::endl;
		os<<" ------------- SHP SURFACE OBJECT -----------"<<std::endl;
		os<<" --------------------------------------------"<<std::endl;
		//os<<s.b<<std::endl;//also print its basis
		os<<"Name: "<<s.name<<std::endl;
		os<<"      xc     yc     zc"<<std::endl<<clks<<std::endl;
		os<<"                    area: "<<s.A<<std::endl;
		os<<"                  volume: "<<s.V<<std::endl;
		os<<"    Total mean curvature: "<<s.h<<std::endl;
		os<<"Total Gaussian curvature: "<<s.T<<std::endl;
		os<<"                      Eb: "<<s.Eb<<std::endl;
		return os;
	};

	friend ostream& operator<<(ostream& os, const shp_surface& s)
	{
		
		os<<*(s.b)<<std::endl;
		os<<*(s.sm)<<std::endl;
		//now export matrices/vectors
		os<<s.xc.rows()<<std::endl;	os<<s.xc.cols()<<std::endl;
		for(int r = 0; r<s.xc.rows();r++){for(int c = 0;c<s.xc.cols();c++){os<<s.xc(r,c)<<std::endl;}};
		
		os<<s.yc.rows()<<std::endl;	os<<s.yc.cols()<<std::endl;
		for(int r = 0; r<s.yc.rows();r++){for(int c = 0;c<s.yc.cols();c++){os<<s.yc(r,c)<<std::endl;}};
		
		os<<s.zc.rows()<<std::endl;	os<<s.zc.cols()<<std::endl;
		for(int r = 0; r<s.zc.rows();r++){for(int c = 0;c<s.zc.cols();c++){os<<s.zc(r,c)<<std::endl;}};
		
		os<<s.sfmx.rows()<<std::endl; os<<s.sfmx.cols()<<std::endl;
		for(int r = 0; r<s.sfmx.rows();r++){for(int c = 0;c<s.sfmx.cols();c++){os<<s.sfmx(r,c)<<std::endl;}};

		os<<s.sf_tags.size()<<std::endl;
		for(int i = 0;i<s.sf_tags.size();i++){os<<s.sf_tags[i]<<std::endl;}
		//
		return os;
	};


	friend istream& operator>>(istream& in, shp_surface& s)
	{
		in>>*(s.b);
		in>>*(s.sm);
		
		//now import matrices/vectors
		int rows, cols, size;
		in>>rows;in>>cols;s.xc.resize(rows,cols);for(int r = 0; r<rows;r++){for(int c = 0;c<cols;c++){in>>s.xc(r,c);}};
		in>>rows;in>>cols;s.yc.resize(rows,cols);for(int r = 0; r<rows;r++){for(int c = 0;c<cols;c++){in>>s.yc(r,c);}};
		in>>rows;in>>cols;s.zc.resize(rows,cols);for(int r = 0; r<rows;r++){for(int c = 0;c<cols;c++){in>>s.zc(r,c);}};
		in>>rows;in>>cols;s.sfmx.resize(rows,cols);for(int r = 0; r<rows;r++){for(int c = 0;c<cols;c++){in>>s.sfmx(r,c);}};
		in>>size;s.sf_tags.resize(size);for(int r = 0; r<size;r++){in>>s.sf_tags[r];};
		// set all other variables and update the object
		s.mx_init();
		s.needsUpdating = true;
		s.surfaceGen();
		s.update();
		//
		return in;
	};



};



class shp_shell
{
public:
	//public member variables
	shp_surface *srf_m;		// mid-plane (current) surface
	shp_surface *srf_u;		// undeformed surface
	shp_surface *srf_g;		// a surface that contains the gene expression data up to gL_max (at the time of loading the data).
	vector<double> sf1_vis, sf2_vis, sf3_vis;	// vectors of the gene expression fields for both visualization and energy calculation
	vector<double> sf1_2_vis, sf2_2_vis, sf3_2_vis;
	
	// generation of the cut surface for calculation of volume, self_intersection and visualization
	vector< vector<double> > Xin, Xout;			// vertex positions of outer and inner shells after cutting
	vector< vector<int> > Fc;		// the faces of the inner and outer shell cut surface (not the cut polygons of the cut plane)
	vector< vector<double> > Xcut;			// vertices of the cut-plane through the thickness of the shell
	vector< vector<int> > Fcp;	// faces of the polygons corresponding to the cut-plane through the thickness of the shell
	MatrixXd mu1_u, mu2_u;	   //contravariant components of the undeformed configuration
	MatrixXd C_o, MP_vis, MP_vis_2;				// preferred curvature of the midplane and visualization of MP region
	MatrixXd Dc;	// holds CLKs that define the thickness of the shell

	/// for fast shell_energy() calculation (initialized in prepare_energy()
	MatrixXd dgt11, dgt12, DGT1, FT,C, J,res, lambda_mx;		// 
	Eigen::MatrixXd mx;
	Eigen::MatrixXi indx_mx;

	// variables specific for the drosophila gastrulation simulation
	vector< vector<double> > Xvm; // geometry of the vittelline membrane
	MatrixXd sf1_2, sf2_2, sf3_2, MP_2; // expression fields as above, but for the second fold

	MatrixXd miu_eff, bulk_eff, kb_eff;	// matrices that get prepared in "prepare_energy_calc()"

	double active_area_fraction, active_area_fraction_vis,active_area_fraction_2, active_area_fraction_vis_2, MP1_count, MP2_count, Total_count;
	double Eb, E_nHook;
	int sfindx[6];	//indices into columns of sfmx that represent the scalar fields that we need
	double vout, vin, aout, ain;
	//MatrixXd test1, test2, test3, test4, test5;
	bool fix_normals_sign_known;
	int fix_normals_sign;
	int version;

	//public methods
	void set_new_gdim(const int D_new)
	{
		this->srf_m->set_new_gdim(D_new);
		this->srf_u->set_new_gdim(D_new);
		this->srf_g->set_new_gdim(D_new);
	}
	void set_L_max(const int L_max)
	{
		this->srf_m->set_L_max(L_max);
		this->srf_u->set_L_max(L_max);
		// set Dc L_max
		MatrixXd tmp;
		tmp = this->Dc;
		
		// 
			if(this->srf_m->xc.rows()>tmp.rows())
			{
				Dc.resize(this->srf_m->xc.rows(),1);
				//Dc.fill(0.0);
				for(int i=0;i<tmp.rows();i++){Dc(i) = tmp(i);}
			}
			if(this->srf_m->xc.rows()<tmp.rows())
			{
				Dc.resize(this->srf_m->xc.rows(),1);
				for(int i=0;i<this->srf_m->xc.rows();i++){Dc(i) = tmp(i);}
			}

	}
	void set_undeformed_shape()
	{
		/// this is a strange function -- we should replace with copy constructor and get rid of this
		/// warning: memory leak
		
		//srf_u = new shp_surface(srf_m->b->L_max, srf_m->b->dim,srf_m->sm->tri_n);
		srf_u->set_clks(srf_m->xc, srf_m->yc, srf_m->zc);// this is the undeformed shape
	}
	void update_scalar_fields()
{			
		// resize and then reset all expression fields to zero
		sf1.resize(this->srf_g->b->YLK.rows(),1);sf1.fill(0.0);
		sf2.resize(this->srf_g->b->YLK.rows(),1);sf2.fill(0.0);
		sf3.resize(this->srf_g->b->YLK.rows(),1);sf3.fill(0.0);
		sf1_2.resize(this->srf_g->b->YLK.rows(),1);sf1_2.fill(0.0);
		sf2_2.resize(this->srf_g->b->YLK.rows(),1);sf2_2.fill(0.0);
		sf3_2.resize(this->srf_g->b->YLK.rows(),1);sf3_2.fill(0.0);
		this->C_o.resize(srf_u->H.rows(), srf_u->H.cols());
		if (sfindx[0] >-1 | sfindx[1] >-1 | sfindx[2] >-1 | sfindx[3] >-1 | sfindx[4] >-1 | sfindx[5] >-1) // i.e. if any of them points to an existing field
		{
			//evaluate the gene expression fields and store them in sf1, sf2 and sf3 --- sosi: deal with index pointing to zero
			MatrixXd ones;ones.resize(this->srf_g->b->YLK.rows(),1);ones.fill(1.0);
			if (sfindx[0]>-1) {sf1 = this->srf_g->b->YLK*(srf_g->sfmx.col(this->sfindx[0]));sf1 = sf1.array()- (ones.array()* sf1.minCoeff());sf1 = sf1.array() * 3/sf1.maxCoeff();}
			if (sfindx[1]>-1) {sf2 = this->srf_g->b->YLK*(srf_g->sfmx.col(this->sfindx[1]));sf2 = sf2.array()-(ones.array()* sf2.minCoeff());sf2 = sf2.array()/sf2.maxCoeff();}
			if (sfindx[2]>-1) {sf3 = this->srf_g->b->YLK*(srf_g->sfmx.col(this->sfindx[2]));sf3 = sf3.array()-(ones.array()* sf3.minCoeff());sf3 = sf3.array()/sf3.maxCoeff();}
			if (sfindx[3]>-1) {sf1_2 = this->srf_g->b->YLK*(srf_g->sfmx.col(this->sfindx[3]));sf1_2 = sf1_2.array()-(ones.array()* sf1_2.minCoeff());sf1_2 = sf1_2.array()/sf1_2.maxCoeff();}
			if (sfindx[4]>-1) {sf2_2 = this->srf_g->b->YLK*(srf_g->sfmx.col(this->sfindx[4]));sf2_2 = sf2_2.array()-(ones.array()* sf2_2.minCoeff());sf2_2 = sf2_2.array()/sf2_2.maxCoeff();}
			if (sfindx[5]>-1) {sf3_2 = this->srf_g->b->YLK*(srf_g->sfmx.col(this->sfindx[5]));sf3_2 = sf3_2.array()-(ones.array()* sf3_2.minCoeff());sf3_2 = sf3_2.array()/sf3_2.maxCoeff();}
            // /////////////////////////////////////////////////
			// calculate relative effects of gene exression pattens
            // for most simulations we only have one field determined by psnl, and
            // psnl is generally set to 1.0. This means that facf then becomes equal to the
            // morphogen gradient determined in sf1.
            // /////////////////////////////////////////////////
			MatrixXd facf, fones;
			facf.resize(this->srf_g->b->YLK.rows(),1);
			fones.resize(this->srf_g->b->YLK.rows(),1);fones.fill(1.0);
			double lambda = 1.0;
			if (sfindx[0] >-1 | sfindx[1] >-1 | sfindx[2] >-1)  // only execute this code if one of the expression fields is defined
			{
				facf.fill(1.0);
				if (sfindx[0]>-1 & psnl!= 0.0){ 
					if (psnl<0.0){
						facf	= facf.array() * (sf1.array()* (psnl*lambda*fones.array()).array().exp());}
					else{
						facf	= facf.array() * (sf1.array().pow(psnl));
					}
				}
				if (sfindx[1]>-1 & ptw != 0.0){
					if (ptw<0.0){
						facf	= facf.array() * (sf2.array()* (ptw*lambda*fones.array()).array().exp());}
					else{
						facf	= facf.array() * (sf2.array().pow(ptw));
					}
				}
				if (sfindx[2]>-1 & phkb!= 0.0){ 
					if (phkb<0.0){
						facf	= facf.array() * (sf3.array()* (phkb*lambda*fones.array()).array().exp());}
					else{
						facf	= facf.array() * (sf3.array().pow(phkb));
					}
				}
				facf = facf/facf.norm();
			}
			else
			{ 
				facf.fill(0.0);
			}

			MatrixXd facf_2;
			facf_2.resize(this->srf_g->b->YLK.rows(),1);
			if (sfindx[3] >-1 | sfindx[4] >-1 | sfindx[5] >-1)
			{
				facf_2.fill(1.0);
				if (sfindx[3]>-1 & psf1_2!= 0.0){
					if (psf1_2<0.0){
						facf_2	= facf_2.array() * (sf1_2.array()* (psf1_2*lambda*fones.array()).array().exp());}
					else{
						facf_2	= facf_2.array() * (sf1_2.array().pow(psf1_2));
					}
				}
				if (sfindx[4]>-1 & psf2_2!= 0.0){ 
					if (psf2_2<0.0){
						facf_2	= facf_2.array() * (sf2_2.array()* (psf2_2*lambda*fones.array()).array().exp());}
					else{
						facf_2	= facf_2.array() * (sf2_2.array().pow(psf2_2));
					}
				}
				if (sfindx[5]>-1 & psf3_2!= 0.0){ 
					if (psf3_2<0.0){
						facf_2	= facf_2.array() * (sf3_2.array()* (psf3_2*lambda*fones.array()).array().exp());}
					else{
						facf_2	= facf_2.array() * (sf3_2.array().pow(psf3_2));
					}
				}
				facf_2 = facf_2/facf_2.norm();
			}
			else
			{ 
				facf_2.fill(0.0);
			}
			
			// initialize C_o according to the current model -- experiment with more options here!!
			if(this->C_o_model==0){this->C_o.fill(0.0);}		// tendency to flat tissue
			if(this->C_o_model==1){this->C_o = srf_u->H.array();} // tendency to follow the curvature of the egg
			//this->C_o = srf_u->H.array();
			//this->C_o = srf_u->H.array() + mp_co * facf.array() + mp_co_2*facf_2.array();

			// determine the MP regions for both folds and set values for C_o and miu_eff, bulk_eff, kb_eff
			
			// initialize activity maps (MP and MP2), as well as the effective miu, bulk and kb
			this->MP.resize(facf.rows(), facf.cols());MP.fill(0.0);
			this->MP_2.resize(facf_2.rows(), facf_2.cols());MP_2.fill(0.0);
			this->miu_eff.resize(facf.rows(), facf.cols());miu_eff.fill(this->miu);
			this->bulk_eff.resize(facf.rows(), facf.cols());bulk_eff.fill(this->miu);
			this->kb_eff.resize(facf.rows(), facf.cols());kb_eff.fill(this->miu);


			double maxfacf = this->mp_cutoff*facf.maxCoeff(); // sets a threshold for what is considered as active region
            double maxfacf_2 = this->mp_cutoff_2*facf_2.maxCoeff();
            
            ////////////////////////////////////////////////////////////////////////////
            ////////////// Note regarding accommodating morphogen gradients ////////////
            ////////////////////////////////////////////////////////////////////////////
			// comment out the     lines of code below to prevent thresholding of the //
            // mesoderm primordium region into an active region that is binary.       //
            ////////////////////////////////////////////////////////////////////////////
            /**/
			for(int ix=0;ix<facf.rows();ix++)
			{
                
                
                if(facf(ix,0)<=(maxfacf)){MP(ix,0) = 0.0;} else{MP(ix,0) = 1.0;}
				if(facf_2(ix,0)<=(maxfacf_2)){MP_2(ix,0) = 0.0;}else{MP_2(ix,0) = 1.0;}
				
				if(MP(ix,0)>0 & MP_2(ix,0)>0)
				{
					
					if(mp_co>mp_co_2){this->C_o(ix,0)= this->C_o(ix,0) + (mp_co * facf(ix,0) + mp_co_2* facf_2(ix,0));} // add the effects of both
					
					// use average material properties
					miu_eff(ix,0) = (this->miuMP + this->miuMP_2)/2.0;
					bulk_eff(ix,0) = (this->bulkMP + this->bulkMP_2)/2.0;
					kb_eff(ix,0) = (this->kbMP + this->kbMP_2)/2.0;

				}
				if(MP(ix,0)>0 & MP_2(ix,0)==0)
				{ 
					this->C_o(ix,0) = this->C_o(ix,0)+ mp_co* facf(ix,0);
					miu_eff(ix,0) = this->miuMP;
					bulk_eff(ix,0) = this->bulkMP;
					kb_eff(ix,0) = this->kbMP;
				}
				if(MP(ix,0)==0 & MP_2(ix,0)>0)
				{ 
					this->C_o(ix,0) = this->C_o(ix,0) +mp_co_2* facf_2(ix,0);
					miu_eff(ix,0) = this->miuMP_2;
					bulk_eff(ix,0) = this->bulkMP_2;
					kb_eff(ix,0) = this->kbMP_2;
				}
			}
            

            ////////////////////////////////////////////////////////////////////////////
            ////////////// Note regarding accommodating morphogen gradients ////////////
            ////////////////////////////////////////////////////////////////////////////
            // uncomment  the lines of code below to use morphogenetic gradient for   //
            // determining C_o,which depends on facf and C_o (the original undeformed)//
            // this only works for one field sf1 with psnl = 1
            ////////////////////////////////////////////////////////////////////////////
            /*for(int ix=0;ix<facf.rows();ix++)
            {

                    //this->C_o(ix,0) = (1-facf(ix,0)) * this->C_o(ix,0)+ mp_co* facf(ix,0);
                    this->C_o(ix,0) = this->C_o(ix,0)+ mp_co* facf(ix,0);
                    miu_eff(ix,0) = this->miuMP;
                    bulk_eff(ix,0) = this->bulkMP;
                    kb_eff(ix,0) = this->kbMP;
               
            }
            */
            

			// calculate the area fraction of the active region
			this->active_area_fraction= double((this->srf_m->b->w.array()*this->srf_m->SSn.array()*this->MP.array()).sum())/this->srf_m->A;
			this->active_area_fraction_2= double((this->srf_m->b->w.array()*this->srf_m->SSn.array()*this->MP_2.array()).sum())/this->srf_m->A;
			
			////////////////// visualization //////////////

			Eigen::Map<Eigen::VectorXd> sf1_v(&sf1_vis[0], this->sf1_vis.size()); //MatrixXd sf1_v, sf2_v, sf3_v;	// unfortunately sf1_vis is a std vector, so we need to make MatrixXd to use eigen's norm
			Eigen::Map<Eigen::VectorXd> sf2_v(&sf2_vis[0], this->sf2_vis.size());
			Eigen::Map<Eigen::VectorXd> sf3_v(&sf3_vis[0], this->sf3_vis.size());
			Eigen::Map<Eigen::VectorXd> sf1_2_v(&sf1_2_vis[0], this->sf1_2_vis.size()); 
			Eigen::Map<Eigen::VectorXd> sf2_2_v(&sf2_2_vis[0], this->sf2_2_vis.size());
			Eigen::Map<Eigen::VectorXd> sf3_2_v(&sf3_2_vis[0], this->sf3_2_vis.size());

//			MatrixXd facf_vis = ((sf1_v/sf1_v.norm()).array().pow(psnl)) * ((sf2_v/sf2_v.norm()).array().pow(ptw)) * ((sf3_v/sf3_v.norm()).array().pow(phkb));
			//MatrixXd facf_vis = ((sf1_v).array().pow(psnl)) * ((sf2_v).array().pow(ptw)) * ((sf3_v).array().pow(phkb));
			//facf_vis = facf_vis/facf_vis.norm();
			MatrixXd facf_vis;
			facf_vis.resize(sf1_v.rows(),sf1_v.cols());
			if (sfindx[0] >-1 | sfindx[1] >-1 | sfindx[2] >-1)
			{
				facf_vis.fill(1.0);
				if (sfindx[0]>-1 & psnl!= 0.0){ facf_vis	= facf_vis.array() * (sf1_v.array().pow(psnl));}
				if (sfindx[1]>-1 & ptw != 0.0){ facf_vis	= facf_vis.array() * (sf2_v.array().pow(ptw));}
				if (sfindx[2]>-1 & phkb!= 0.0){ facf_vis	= facf_vis.array() * (sf3_v.array().pow(phkb));}
				facf_vis = facf_vis/facf_vis.norm();
			}
			else
			{ 
				facf_vis.fill(0.0);
			}



			MP_vis.resize(facf_vis.size(),1);MP_vis.fill(0.0);
			double maxfacf_vis = facf_vis.maxCoeff();
			int cnt = 0;
			for(int ix=0;ix<facf_vis.size();ix++)
			{
				if(facf_vis(ix)<=(this->mp_cutoff*maxfacf_vis))	// 
				{
					MP_vis(ix,0) = 0.0;
				}
				else
				{
					MP_vis(ix,0) = 1.0;
					cnt++;
				}
			}
			// to cross reference with the Gaussian quadrature calcuation, we need the relative area here too.
			this->active_area_fraction_vis = double(cnt)/double(facf.size());
			this->MP1_count = double(cnt);
			this->Total_count = double(facf.size());

			// do the same for fold 2
			//MatrixXd facf_vis_2 = ((sf1_2_v).array().pow(psf1_2)) * ((sf2_2_v).array().pow(psf2_2)) * ((sf3_2_v).array().pow(psf3_2));
			//facf_vis_2 = facf_vis_2/facf_vis_2.norm();
			MatrixXd facf_vis_2;
			facf_vis_2.resize(sf1_2_v.rows(),sf1_2_v.cols());
			if (sfindx[3] >-1 | sfindx[4] >-1 | sfindx[5] >-1)
			{
				facf_vis_2.fill(1.0);
				if (sfindx[3]>-1 & psf1_2!= 0.0){ facf_vis_2	= facf_vis_2.array() * (sf1_2_v.array().pow(psf1_2));}
				if (sfindx[4]>-1 & psf2_2!= 0.0){ facf_vis_2	= facf_vis_2.array() * (sf2_2_v.array().pow(psf2_2));}
				if (sfindx[5]>-1 & psf3_2!= 0.0){ facf_vis_2	= facf_vis_2.array() * (sf3_2_v.array().pow(psf3_2));}
				facf_vis = facf_vis/facf_vis.norm();
			}
			else
			{ 
				facf_vis_2.fill(0.0);
			}
			
			
			MP_vis_2.resize(facf_vis_2.size(),1);MP_vis_2.fill(0.0);
			double maxfacf_vis_2 = facf_vis_2.maxCoeff();
			int cnt_2 = 0;
			for(int ix=0;ix<facf_vis_2.size();ix++)
			{
				if(facf_vis_2(ix)<=(this->mp_cutoff_2*maxfacf_vis_2))	// 
				{
					MP_vis_2(ix,0) = 0.0;
				}
				else
				{
					MP_vis_2(ix,0) = 1.0;
					cnt_2++;
				}
			}
			this->active_area_fraction_vis_2 = double(cnt_2)/double(facf_vis_2.size());
			this->MP2_count = cnt_2;
		}
	}
	void prepare_energy_calc()
	{
		// this function precedes the calculation of the shell energy
		// it calculates the following essential quantities (all are member variables of shp_shell)
		// mu1_u, mu2_u, C_o, miu_eff, bulk_eff, kb_eff
		// it calculates auxilliary quantities for visualization
		// sf1_v, sf2_v, sf3_v, active_area_fraction
		////////////////////////////////////////////////////////////////////////////////
		
		// prepare quantities needed from the undeformed shape
		// specifically we need mu1_u and mu2_u for calculating the deformation gradient field
		// and we need the preferred curvature set equal to the curvature at every point of the undeformed configuration

		

		//// prepare constraints surface (vitteline membrane)
		this->update_scalar_fields();
		this->Xin.clear();
		this->Xout.clear();
		this->srf_u->update();
		this->srf_u->surfaceGen();// update the triangular mesh surface of the midplane
		this->srf_u->sm->_curv_calc = 1;	
		this->srf_u->update_tri();// Calculate the triangular mesh properties (importantly the normals in the double**n)
		this->srf_u->sm->sfGen(this->Dc);	// generates the scalar field corresponding to the thickness

		// generate the inner and outer planes
		//double fac = 10000;
		int n_points = this->srf_u->sm->n_points;
		Xin.resize(n_points);
		Xout.resize(n_points);
		for(int i=0;i<n_points;i++)
		{
			Xout[i].resize(3);
			Xout[i][0] = this->srf_u->sm->X[i][0]-fix_normals_sign*this->srf_u->sm->vN[i][0]*this->srf_u->sm->sf[i]/2;
			Xout[i][1] = this->srf_u->sm->X[i][1]-fix_normals_sign*this->srf_u->sm->vN[i][1]*this->srf_u->sm->sf[i]/2;
			Xout[i][2] = this->srf_u->sm->X[i][2]-fix_normals_sign*this->srf_u->sm->vN[i][2]*this->srf_u->sm->sf[i]/2;
			Xin[i].resize(3);
			Xin[i][0] = this->srf_u->sm->X[i][0]+fix_normals_sign*this->srf_u->sm->vN[i][0]*this->srf_u->sm->sf[i]/2;
			Xin[i][1] = this->srf_u->sm->X[i][1]+fix_normals_sign*this->srf_u->sm->vN[i][1]*this->srf_u->sm->sf[i]/2;
			Xin[i][2] = this->srf_u->sm->X[i][2]+fix_normals_sign*this->srf_u->sm->vN[i][2]*this->srf_u->sm->sf[i]/2;
		}
		// check volumes of Xin and Xout (and determine which is outside)
		// swapping should be made more efficient in future, and can be eliminated by looking at the normals sign
		//double vout, vin, aout, ain;
		if(fix_normals_sign_known ==0)		//i.e. if we don't know if we need to test for the sign
		{
			this->srf_u->sm->tri_props(Xout,aout, vout);
			this->srf_u->sm->tri_props(Xin,ain, vin);
			
			if(vout<vin)	// then we need to swap inner and outer surfaces
			{
				fix_normals_sign = -1;
				for(int i=0;i<n_points;i++)
				{
					Xout[i][0] = this->srf_u->sm->X[i][0]-fix_normals_sign*this->srf_u->sm->vN[i][0]*this->srf_u->sm->sf[i]/2;
					Xout[i][1] = this->srf_u->sm->X[i][1]-fix_normals_sign*this->srf_u->sm->vN[i][1]*this->srf_u->sm->sf[i]/2;
					Xout[i][2] = this->srf_u->sm->X[i][2]-fix_normals_sign*this->srf_u->sm->vN[i][2]*this->srf_u->sm->sf[i]/2;
					Xin[i][0] = this->srf_u->sm->X[i][0]+fix_normals_sign*this->srf_u->sm->vN[i][0]*this->srf_u->sm->sf[i]/2;
					Xin[i][1] = this->srf_u->sm->X[i][1]+fix_normals_sign*this->srf_u->sm->vN[i][1]*this->srf_u->sm->sf[i]/2;
					Xin[i][2] = this->srf_u->sm->X[i][2]+fix_normals_sign*this->srf_u->sm->vN[i][2]*this->srf_u->sm->sf[i]/2;
				}
				/*double vtmp = vout;vout = vin;vin = vtmp;
				double atmp = aout;aout = ain;ain = atmp;
				vector< vector<double> > tmp;
				tmp.resize(n_points);
				for(int i=0;i<n_points;i++)
				{
					tmp[i].resize(3);
					tmp[i][0] = Xout[i][0];
					tmp[i][1] = Xout[i][1];
					tmp[i][2] = Xout[i][2];

					Xout[i][0] = Xin[i][0];
					Xout[i][1] = Xin[i][1];
					Xout[i][2] = Xin[i][2];

					Xin[i][0] = tmp[i][0];
					Xin[i][1] = tmp[i][1];
					Xin[i][2] = tmp[i][2];
				}*/
			}
			
		}
		fix_normals_sign_known = true;
		// prepare the morphology of the vittelline membrane
		Xvm.clear();
		Xvm.resize(Xout.size());
		for(int i = 0;i<Xout.size();i++)
		{
			Xvm[i].resize(3);
			Xvm[i][0] = this->srf_u->sm->X[i][0]-fix_normals_sign*this->srf_u->sm->vN[i][0]*(this->srf_u->sm->sf[i]/2+this->VMscale);
			Xvm[i][1] = this->srf_u->sm->X[i][1]-fix_normals_sign*this->srf_u->sm->vN[i][1]*(this->srf_u->sm->sf[i]/2+this->VMscale);
			Xvm[i][2] = this->srf_u->sm->X[i][2]-fix_normals_sign*this->srf_u->sm->vN[i][2]*(this->srf_u->sm->sf[i]/2+this->VMscale);

		}

		/// carry on with calculating the starting values

		this->update();		// among other things sets the current Xout as represented by srf_m 
		
		this->mu1_u.resize(srf_u->Xp.rows(),3);mu1_u.fill(0.0);
		this->mu2_u.resize(srf_u->Xp.rows(),3);mu2_u.fill(0.0);
		
		kk_cross(srf_u->Xp, srf_u->n, mu1_u);
		this->mu1_u.col(0) = mu1_u.col(0).array()/srf_u->SSn.array();
		this->mu1_u.col(1) = mu1_u.col(1).array()/srf_u->SSn.array();
		this->mu1_u.col(2) = mu1_u.col(2).array()/srf_u->SSn.array();
		
		kk_cross(srf_u->n, srf_u->Xt, mu2_u);
		this->mu2_u.col(0) = mu2_u.col(0).array()/srf_u->SSn.array();
		this->mu2_u.col(1) = mu2_u.col(1).array()/srf_u->SSn.array();
		this->mu2_u.col(2) = mu2_u.col(2).array()/srf_u->SSn.array();

		// prepare quantities for fast shell_energy() calculation
		dgt11.resize(this->srf_m->Xt.rows(),9);
		dgt12.resize(this->srf_m->Xp.rows(),9);
		DGT1.resize(this->srf_m->Xt.rows(),9);
		FT.resize(DGT1.rows(), DGT1.cols());
		C.resize(FT.rows(), FT.cols());
		lambda_mx.resize(C.rows(),2);

		
	}
	double shell_energy()
	{
		// before calling this function "prepare_energy_calc()" must have been called at some point
		
		// "shell_energy" function should be improved for speed, because it is the main energy calculation function for optimization
		//this->srf_m->update_fast();	// call this also before shell_energy for performing the Gaussian quadrature-based calculation to calculate
								// m_l1, m_l2 and H
		//%% Principal stretches: in-plane deformations (stretch and shear)vectorized version
		kk_kron_fast(this->srf_m->Xt, mu1_u, dgt11);
		kk_kron_fast(this->srf_m->Xp, mu2_u, dgt12);
		DGT1  = dgt11 + dgt12;
		kk_transpose_fast(DGT1, FT);
		kk_mx_mult_fast(FT, DGT1, C); //% right Cauchy-Green deformation tensor
		//lambda_mx.resize(C.rows(),2);
		//EigenSolver<Eigen::Matrix3d> eig(mx); // [vec val] = eig(A)
		for(int ix = 0;ix<int(C.rows());ix++)
		{
			
			//mx(0,0) = C(ix,0);	mx(1,0) = C(ix,1); mx(2,0) = C(ix,2);
			//mx(0,1) = C(ix,3);	mx(1,1) = C(ix,4); mx(2,1) = C(ix,5);
			//mx(0,2) = C(ix,6);	mx(1,2) = C(ix,7); mx(2,2) = C(ix,8);
			
			mx = C.row(ix);
			mx.resize(3,3);
			EigenSolver<Eigen::Matrix3d> eig(mx); // [vec val] = eig(A)	
			res = eig.pseudoEigenvalueMatrix().diagonal(); // diag(val)
			sort_eig_vec(res, indx_mx); // res gets sorted
			lambda_mx(ix,0) = res(1);
			lambda_mx(ix,1) = res(2);
		}
		lambda_mx = lambda_mx.array().sqrt();

		///////////// shear and stretch energy models based on principal stretches go here //////////////////////////////////
		/// Calculate the energy according to the model: Generalized neo-Hookean
		J = lambda_mx.col(0).array()*lambda_mx.col(1).array();
		//MatrixXd lamsq = lambda_mx.array().pow(2.0);
		//MatrixXd ones(lamsq.rows(),1);
		MatrixXd I1_bar = ((lambda_mx.col(0).array().pow(2))+ (lambda_mx.col(1).array().pow(2))+ 1).array()/J.array().pow(2.0/3.0);// % J = 1 for incompressible materials
		MatrixXd PHI = Dc(0)*(miu_eff.array()/2*(I1_bar.array()-3).array() + bulk_eff.array()/2*(J.array()-1).array().pow(2));
		//////////////////////
		E_nHook = (PHI.array()*srf_u->SSn.array()*srf_u->b->w.array()).array().sum() ; // sum over the undeformed shape
        Eb      =  (kb_eff.array()*(2*srf_m->H.array()-2*C_o.array()).array().pow(2)*srf_m->b->w.array()*srf_m->SSn.array()).array().sum();//   % calculate bending energy as in RBC paper
		
		/*
		// return memory to the system
		dgt11.resize(0,0);
		dgt12.resize(0,0);
		DGT1.resize(0,0);
		FT.resize(0,0);
		C.resize(0,0);
		J.resize(0,0);
		res.resize(0,0);
		mx.resize(0,0);
		lambda.resize(0,0);
		lamsq.resize(0,0);
		ones.resize(0,0);
		I1_bar.resize(0,0);
		PHI.resize(0,0);
		*/
		// return the energy value
		return E_nHook + Eb;
	}
	double get_energy(int flag = 0)
	{
		if(flag==1)
		{
			this->prepare_energy_calc();
			this->srf_u->update();
			this->srf_g->update();
			this->update_scalar_fields();
			this->srf_m->update_fast();
		}
		return this->shell_energy();
		//return 0.0;
	}
	void update()	// 
	{
		this->Xin.clear();
		this->Xout.clear();
		this->srf_m->update();
		this->srf_m->surfaceGen();// update the triangular mesh surface of the midplane
		this->srf_m->sm->_curv_calc = 1;	
		this->srf_m->update_tri();// Calculate the triangular mesh properties (importantly the normals in the double**n)
		this->srf_m->sm->sfGen(this->Dc);	// generates the scalar field corresponding to the thickness
		// generate the inner and outer planes
		//double fac = 10000;
		int n_points = this->srf_m->sm->n_points;
		Xin.resize(n_points);
		Xout.resize(n_points);
		for(int i=0;i<n_points;i++)
		{
			Xout[i].resize(3);
			Xout[i][0] = this->srf_m->sm->X[i][0]-fix_normals_sign*this->srf_m->sm->vN[i][0]*this->srf_m->sm->sf[i]/2;
			Xout[i][1] = this->srf_m->sm->X[i][1]-fix_normals_sign*this->srf_m->sm->vN[i][1]*this->srf_m->sm->sf[i]/2;
			Xout[i][2] = this->srf_m->sm->X[i][2]-fix_normals_sign*this->srf_m->sm->vN[i][2]*this->srf_m->sm->sf[i]/2;
			Xin[i].resize(3);
			Xin[i][0] = this->srf_m->sm->X[i][0]+fix_normals_sign*this->srf_m->sm->vN[i][0]*this->srf_m->sm->sf[i]/2;
			Xin[i][1] = this->srf_m->sm->X[i][1]+fix_normals_sign*this->srf_m->sm->vN[i][1]*this->srf_m->sm->sf[i]/2;
			Xin[i][2] = this->srf_m->sm->X[i][2]+fix_normals_sign*this->srf_m->sm->vN[i][2]*this->srf_m->sm->sf[i]/2;
		}
		// check volumes of Xin and Xout (and determine which is outside)
		// swapping should be made more efficient in future, and can be eliminated by looking at the normals sign
		//double vout, vin, aout, ain;
		if(fix_normals_sign_known ==0)		//i.e. if we don't know if we need to test for the sign
		{
			this->srf_m->sm->tri_props(Xout,aout, vout);
			this->srf_m->sm->tri_props(Xin,ain, vin);
			if(vout<vin)	// then we need to swap inner and outer surfaces
			{
				fix_normals_sign = -1;
				for(int i=0;i<n_points;i++)
				{
					Xout[i][0] = this->srf_m->sm->X[i][0]-fix_normals_sign*this->srf_m->sm->vN[i][0]*this->srf_m->sm->sf[i]/2;
					Xout[i][1] = this->srf_m->sm->X[i][1]-fix_normals_sign*this->srf_m->sm->vN[i][1]*this->srf_m->sm->sf[i]/2;
					Xout[i][2] = this->srf_m->sm->X[i][2]-fix_normals_sign*this->srf_m->sm->vN[i][2]*this->srf_m->sm->sf[i]/2;
					Xin[i][0] = this->srf_m->sm->X[i][0]+fix_normals_sign*this->srf_m->sm->vN[i][0]*this->srf_m->sm->sf[i]/2;
					Xin[i][1] = this->srf_m->sm->X[i][1]+fix_normals_sign*this->srf_m->sm->vN[i][1]*this->srf_m->sm->sf[i]/2;
					Xin[i][2] = this->srf_m->sm->X[i][2]+fix_normals_sign*this->srf_m->sm->vN[i][2]*this->srf_m->sm->sf[i]/2;
				}
				/*double vtmp = vout;vout = vin;vin = vtmp;
				double atmp = aout;aout = ain;ain = atmp;
				vector< vector<double> > tmp;
				tmp.resize(n_points);
				for(int i=0;i<n_points;i++)
				{
					tmp[i].resize(3);
					tmp[i][0] = Xout[i][0];
					tmp[i][1] = Xout[i][1];
					tmp[i][2] = Xout[i][2];

					Xout[i][0] = Xin[i][0];
					Xout[i][1] = Xin[i][1];
					Xout[i][2] = Xin[i][2];

					Xin[i][0] = tmp[i][0];
					Xin[i][1] = tmp[i][1];
					Xin[i][2] = tmp[i][2];
				}*/
			}
		}
		fix_normals_sign_known = true;
	}

	void update_fast()	// 
	{
		// assumes that prepare energy calc has been called already
		//this->Xin.clear();
		//this->Xout.clear();
		//this->srf_m->update();
		this->srf_m->surfaceGen();// update the triangular mesh surface of the midplane
		this->srf_m->sm->_curv_calc = 0;	
		this->srf_m->update_tri();// Calculate the triangular mesh properties (importantly the normals in the double**n)
		this->srf_m->sm->sfGen(this->Dc);	// generates the scalar field corresponding to the thickness
		// generate the inner and outer planes
		//double fac = 10000;
		int n_points = this->srf_m->sm->n_points;
		for(int i=0;i<n_points;i++)
		{
			Xout[i][0] = this->srf_m->sm->X[i][0]-fix_normals_sign*this->srf_m->sm->vN[i][0]*this->srf_m->sm->sf[i]/2;
			Xout[i][1] = this->srf_m->sm->X[i][1]-fix_normals_sign*this->srf_m->sm->vN[i][1]*this->srf_m->sm->sf[i]/2;
			Xout[i][2] = this->srf_m->sm->X[i][2]-fix_normals_sign*this->srf_m->sm->vN[i][2]*this->srf_m->sm->sf[i]/2;

			Xin[i][0] = this->srf_m->sm->X[i][0]+fix_normals_sign*this->srf_m->sm->vN[i][0]*this->srf_m->sm->sf[i]/2;
			Xin[i][1] = this->srf_m->sm->X[i][1]+fix_normals_sign*this->srf_m->sm->vN[i][1]*this->srf_m->sm->sf[i]/2;
			Xin[i][2] = this->srf_m->sm->X[i][2]+fix_normals_sign*this->srf_m->sm->vN[i][2]*this->srf_m->sm->sf[i]/2;
		}
		
	}

	void update_cut_plane()
	{
		double yval = 0.0;
		this->Xin.clear();
		this->Xout.clear();
		this->Fc.clear();
		this->Fcp.clear();
		this->Xcut.clear();
		

		this->srf_m->update();
		this->srf_m->surfaceGen();// update the triangular mesh surface of the midplane
		this->srf_m->sm->_curv_calc = 1;	
		this->srf_m->update_tri();// Calculate the triangular mesh properties (importantly the normals in the double**n)
		// calculate thickness (as a scalar field) at all the points defined in the spherical mesh object
		this->srf_m->sm->sfGen(this->Dc);	// generates the scalar field corresponding to the thickness
		// Generate the cut-plane surface
		tri_cut_plane(this->srf_m->sm->X, this->srf_m->sm->n_points, this->srf_m->sm->f0, this->srf_m->sm->f1, this->srf_m->sm->f2,this->srf_m->sm->n_faces, this->srf_m->sm->sf, this->srf_m->sm->vN, yval, Xin,Xout,Xcut,Fc,Fcp);
		
		// check volumes of Xin and Xout (and determine which is outside)
		// swapping should be made more efficient in future, and can be eliminated by looking at the normals sign

		if(fix_normals_sign==-1)	// then we need to swap inner and outer surfaces
		{
			vector< vector<double> > tmp;
			tmp.resize(this->srf_m->sm->n_points);
			for(int i=0;i<this->srf_m->sm->n_points;i++)
			{
				tmp[i].resize(3);
				tmp[i][0] = Xout[i][0];
				tmp[i][1] = Xout[i][1];
				tmp[i][2] = Xout[i][2];

				Xout[i][0] = Xin[i][0];
				Xout[i][1] = Xin[i][1];
				Xout[i][2] = Xin[i][2];

				Xin[i][0] = tmp[i][0];
				Xin[i][1] = tmp[i][1];
				Xin[i][2] = tmp[i][2];
			}
		}

	}

	int VM_intersect()
	{ // buggy use vtk instead
		// intersection test between Xout and Xvm
		double* V0;
		double* V1;
		double* V2;
		double* U0;
		double* U1;
		double* U2;
		for(int i = 0;i<this->srf_m->sm->n_faces;i++)
		{
			V0 = (&Xout[int(this->srf_m->sm->f0[i]-1)][0]);
			V1 = (&Xout[int(this->srf_m->sm->f1[i]-1)][0]);
			V2 = (&Xout[int(this->srf_m->sm->f2[i]-1)][0]);
			for(int j = 0;j<this->srf_m->sm->n_faces;j++)
			{
				U0 = (&Xvm[int(this->srf_m->sm->f0[j]-1)][0]);
				U1 = (&Xvm[int(this->srf_m->sm->f1[j]-1)][0]);
				U2 = (&Xvm[int(this->srf_m->sm->f2[j]-1)][0]);
				if(tri_tri_intersect(V0,V1,V2,U0,U1,U2))
				{
					int sv = share_vert((this->srf_m->sm->f0[i]-1),(this->srf_m->sm->f1[i]-1),(this->srf_m->sm->f2[i]-1),(this->srf_m->sm->f0[j]-1),(this->srf_m->sm->f1[j]-1),(this->srf_m->sm->f2[j]-1));
					if(sv==0)
					{
						return 1;
					}
				}
			}
		}
		return 0;
	}
	int self_intersect()
	{
		if(self_intersect_individual(this->Xin) == 0 && self_intersect_individual(this->Xout)==0)
		{
			return 0;
		}else
		{
			return 1;
		}
	}
	int self_intersect_individual(vector< vector<double> >X)
	{  ///buggy use vtk instead
		double* V0;
		double* V1;
		double* V2;
		double* U0;
		double* U1;
		double* U2;
		for(int i = 0;i<this->srf_m->sm->n_faces;i++)
		{
			V0 = (&X[int(this->srf_m->sm->f0[i]-1)][0]);
			V1 = (&X[int(this->srf_m->sm->f1[i]-1)][0]);
			V2 = (&X[int(this->srf_m->sm->f2[i]-1)][0]);
			for(int j = i+1;j<this->srf_m->sm->n_faces;j++)
			{
				U0 = (&X[int(this->srf_m->sm->f0[j]-1)][0]);
				U1 = (&X[int(this->srf_m->sm->f1[j]-1)][0]);
				U2 = (&X[int(this->srf_m->sm->f2[j]-1)][0]);
				if(tri_tri_intersect(V0,V1,V2,U0,U1,U2))
				{
					int sv = share_vert((this->srf_m->sm->f0[i]-1),(this->srf_m->sm->f1[i]-1),(this->srf_m->sm->f2[i]-1),(this->srf_m->sm->f0[j]-1),(this->srf_m->sm->f1[j]-1),(this->srf_m->sm->f2[j]-1));
					if(sv==0)
					{
						return 1;
					}
				}
			}
		}
		return 0;
	}
	// get and set
	double get_Young(){return Young;}
	double get_Poisson(){return Poiss;}
	double get_D(){return Dc(0);}
	double get_mp_co(){return mp_co;}
	double get_mp_co_2(){return mp_co_2;}
	double get_mp_fac(){return mp_fac;}
	double get_mp_fac_2(){return mp_fac_2;}
	double get_mp_cutoff(){return mp_cutoff;}
	double get_mp_cutoff_2(){return mp_cutoff_2;}
	double get_VMscale(){return VMscale;}
	int    get_gL_max(){return gL_max;}
	double get_psnl(){return psnl;};
	double get_ptw(){return ptw;};
	double get_phkb(){return phkb;};

	double get_psf1_2(){return psf1_2;}
	double get_psf2_2(){return psf2_2;}
	double get_psf3_2(){return psf3_2;}

	void set_Young(double Y)
	{
		this->Young = Y;
		this->YoungMP = this->mp_fac * this->Young;
		this->YoungMP_2 = this->mp_fac_2 * this->Young;
		this->lambda = this->Poiss*this->Young/(1.0-2.0*this->Poiss)/(1+this->Poiss);
		this->miu = this->Young/2.0/(1.0+this->Poiss);
		this->bulk = this->Young/3.0/(1.0-2.0*this->Poiss);
		this->kb = this->Dc(0)*this->Dc(0)*this->Dc(0)*this->Young/12.0/(1.0- (this->Poiss) * (this->Poiss));//phys.kb = phys.Young * phys.D^3/12/(1-phys.Poiss^2); %
		PoissMP    = Poiss;
		PoissMP_2  = Poiss;
		lambdaMP   = PoissMP*YoungMP/(1.0-2.0*PoissMP)/(1.0+PoissMP);
		lambdaMP_2 = PoissMP_2*YoungMP_2/(1.0-2.0*PoissMP_2)/(1.0+PoissMP_2);

		miuMP      = YoungMP/2.0/(1.0+PoissMP);
		bulkMP     = YoungMP/3.0/(1.0-2.0*PoissMP);
		kbMP       = mp_fac*kb;

		miuMP_2    = YoungMP_2/2.0/(1.0+PoissMP_2);
		bulkMP_2   = YoungMP_2/3.0/(1.0-2.0*PoissMP_2);
		kbMP_2     = mp_fac_2*kb;
	}
	void set_Poisson(double P)
	{
		this->Poiss = P;
		this->lambda = this->Poiss*this->Young/(1.0-2.0*this->Poiss)/(1+this->Poiss);
		this->miu = this->Young/2.0/(1.0+this->Poiss);
		this->bulk = this->Young/3.0/(1.0-2.0*this->Poiss);
		this->kb = this->Dc(0)*this->Dc(0)*this->Dc(0)*this->Young/12.0/(1.0- (this->Poiss) * (this->Poiss));
		PoissMP    = Poiss;
		PoissMP_2  = Poiss;
		lambdaMP   = PoissMP*YoungMP/(1.0-2.0*PoissMP)/(1.0+PoissMP);
		lambdaMP_2 = PoissMP_2*YoungMP_2/(1.0-2.0*PoissMP_2)/(1.0+PoissMP_2);

		miuMP      = YoungMP/2.0/(1.0+PoissMP);
		bulkMP     = YoungMP/3.0/(1.0-2.0*PoissMP);
		kbMP       = mp_fac*kb;

		miuMP_2    = YoungMP_2/2.0/(1.0+PoissMP_2);
		bulkMP_2   = YoungMP_2/3.0/(1.0-2.0*PoissMP_2);
		kbMP_2     = mp_fac_2*kb;

	}

	void set_D(double d)
	{
		this->Dc(0) = d;
		set_Poisson(this->Poiss);
		this->update();
	}

	// scalar fields specific
	void set_mp_co_2(double d){this->mp_co_2 = d;}
	void set_mp_fac_2(double d)
	{
		this->mp_fac_2 = d;
		this->set_Young(this->Young);
		this->set_Poisson(this->Poiss);
	}
	void set_mp_cutoff_2(double d){this->mp_cutoff_2 = d;}
	void set_psf1_2(double d){this->psf1_2 = d;}
	void set_psf2_2(double d){this->psf2_2 = d;}
	void set_psf3_2(double d){this->psf3_2 = d;}
	void set_mp_co(double d){this->mp_co = d;}
	void set_mp_fac(double d)
	{
		this->mp_fac = d;
		this->set_Young(this->Young);
		this->set_Poisson(this->Poiss);
	}
	void set_mp_cutoff(double d){this->mp_cutoff = d;}
	void set_psnl(double d){this->psnl = d;}
	void set_ptw(double d){this->ptw = d;}
	void set_phkb(double d){this->phkb = d;}
	void set_VMscale(double d){this->VMscale = d;}
	void set_gL_max(const int d)
	{
		srf_g->set_L_max(d);
		this->gL_max = d;
		/*if(this->gL_max>d)
		{
			this->srf_g->flush_after_L(d);
			this->srf_g->flush_fields_after_L(d);
			this->gL_max = d;

			shp_surface* tmp = this->srf_g;
			srf_g = new shp_surface(d, tmp->b->dim);
			srf_g->sfmx = tmp->sfmx;	// I hope this is a deep copy
			srf_g->xc = tmp->xc;
			srf_g->yc = tmp->yc;
			srf_g->zc = tmp->zc;
			delete tmp;
		}*/
	}
	void set_C_o_model(const int d){C_o_model = d;}
	//constructor and destructor
	shp_shell(const int L_max = 6, const int dim = 30, const int tri_n = 2, const int gL = 16)
	{
		//initialize the midplane and undeformed morphologies
		version = 2000;
		srf_m = new shp_surface(L_max, dim, tri_n);
		srf_u = new shp_surface(L_max, dim, tri_n);
		gL_max = gL;
		srf_g = new shp_surface(gL_max, dim, tri_n);
		this->Young = 100;		// in kP
		this->mp_fac = 20.0;
		this->mp_fac_2 = 10.0;
		this->YoungMP = this->mp_fac * this->Young;
		this->YoungMP_2 = this->mp_fac_2 * this->Young;

		this->Dc = MatrixXd((L_max + 1)*(L_max + 1), 1);		// matrix of coefficients that determine the shell thickness as function of theta and phi
		Dc.fill(0.0);
		Dc(0) = 0.5;		// in micron
		this->set_Poisson(0.45);		// also calculates and sets the other physical parameters automatically

		psnl = 0.0;
		ptw = 0.0;
		phkb = 0.0;
		psf1_2 = 0.0;
		psf2_2 = 0.0;
		psf3_2 = 0.0;

		VMscale = 5.0;
		sfindx[0] = -1;
		sfindx[1] = -1;
		sfindx[2] = -1;
		sfindx[3] = -1;
		sfindx[4] = -1;
		sfindx[5] = -1;
		mp_co = 20.0;
		mp_co_2 = 10.0;
		mp_cutoff = 0.50;
		mp_cutoff_2 = 0.80;

		active_area_fraction = 0;
		active_area_fraction_vis = 0;
		MP1_count = 0;
		MP2_count = 0;
		Total_count = 0;
		fix_normals_sign_known = 0;
		fix_normals_sign = 1;
		C_o_model = 1;
	}
	~shp_shell()
	{
		delete srf_m;
		delete srf_u;
		delete srf_g;
	}
private:
	MatrixXd sf1, sf2, sf3, MP;     // stores the actual expression field at the Gaussian weight base-points and the mesoderm primordium

	//physical variables
	double Young, YoungMP, YoungMP_2;
	double Poiss, PoissMP, PoissMP_2;
	double bulk, bulkMP, bulkMP_2; 
	double kb, kbMP, kbMP_2;
	double lambda, lambdaMP, lambdaMP_2;		// Lame's first parameter% phys.lambda = 577;   % N m^-2
	double miu, miuMP, miuMP_2;			// Shear modulus, Lame's second parameter (also G sometimes)% phys.miu = 385;      % N m^-2

	double mp_fac, mp_fac_2;				// active region stiffness factor
	double mp_co, mp_co_2;				// preferred curvature of the mp region
	double mp_cutoff, mp_cutoff_2;			// factor for thresholding the MP region (thereby controlling its area fraction)
	double psnl, ptw, phkb;		// parameters that specify the weighting of gene expressions
	double psf1_2, psf2_2, psf3_2; // parameters that specify the weighting of scalar fields for fold 2

	double VMscale;				// scale factor for the spacing between embryo tissue and vittelline membrane
	int gL_max,C_o_model;
public:
	ostream& disp(ostream& os, const shp_shell& s)
	{
		
		os<<" ------------------------------------------"<<std::endl;
		os<<" ------------- SHP SHELL OBJECT -----------"<<std::endl;
		os<<"Physical Parameters:------------------"<<std::endl;
		os<<"\tYoung's modulus: "<<s.Young<<std::endl;
		os<<"\tPoisson's ratio: "<<s.Poiss<<std::endl;
		os<<"\tlambda (Lame's first parameter): "<<s.lambda<<std::endl;
		os<<"\tmiu (shear modulus, Lame's second parameter): "<<s.miu<<std::endl;
		os<<"\tbulk modulus: "<<s.bulk<<std::endl;
		os<<"\tbending modulus: "<<s.kb<<std::endl;
		os<<"\tAverage thickness: "<<s.Dc(0)<<std::endl;
		os<<"Other Parameters:------------------"<<std::endl;
		os<<"\tSF1 exponent: "<<s.psnl<<std::endl;
		os<<"\tSF2 exponent: "<<s.ptw<<std::endl;
		os<<"\tSF3 exponent: "<<s.phkb<<std::endl;
		os<<"\tActive region threshold: "<<s.mp_cutoff<<std::endl;
		os<<"\tActive region stiffness factor: "<<s.mp_fac<<std::endl;
		os<<"\tSpacing factor for outer constraint on shell (>1.0): "<<s.VMscale<<std::endl;
		os<<"Geometric properties:----------------------"<<std::endl;
		os<<"\tmid-plane Area: "<<s.srf_m->A<<std::endl;
		os<<"\tmid-plane enclosed volume: "<<s.srf_m->V<<std::endl;
		os<<"\tTotal Gaussian curvature /(4*pi): "<<s.srf_m->T<<std::endl;
		return os;
	};


friend ostream& operator<<(ostream& os, const shp_shell& s)
	{
		
		os<<"version"<<std::endl;os<<s.version<<std::endl;
		os<<"midplane_surface>>>>>>>>>>>>>>>>>>>"<<std::endl;os<<*(s.srf_m)<<std::endl;
		os<<"undeformed_surface>>>>>>>>>>>>>>>>>>>"<<std::endl;os<<*(s.srf_u)<<std::endl;
		os<<"expression_field_surface>>>>>>>>>>>>>>>>>>>"<<std::endl;os<<*(s.srf_g)<<std::endl;

		os<<s.Young<<std::endl;
		os<<s.YoungMP<<std::endl;
		os<<s.Poiss<<std::endl;
		os<<s.PoissMP<<std::endl;
		os<<s.bulk<<std::endl;
		os<<s.bulkMP<<std::endl;
		os<<s.kb<<std::endl;
		os<<s.kbMP<<std::endl;
		os<<s.lambda<<std::endl;
		os<<s.lambdaMP<<std::endl;
		os<<s.miu<<std::endl;
		os<<s.miuMP<<std::endl;
		os<<s.mp_fac<<std::endl;
		os<<s.mp_co<<std::endl;
		os<<s.mp_cutoff<<std::endl;
		os<<s.psnl<<std::endl;
		os<<s.ptw<<std::endl;
		os<<s.phkb<<std::endl;
		os<<s.VMscale<<std::endl;
		os<<s.gL_max<<std::endl;
		os<<s.sfindx[0]<<std::endl;
		os<<s.sfindx[1]<<std::endl;
		os<<s.sfindx[2]<<std::endl;

		// for the second fold
		os<<s.YoungMP_2<<std::endl;
		os<<s.PoissMP_2<<std::endl;
		os<<s.bulkMP_2<<std::endl;
		os<<s.kbMP_2<<std::endl;
		os<<s.lambdaMP_2<<std::endl;
		os<<s.miuMP_2<<std::endl;
		os<<s.mp_fac_2<<std::endl;
		os<<s.mp_co_2<<std::endl;
		os<<s.mp_cutoff_2<<std::endl;
		os<<s.psf1_2<<std::endl;
		os<<s.psf2_2<<std::endl;
		os<<s.psf3_2<<std::endl;
		
		os<<s.sfindx[3]<<std::endl;
		os<<s.sfindx[4]<<std::endl;
		os<<s.sfindx[5]<<std::endl;

		//now export matrices/vectors
		os<<"Thickness_field_>>>>>>>>>>>>>>>>>>>"<<std::endl;
		os<<s.Dc.rows()<<std::endl;
		os<<s.Dc.cols()<<std::endl;
		for(int r = 0; r<s.Dc.rows();r++)
		{
			for(int c = 0;c<s.Dc.cols();c++)
			{
				os<<double(s.Dc(r,c))<<std::endl;
			}
		};
		//
		return os;
	};


friend istream& operator>>(istream& in, shp_shell& s)
	{
		std::string tstr = "";
		int version = 0;
		in>>tstr;in>>version;
		in>>tstr;in>>*(s.srf_m);
		in>>tstr;in>>*(s.srf_u);
		in>>tstr;in>>*(s.srf_g);

		/// read the essential variables
		in>>s.Young;
		in>>s.YoungMP ;
		in>>s.Poiss ;
		in>>s.PoissMP ;
		in>>s.bulk ;
		in>>s.bulkMP ;
		in>>s.kb ;
		in>>s.kbMP ;
		in>>s.lambda ;
		in>>s.lambdaMP ;
		in>>s.miu ;
		in>>s.miuMP ;
		in>>s.mp_fac ;
		in>>s.mp_co ;
		in>>s.mp_cutoff ;
		in>>s.psnl ;
		in>>s.ptw ;
		in>>s.phkb ;
		in>>s.VMscale ;
		in>>s.gL_max ;

		in>>s.sfindx[0];
		in>>s.sfindx[1];
		in>>s.sfindx[2];

		// for the second fold
		if(version>1000)
		{
			in>>s.YoungMP_2;
			in>>s.PoissMP_2;
			in>>s.bulkMP_2;
			in>>s.kbMP_2;
			in>>s.lambdaMP_2;
			in>>s.miuMP_2;
			in>>s.mp_fac_2;
			in>>s.mp_co_2;
			in>>s.mp_cutoff_2;
			in>>s.psf1_2;
			in>>s.psf2_2;
			in>>s.psf3_2;
			
			in>>s.sfindx[3];
			in>>s.sfindx[4];
			in>>s.sfindx[5];
		}


		//now import matrices/vectors
		int rows = 0;
		int cols = 0;
		in>>tstr;
		in>>rows;
		in>>cols;
		s.Dc.resize(rows,cols);
		double val = 0.0;
		for(int r = 0; r<rows;r++)
		{
			for(int c = 0;c<cols;c++)
			{
				in>>val;
				s.Dc(r,c) = val;
			}
		}

		// update the shell object accordingly
		s.update();
		s.update_scalar_fields();

		//
		return in;
	};



};



//////////////////////////////////////////
