/* omit map calculation by direct density optimization using clipper ccp4 and NLopt */
/* Copyright 2018 Masato Yoshimura  all rights reserved */

#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/core/coords.h>
#include <clipper/core/container.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <math.h>
#include <nlopt.h>
#include <nlopt.hpp>


class OmitCoordinates {
 public:
  OmitCoordinates() {}
  OmitCoordinates( clipper::Spacegroup spgr, clipper::Cell cell, int nomit );
  const std::vector<clipper::Coord_orth>& coordinates() { return coords; }
  double radius() const { return rad; }
 private:
  std::vector<clipper::Coord_orth> coords;
  double rad;
};


class MapFilterFn_smooth : public clipper::MapFilterFn_base {
public:
  MapFilterFn_smooth( const clipper::ftype& r0, const clipper::ftype& r1 ) : r0_(r0), r1_(r1) {}
    clipper::ftype operator() ( const clipper::ftype& radius ) const 
    {
      if ( radius > r1_ ) return 0.0;
      if ( radius < r0_ ) return 1.0;
      clipper::ftype x = ( radius - r0_ ) / ( r1_ - r0_ );
      return (2.0*x*x*x-3.0*x*x+1.0);
    }
private:
    clipper::ftype r0_, r1_;
};


OmitCoordinates::OmitCoordinates( clipper::Spacegroup spgr, clipper::Cell cell, int nomit )
{
  typedef clipper::Xmap<char>::Map_reference_index MRI;

  // calculate sphere radius
  double vasu = cell.volume() / spgr.num_symops();
  //  double r0 = 1.10*pow( vasu/double(nomit), 0.333 );
  double r0 = 1.00*pow( vasu/double(nomit), 0.333 );

  // make a map
  clipper::Resolution reso( 1.0 );
  clipper::Grid_sampling grid( spgr, cell, reso );
  clipper::Xmap<float> xmap( spgr, cell, grid );
  xmap = 0.0;

  clipper::Grid_range gm( cell, grid, r0 );
  clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw ;

  // now start packing spheres in ASU
  double cut2 = r0*r0;
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    if ( xmap[ix] == 0.0 ) {
      clipper::Coord_orth x0 = ix.coord_orth();
      coords.push_back( x0 );
      clipper::Coord_grid cent = x0.coord_frac( cell ).coord_grid( grid );
      clipper::Coord_grid g0 = cent + gm.min();
      clipper::Coord_grid g1 = cent + gm.max();
      i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
	for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    double r2 = ( iw.coord_orth() - x0 ).lengthsq();
	    if ( r2 < cut2 ) {
	      float r = 1.0-r2/cut2;
	      xmap[iw] = clipper::Util::max( xmap[iw], r );
	    }
	  }
    }
  }
  std::cout << " find the low coords.size() " << coords.size() << std::endl;
  // find the lowest points remaining
  float fcut = 0.40;
  for ( int i = 0; i < 2*nomit; i++ ) {
    MRI iz = xmap.first();
    for ( MRI iy = xmap.first(); !iy.last(); iy.next() )
      if ( xmap[iy] < xmap[iz] ) iz = iy;
    if ( xmap[iz] < fcut ) {
      clipper::Coord_orth x0 = iz.coord_orth();
      coords.push_back( x0 );
      clipper::Coord_grid cent = x0.coord_frac( cell ).coord_grid( grid );
      clipper::Coord_grid g0 = cent + gm.min();
      clipper::Coord_grid g1 = cent + gm.max();
      i0 = clipper::Xmap_base::Map_reference_coord( xmap, g0 );
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
	for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
	  for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    double r2 = ( iw.coord_orth() - x0 ).lengthsq();
	    if ( r2 < cut2 ) {
	      float r = 1.0-r2/cut2;
	      xmap[iw] = clipper::Util::max( xmap[iw], r );
	    }
	  }
    } else {
      break;
    }
  }


  // check the map for unfilled points
  int n0, n1, n2;
  n0 = n1 = n2 = 0;
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    if      ( xmap[ix] < 0.01 ) n0++;
    else if ( xmap[ix] < fcut ) n1++;
    else                        n2++;
  }
  //std::cout << coords.size() << "\n";
  //std::cout << n0 << "\t" << n1 << "\t" << n2 << "\n";
  //for ( int i = 0; i < coords.size(); i++ ) std::cout << i << " " << coords[i].format() << "\n";

  rad = r0*sqrt(1.0-fcut);
}

struct function_params
{

  clipper::HKL_data<clipper::data32::F_sigF> fo;
  clipper::Xmap<float> map;
  clipper::Xmap<float> mask;
  clipper::MTZcrystal cxtl;
  clipper::HKL_info hkls;
  clipper::Vec3<> *fxyz;
  float rf;
  //  std::vector<double> bulks;
  clipper::Coord_grid mcmin;
  clipper::Coord_grid mcmax;
  double maxderiv;
};

int count = 0;
std::vector<double> bulkscale(1500000,1.);
std::vector<double> sigmaweight(1500000,1.);
double lamda = 0.; 
bool rayment = false;
int  outfreq = 2000;
double 
rfactors (const std::vector<double> &v, std::vector<double> &grad, void *p )
{
  ++count;

  struct function_params * params = (struct function_params *)p;

  using clipper::data32::F_sigF; using clipper::data32::F_phi;

  clipper::Xmap<float> xmap = (params->map);
  clipper::Xmap<float> mask = (params->mask);
  clipper::MTZcrystal cxtl = (params->cxtl);
  clipper::HKL_info hkls = (params->hkls);
  clipper::HKL_data<F_sigF> fo = (params->fo);
  clipper::Vec3<> *fxyz = (params->fxyz);
  clipper::Resolution reso = hkls.resolution();
  float resolimit = reso.limit();
  clipper::Coord_grid gg0 = (params->mcmin);
  clipper::Coord_grid gg1 = (params->mcmax);
  double maxderiv = (params->maxderiv);

  clipper::HKL_data<clipper::data32::F_phi> fc( hkls, cxtl );
  typedef clipper::Xmap<float>::Map_reference_index MRI;
  typedef clipper::HKL_data_base::HKL_reference_index HRI;



  int i = 0;
  double sumrho = 0.;
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
    {
      if(mask[ix] > 0.0 ) {
	xmap[ix] = v[i];
        sumrho += xmap[ix];
	//        int u  = ix.coord().u();
        // int v  = ix.coord().v();
        // int w  = ix.coord().w();

        i++;
      }

    }
  

  double sumderiv = 0.;
  clipper::Xmap_base::Map_reference_coord j0, ju, jv, jw ;
  j0 = clipper::Xmap_base::Map_reference_coord( xmap, gg0);
  for ( ju = j0; ju.coord().u() < gg1.u(); ju.next_u() )
        for ( jv = ju; jv.coord().v() < gg1.v(); jv.next_v() )
          for ( jw = jv; jw.coord().w() < gg1.w(); jw.next_w() ){
	    sumderiv += fabs(xmap[jw]-xmap[jw.next_u()]) ;  
	    sumderiv += fabs(xmap[jw]-xmap[jw.next_v()]) ;  
	    sumderiv += fabs(xmap[jw]-xmap[jw.next_w()]) ;  
            //if(abs(xmap[jw]-xmap[jw.next_u()]) > maxderiv )sumderiv += 1.0 ;  
            //if(abs(xmap[jw]-xmap[jw.next_v()]) > maxderiv )sumderiv += 1.0 ;  
            //if(abs(xmap[jw]-xmap[jw.next_w()]) > maxderiv )sumderiv += 1.0 ;  
	  }
  
  
  //  if(count%1000 == 1) printf("sum deriv %f \n",sumderiv);

  int ngrid = i++;
  
  xmap.fft_to( fc  );




  int j , j2;
  j = j2 = 0;



  
   for ( HRI ih = fo.first(); !ih.last(); ih.next() ){
     if ( !fo[ih].missing() ) {
      fc[ih].f() = fc[ih].f()*bulkscale[j];
      j++;
     }}
   
  

  //    printf("fo %d fc %d  reso %f \n", j,j2,resolimit);
  int n_param = 20;

 

  std::vector<double> para(n_param, 1.0 );
  clipper::BasisFn_spline basisfn( fo, n_param, 1.0 );
  clipper::TargetFn_scaleF1F2<F_phi,F_sigF> targetfn( fc, fo );
  clipper::ResolutionFn rfn( hkls, basisfn, targetfn, para );
  double r1w, f1w, r1f, f1f, Fo, Fc, rf, fval ;
  r1w = f1w = r1f = f1f = rf = fval = 0.0;
  int ii = 0;
  int jj = 0;

  double scalefirst = sqrt( rfn.f(fo.first()) );

    fval += sumderiv * 100000. ;
  for ( HRI ih = fo.first(); !ih.last(); ih.next() )
    if ( !fo[ih].missing() ) {
       

      Fo = fo[ih].f();
      double sigF = fo[ih].sigf();
      double scalef = sqrt( rfn.f(ih) );
      Fc = scalef * fc[ih].f();


        r1w  += fabs( Fo - Fc );
	//	fval +=  ( Fo - Fc )*( Fo - Fc ) *sigmaweight[ii];
        double dif2 = ( Fo - Fc )*( Fo - Fc );
        
        double raymw;
	if (rayment ){
          raymw = exp(-fabs( Fo - Fc )/Fo);
	}else{
	  raymw = 1.0;
	}
	//	fval +=  dif2/(1. + 0.5*dif2)*sigmaweight[ii] ;
	//fval +=  dif2*raymw*sigmaweight[ii] ;
	//fval +=  dif2*raymw*sigmaweight[ii]/sigF ; 
        fval +=  dif2*raymw ;

	//      fval += dif2;
        f1w += Fo;
        ii++;

	//        jj++;
      

    }

  //  if(count == 1 )printf("reject outlier %d number of reflection used %d  hireso j2 = %d  sum fc_hireso %f\n",jj,ii,j2,sumderiv);

  //    fval = r1w;
    rf = r1w / clipper::Util::max( f1w, 0.1 );
    /// std::cout << "rfactor " << rf << "\n";   
  (params->rf) = rf;  

  //  if (count%1000 == 1)printf("count = %d  f = %8.5e rfac = %8.4f \n",count,r1w,rf);

  if (count%outfreq == 1)printf("count = %d scale[1]= %f f = %8.5e nummax deriv %f rfac = %8.4f x1 %8.4f %8.3f %8.3f %8.3f  \n",count,scalefirst,r1w,sumderiv,rf,v[0],v[1],v[2],v[3]);

  if (!grad.empty()){
    for(int i = 0 ; i< ngrid ; i++ ) {
      grad[i] = 0.;

      for( HRI ih = fo.first(); !ih.last(); ih.next() )
      if ( !fo[ih].missing() ) {
        Fo = fo[ih].f();
        double scalef = sqrt( rfn.f(ih) );
        Fc = scalef * fc[ih].f();
	///	     printf("scalef = %f invers = %f \n",scalef,ih.invresolsq() );

        if (Fc == 0.0000) Fc = 0.001;
        double FcR = scalef * fc[ih].a();
        double FcI = scalef * fc[ih].b();
        
	clipper::HKL hkls = ih.hkl();

      
          float h = hkls.h();
          float k = hkls.k();
          float l = hkls.l();

	  double fx = fxyz[i][0];
	  double fy = fxyz[i][1];
	  double fz = fxyz[i][2];


          double phi = clipper::Util::twopi() *(h*fx+k*fy+l*fz);
         
	  grad[i] += -2.*( Fo - Fc )*(FcR*cos(phi) + FcI*sin(phi))/Fc  ;

	  /*          double sign;
          sign = 1.0;
          if(Fo > Fc)sign = -1.0;
	  grad[i] +=  sign  *(FcR*cos(phi) + FcI*sin(phi))/Fc  ;
	  */


     }

      grad[i] += lamda * v[i];
      grad[i] *= sigmaweight[i];
    }

  }


  return fval + lamda * sumrho ;
}


int main( int argc, char** argv )
{
  //  CCP4Program prog( "cnlomit", "0.1.0", "$Date: 2018/01/23" );

  std::cout << "\nCopyright 2018 Masato Yoshimura \n\n";
  std::cout << " \n\n";
 

  // defaults
  clipper::String title;
  clipper::String ippdb = "NONE";
  clipper::String opmap = "NONE";
  clipper::String ipmtz = "NONE";
  clipper::String opmtz = "NONE";
  clipper::String ipcol_fo = "NONE";
  clipper::String ipcol_fc = "NONE";
  clipper::String opcol = "omit";
  clipper::String prefix = "NONE";
  bool bulk = true;
  bool sigmaw = true;
  int maxeval = 20000;
  // int nomit = 50;
  int nomit = 50;
  //  double rpad = 3.0;
  double rpad = 0.0;
  int n_refln = 200000;
  int n_param = 20;
  int verbose = 0;
  int n_free_para = 10000;
  double rx = 0.;
  double ry = 0.;
  double rz = 0.;


  // command input
  CCP4CommandInput args( argc, argv, true );
  int arg = 0;
  while ( ++arg < args.size() ) {
    if        ( args[arg] == "-title" ) {
      if ( ++arg < args.size() ) title = args[arg];
    } else if ( args[arg] == "-pdbin" ) {
      if ( ++arg < args.size() ) ippdb = args[arg];
    } else if ( args[arg] == "-mapout" ) {
      if ( ++arg < args.size() ) opmap = args[arg];
    } else if ( args[arg] == "-mtzin" ) {
      if ( ++arg < args.size() ) ipmtz = args[arg];
    } else if ( args[arg] == "-mtzout" ) {
      if ( ++arg < args.size() ) opmtz = args[arg];
    } else if ( args[arg] == "-prefix" ) {
      if ( ++arg < args.size() ) prefix = args[arg];
    } else if ( args[arg] == "-colin-fo" ) {
      if ( ++arg < args.size() ) ipcol_fo = args[arg];
    } else if ( args[arg] == "-colin-fc" ) {
      if ( ++arg < args.size() ) ipcol_fc = args[arg];
    } else if ( args[arg] == "-no-bulk" ) {
      bulk = false;
    } else if ( args[arg] == "-no-sigmaw" ) {
      sigmaw = false;
    } else if ( args[arg] == "-colout" ) {
      if ( ++arg < args.size() ) opcol = args[arg]; 
    } else if ( args[arg] == "-nomit" ) {
      if ( ++arg < args.size() ) nomit = clipper::String(args[arg]).i();
    } else if ( args[arg] == "-orth" ) {
      if ( ++arg < args.size()  ) rx = clipper::String(args[arg]).f();
      if ( ++arg < args.size()  ) ry = clipper::String(args[arg]).f();
      if ( ++arg < args.size()  ) rz = clipper::String(args[arg]).f();

    } else if ( args[arg] == "-maxeval" ) {
      if ( ++arg < args.size() ) maxeval = clipper::String(args[arg]).i();;
    } else if ( args[arg] == "-verbose" ) {
      if ( ++arg < args.size() ) verbose = clipper::String(args[arg]).i();
    } else {
      std::cout << "\nUnrecognized:\t" << args[arg] << std::endl;
      args.clear();
    }
  }
  if ( args.size() <= 1 )
    clipper::Message::message( clipper::Message_fatal( "\nUsage: cnlomit\n\t-pdbin <filename>\n\t-mtzin <filename>\n\t-mtzout <filename>\n\t-no-bulk \n\t-maxeval <maximum number of evalution to stop>\n\t-mapout <filename>\n\t-colin-fo <colpath>\n\t-nomit <number of spheres>\n.\n" ) );

  const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
  typedef clipper::HKL_info::HKL_reference_index HRI;
  typedef clipper::Xmap<float>::Map_reference_index MRI;


  if ( prefix == "NONE" ) { 


    // read MTZ
    clipper::HKL_data<clipper::data32::F_sigF> fobs;

    clipper::CCP4MTZfile mtz;
    clipper::MTZcrystal cxtl;
    clipper::HKL_info hkls;
    mtz.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
    mtz.open_read( ipmtz );
    mtz.import_hkl_info(hkls);
    mtz.import_hkl_data( fobs, ipcol_fo );
    mtz.import_crystal(  cxtl, ipcol_fo );


    if ( opcol[0] != '/' ) opcol = mtz.assigned_paths()[0].notail()+"/"+opcol;
    mtz.close_read();

    /// read PDB
    clipper::MMDBManager mmdb;
    mmdb.ReadPDBASCII((char*)ippdb.c_str());
    int hndl = mmdb.NewSelection();
    mmdb.SelectAtoms(hndl,0,0, ::mmdb::SKEY_NEW );
    clipper::mmdb::PPCAtom psel;
    int nsel;
    mmdb.GetSelIndex(hndl,psel,nsel);
    clipper::MMDBAtom_list atoms( psel, nsel ); 
    mmdb.DeleteSelection( hndl );


    clipper::HKL_data<clipper::data32::F_phi> fcpdb_bulk(hkls,cxtl);
    clipper::HKL_data<clipper::data32::F_phi> fcpdb_no(hkls,cxtl);
    clipper::HKL_data<clipper::data32::F_phi>  fphi(hkls,cxtl);
    clipper::SFcalc_obs_bulk<float> sfcb;
    sfcb(fcpdb_bulk,fobs,atoms);
    clipper::SFcalc_iso_fft<float> sfc;
    sfc(fcpdb_no,atoms);




    int i = 0;

    if ( bulk ){
    for( HRI ih = fobs.first(); !ih.last(); ih.next() ){
      if ( !fobs[ih].missing() ) {
        if ( fcpdb_bulk[ih].f() < 0.01 ) fcpdb_bulk[ih].f() = 0.01;
	bulkscale[i] = (fcpdb_bulk[ih].f()/fcpdb_no[ih].f()); 
        i++;
      }}
    }


    clipper::HKL_data<clipper::data32::F_phi> fphic( fphi );

    // set up sigmaa calc
        clipper::HKL_data<clipper::data32::Flag> modeflag( fobs );
    for ( HRI ih = modeflag.first(); !ih.last(); ih.next() )
      if ( !fobs[ih].missing() )
    	modeflag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
      else
    	modeflag[ih].flag() = clipper::SFweight_spline<float>::NONE;


      clipper::HKL_data<clipper::data32::F_phi> fb( fphi ), fd( fphi );

      clipper::HKL_data<clipper::data32::Phi_fom> phiw( fobs );
      clipper::SFweight_spline<float> sfw( n_refln, n_param );
      sfw( fb, fd, phiw, fobs, fcpdb_bulk, modeflag );

      // Is it corrrect using fom as weight? 
    bool sigmw = true;
    i = 0;
    if ( sigmw ){
    for( HRI ih = fobs.first(); !ih.last(); ih.next() ){
      if ( !fobs[ih].missing() ) {
	sigmaweight[i] = phiw[ih].fom();
        i++;
      }}
    }




    // prepare omit coords
    clipper::Spacegroup spgr = fobs.spacegroup();
    clipper::Cell       cell = fobs.cell();


    // calculate map
    //    clipper::Grid_sampling grid( spgr, cell, fobs.resolution() ); // original  
    //    clipper::Grid_sampling grid( spgr, cell, fobs.resolution(),0.781592641796 ); // Shannon factor? (3/2/Pi)**(1/3) 
        clipper::Grid_sampling grid( spgr, cell, fobs.resolution(),1.0 ); 
    //    clipper::Grid_sampling grid( spgr, cell, fobs.resolution() ); 
    int n_grid_asu = grid.size();
    clipper::Xmap<float> xmap( spgr, cell, grid );

    int n_grid_omit = n_grid_asu/nomit;

    if ( n_grid_omit > 800 ) nomit = n_grid_asu/800/2;

    std::cout << " nomit= " << nomit << " total ASU grid = " << n_grid_asu <<  std::endl; 

    clipper::Grid_range gd( fobs.cell(), grid, 3.0 );
    clipper::Xmap<float>::Map_reference_coord i0, iu, iv, iw;

    for ( int i = 0; i < atoms.size(); i++ ) if ( !atoms[i].is_null() ) {
      clipper::AtomShapeFn sf( atoms[i] );  // get atom shape fn
      clipper::Coord_frac uvw = atoms[i].coord_orth().coord_frac( hkls.cell() );
      clipper::Coord_grid g0 = uvw.coord_grid( grid ) + gd.min();
      clipper::Coord_grid g1 = uvw.coord_grid( grid ) + gd.max();
      i0 = clipper::Xmap<float>::Map_reference_coord( xmap, g0 );
      // sum all map contributions from this atoms
      for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() )
        for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() )
          for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() )
            xmap[iw] += sf.rho( iw.coord_orth() );
      }

    // now correct for multiplicity of points on special positions
    for ( clipper::Xmap<float>::Map_reference_index ix = xmap.first();
      	!ix.last(); ix.next() )
        xmap[ix] *= xmap.multiplicity( ix.coord() );







    clipper::Xmap<float> xrslt( spgr, cell, grid ), xwght( spgr, cell, grid );
    xrslt = xwght = 0.0;
    float r2 = 5.;
    MapFilterFn_smooth flt_sml( 5., 5. );
    MapFilterFn_smooth flt_lrg( 5., 5. );



     double maxderiv = 0.;
     double sumderiv = 0.;
     int nderiv=0;

    clipper::Xmap<float>::Map_reference_index iix2;
    for ( clipper::Xmap<float>::Map_reference_index iix = xmap.first();
	  !iix.last(); iix.next() ){
      int dw = 1; int du = 0; int dv = 0;
      //      clipper::Xmap<float>::Map_reference_index iix2 = iix.index_offset(du,dv,dw);
      double derivu = fabs(xmap[iix] - xmap.get_data(iix.index_offset(1,0,0)));
      if( maxderiv < derivu )maxderiv = derivu;
      double derivv = fabs(xmap[iix] - xmap.get_data(iix.index_offset(0,1,0)));
      if( maxderiv < derivv )maxderiv = derivv;
      double derivw = fabs(xmap[iix] - xmap.get_data(iix.index_offset(0,0,1)));
      if( maxderiv < derivw )maxderiv = derivw;

      sumderiv += derivu + derivv + derivw;
      nderiv ++;
    }

    printf("maxderivw = %f numderiv %d sum %f ave %f \n",maxderiv,nderiv,sumderiv,sumderiv/nderiv/3);

    maxderiv = 3.0;

    clipper::Xmap<float> xmask(spgr, cell, grid);
    xmask = 0.;


     clipper::Coord_orth x0(rx,ry,rz ) ;
     // make the mask parameters
         clipper::Grid_range gm2( cell, grid, r2 );
         clipper::Xmap_base::Map_reference_coord j0, ju, jv, jw ;
         clipper::Coord_grid cent2 = x0.coord_frac( cell ).coord_grid( grid );
         clipper::Coord_grid gg0 = cent2 + gm2.min();
         clipper::Coord_grid gg1 = cent2 + gm2.max();


          j0 = clipper::Xmap_base::Map_reference_coord( xmap, gg0 );
      int ngrid = 0;





      double aveden = 0.;

     for ( ju = j0; ju.coord().u() <= gg1.u(); ju.next_u() )
	for ( jv = ju; jv.coord().v() <= gg1.v(); jv.next_v() )
	  for ( jw = jv; jw.coord().w() <= gg1.w(); jw.next_w() ) {
	    const double r2 = ( jw.coord_orth() - x0 ).lengthsq();
	    const float wgt = ( 1.0 - flt_lrg( sqrt( r2 ) ) );

            if( wgt < 1.0) {
            xmask[jw] = 1.0;
            aveden += xmap[jw];
   
           ngrid++; }
	  }

      aveden /= ngrid;


      std::vector<double> x1(2000,aveden);

      clipper::Vec3<> fxyzf[10000];

      ngrid = 0;
      int n_all =0;
      float sumini = 0., sum2ini = 0. ;
      for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
        sumini += xmap[ix];
        sum2ini += xmap[ix]*xmap[ix];        
        n_all++;
        if ( xmask[ix] > 0. )
	  {
	       x1[ngrid] = xmap[ix];
		clipper::Coord_grid fxyz = ix.coord();

		fxyzf[ngrid][0] = fxyz.u()/(1.0*grid.nu());
		fxyzf[ngrid][1] = fxyz.v()/(1.0*grid.nv());
        	fxyzf[ngrid][2] = fxyz.w()/(1.0*grid.nw());
        	ngrid++;
	  }
      }

      float sig_ini = sum2ini/n_all - sumini *sumini /n_all/n_all;

      std::cout << "numer of grid  = "  << ngrid << " average density =  " << aveden << std::endl;


      clipper::CCP4MAPfile mskout;
      mskout.open_write( "msk.map" );
      mskout.export_xmap( xmask );
      mskout.close_write();
 



      double rf = 0.0;
      struct function_params inputp = { fobs,xmap,xmask,cxtl,hkls,fxyzf,rf,gg0,gg1,maxderiv };

      nlopt::opt opt(nlopt::LN_SBPLX,ngrid);

      // nlopt::opt local_opt(nlopt::LD_MMA,ngrid);

      // opt.set_local_optimizer(local_opt);


      int j ;
      std::vector<double> lb(ngrid);
      std::vector<double> ub(ngrid);
            for (int i = 0; i < ngrid ;i++) lb[i] = 0.0;
      opt.set_lower_bounds(lb);


            for (int i = 0; i < ngrid ;i++) ub[i] = aveden*30.;
      opt.set_upper_bounds(ub);
            std::cout << "lb0= " << lb[0] << " ub0 = " << ub[0] << " \n"; 

      opt.set_ftol_rel(1.0E-14);
      opt.set_maxtime(140000);
      opt.set_maxeval(maxeval);
          opt.set_min_objective(rfactors,&inputp);
      std::vector<double> x2(ngrid);
           for (int i = 0; i < ngrid ;i++) x2[i] = aveden;

      double rf2;
      nlopt::result result = opt.optimize(x2,rf2);

      std::cout << "result = " << result << " \n";

      size_t iter = 0;
      int status;
      double size;



      printf ("%5d  f()= %7.3f   rfac = %8.3f  count = %d  \n", iter, rf2, inputp.rf, count );

      float sumfinal = 0., sum2final = 0. ;

      float sumdiff = 0.;
      for (int i = 0 ; i < ngrid; i++) {

      float initialx = x1[i];
      float finalx = x2[i];

      sumfinal += finalx;
      sum2final += finalx*finalx;

      sumdiff += (initialx - finalx)*(initialx - finalx);

      if (i < 10 || i > ngrid-10)printf("i = %d ini x %10.3f   final x %10.3f  diff %8.5e  ave_diff %8.4e \n",initialx,finalx,sumdiff,sumdiff/(i+1.));

      }

      float sig_final = sum2final/ngrid - sumfinal*sumfinal/ngrid/ngrid;


      double rf3;
      rayment = true;
      outfreq = 1000;
      count  = 0;
      nlopt::result result2 = opt.optimize(x2,rf3);
      std::cout << "result = " << result2 << " \n";

      printf ("%5d  f()= %7.3f   rfac = %8.3f \n", iter, rf3, inputp.rf );

       sumfinal = 0., sum2final = 0. ;

        sumdiff = 0.;
      for (int i = 0 ; i < ngrid; i++) {

      float initialx = x1[i];
      float finalx = x2[i];

      sumfinal += finalx;
      sum2final += finalx*finalx;

      sumdiff += (initialx - finalx)*(initialx - finalx);

      if (i < 10 || i > ngrid-10)printf("i = %d ini x %10.3f   final x %10.3f  diff %8.5e  ave_diff %8.4e \n",initialx,finalx,sumdiff,sumdiff/(i+1.));

      }

       sig_final = sum2final/ngrid - sumfinal*sumfinal/ngrid/ngrid;







      printf("initial all sigma %10.5e   sigma at mask  %10.5e   ratio %8.3f \n",sig_ini,sig_final,sig_final/sig_ini);


      std::cout << " rfactor = " << inputp.rf << "\n";
      std::cout << " algorithm = " <<  opt.get_algorithm_name() << "\n";

      int ngrid2=0;
    for ( MRI ix = xmap.first(); !ix.last(); ix.next() )
      {
      if (xmask[ix]> 0.){
	xrslt[ix] = x2[ngrid2];
        ngrid2++;
      }else{
      xrslt[ix] = xmap[ix];
      }
      }


    // write map
    if ( opmap != "NONE" ) {
      clipper::CCP4MAPfile mapout;
      mapout.open_write( opmap );
      mapout.export_xmap( xrslt );
      mapout.close_write();
   }



    // write mtz
    if ( opmtz != "NONE" ) {
      // calculate
      xrslt.fft_to( fphi );

      // write results
      mtz.open_append( ipmtz, opmtz );
      mtz.export_hkl_data( fphi, opcol );
      //    mtz.export_hkl_data( abcd, opcol );
      mtz.close_append();
    }

  }


}
