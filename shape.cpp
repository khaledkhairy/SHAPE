#include <qapplication>
#include <qfiledialog>
#include <qmenubar>
#include <qdebug>
#include <QApplication>
#include <QFont>
#include <QPushButton>
#include <QWidget>
#include <QSlider>
#include <QtGui>

#include "shape.h"

#include <vtkActor.h>
#include <vtkScalarBarActor.h>
#include <vtkAxesActor.h>
#include <vtkCubeAxesActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include "vtkCylinderSource.h"
#include <vtkDataSetMapper.h>
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkLookupTable.h"
#include "vtkPolyDataMapper.h"
#include "vtkCurvatures.h"
#include "vtkCleanPolyData.h"
#include "vtkCamera.h"
#include "vtkOBJExporter.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPLYWriter.h"
#include "vtkSTLWriter.h"
#include "vtkTextProperty.h"

#include <iostream>
#include "shape_tools.h"


/////
void shape2polydata(const shp_surface* s, vtkPolyData* polydata);

shape::shape(QWidget *parent)
	: QMainWindow(parent)
{
	Ui_MainWindow::setupUi(this);
	//this->setFixedSize(this->size());
	_minH = 0.0;
	_maxH = 0.0;
	// QT/VTK
	_ren1 = vtkRenderer::New();
	vtkWidget1->GetRenderWindow()->AddRenderer(_ren1);
	vtkWidget1->GetRenderWindow()->SetAAFrames(3);
	if(KKDEBUG)
	{
		std::cout<<"Spherical Harmonics Parameterization Object Viewer (v.1.0 x64 release)"<<std::endl<<"Created by Khaled Khairy"<<std::endl;
		std::cout<<"Copyright 2012 Howard Hughes Medical Institute"<<std::endl<<"---------------------------------------"<<std::endl<<std::endl;
	}
	// initialize the shape object
	int tri_n = 4;
	int L_max = 24;
	int dim = 30;
	if(KKDEBUG){std::cout<<"Using basis with Lmax = "<<L_max<<" and Gaussian quadrature base points: "<<dim*dim<<std::endl;}
if(KKDEBUG){std::cout<<"Initializing spherical harmonics surface object ... ";}
	_s = new shp_surface(L_max, dim, tri_n);
if(KKDEBUG){std::cout<<"Done!"<<std::endl;}
	check_self_intersection = YES;
	cscale = 1.0/cscale_calc();
	_scale_fac = 3;
	_axis_on = true;
	scalar_bar_vis = true;

	_bgrdr = 0.0;
	_bgrdg = 0.0;
	_bgrdb = 0.0;
	// intialize the triangular mesh 
	// NOTE: the MINUS ONE because the triangulations in trifaces.h come from Matlab
if(KKDEBUG){std::cout<<"Initializing mesh faces ... ";}
	_polys = vtkCellArray::New();
	for (int i=0; i<_s->sm->n_faces; i++)
	{
		_polys->InsertNextCell(3);
		_polys->InsertCellPoint(_s->sm->f0[i]-1);
		_polys->InsertCellPoint(_s->sm->f1[i]-1);
		_polys->InsertCellPoint(_s->sm->f2[i]-1);
	}
if(KKDEBUG){std::cout<<"Done!"<<std::endl;}
if(KKDEBUG){std::cout<<"Initializing axis ... ";}
	_axes_actor	= vtkAxesActor::New();
	_axes_actor->SetShaftTypeToLine();
	_axes_actor->SetConeRadius(0.1);
	_axes_actor->SetXAxisLabelText("x");
	_axes_actor->SetYAxisLabelText("y");
	_axes_actor->SetZAxisLabelText("z");
	_axes_actor->SetTotalLength(5, 5, 5);
 _ren1->AddActor(_axes_actor);		//  to see the axis
if(KKDEBUG){std::cout<<"Done!"<<std::endl;}

// setting up the surface look and look up table
	lut= vtkLookupTable::New();
	lut->SetNumberOfTableValues(400);
    lut->SetHueRange(0.667,0.0);
    lut->SetSaturationRange(1,1);
    lut->SetValueRange(1,1);
	lut->SetRampToSCurve();

	Mapper		= vtkPolyDataMapper::New();	
	//Mapper->ScalarVisibilityOn();
	//Mapper->SetScalarModeToUsePointData();
	//Mapper->SetColorModeToMapScalars();
	Mapper->SetLookupTable(lut);

//	Mapper->SetScalarMaterialModeToAmbient();
	Actor		= vtkActor::New();
	_ren1->SetBackground(_bgrdr,_bgrdg,_bgrdb);
	cube_axis_actor1_on = false;
	_cube_axes_actor1 = vtkCubeAxesActor::New();
	
	vtkTextProperty* tprop = vtkTextProperty::New();
    tprop->SetColor(1,1,1);
    tprop->SetFontFamilyToArial();
	tprop->SetFontSize(32);

	scalarBar = vtkScalarBarActor::New();
	scalarBar->SetNumberOfLabels(4);
	//scalarBar->SetTitle("Local mean curvature");
	scalarBar->SetMaximumWidthInPixels(80);
	scalarBar->SetLabelTextProperty(tprop);
	scalarBar->SetTitleTextProperty(tprop);
	_ren1->AddActor(scalarBar);

	// update vtk view and shp_surface object
	 //update_vtkwindow();
	  on_initializeButton_clicked();

	// update gui objects
	//  synchronize_shape2sliders();
	comboBox_gdim->blockSignals(true);
	comboBox_gdim->setCurrentIndex(4);
	comboBox_gdim->blockSignals(false);
	comboBox_spherical_mesh->blockSignals(true);
	comboBox_spherical_mesh->setCurrentIndex(tri_n-1);
	comboBox_spherical_mesh->blockSignals(false);

	checkBox_self_intersection->blockSignals(true);
	checkBox_self_intersection->setChecked(true);
	checkBox_self_intersection->blockSignals(false);

	if(tri_n>2)
	{
		if(KKDEBUG){std::cout<<".... self-intersection test disabled! "<<std::endl;}
		checkBox_self_intersection->setChecked(false);
		this->checkBox_self_intersection->setEnabled(false);
		this->check_self_intersection= NO;
	}

}

shape::~shape()
{

}
void shape::on_initializeButton_clicked(){
	//_s->set_new_basis(6, _s->b->dim);
	_s->initializeToSphere();
	cscale = 1.0/cscale_calc();
	_axes_actor->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
	this->comboBox_sf->setCurrentIndex(0);
	synchronize_shape2sliders();
	update_vtkwindow();
	this->update_comboBox_sf();
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
	this->update();
}
void shape::on_axes_toggle_Button_clicked()
{
	if(_axis_on)
	{
		_axes_actor->SetVisibility(false);
		_ren1->GetRenderWindow()->Render();
		_axis_on = false;
		
	}
	else
	{
		_axes_actor->SetVisibility(true);
		_ren1->GetRenderWindow()->Render();
		_axis_on = true;
	}
	
}
void shape::on_pushButton_cube_axis_toggle_clicked()
{
	if(this->cube_axis_actor1_on){this->cube_axis_actor1_on = 0;}
	else
	{this->cube_axis_actor1_on = 1;}
	this->update_vtkwindow();
}
void shape::on_actionOpen_triggered()// Action to be taken upon file open 
{
	//QFileDialog dialog(this);
	//dialog.setFileMode(QFileDialog::ExistingFile);
	QString filename = QFileDialog::getOpenFileName(0, tr("Open shape"), QDir::currentPath(), tr("SHP Files (*.shp3)"));
	qDebug() << filename;
	
	if ( filename.isEmpty() == false )
	{
		std::string filestr = (const char*)filename.toLatin1(); // do it like this because microscoft stl and qt stl clash
		std::cout<<"Opening: "<<filestr<<std::endl;
		if(!(filestr.length()==0)){	this->import_shape_from_disc(filestr);}
	}
}
//functions to handle the file-drop event
QString shape::mimeData_2_fileName( const QMimeData * mimeData )
{
	if ( mimeData->hasUrls() )
	{
		foreach ( const QUrl & url, mimeData->urls() )
		{
			QString str = url.toLocalFile();
			if ( str.isEmpty() == false )
			{
				return str;
			}
		}
	}
	return QString();
}
void shape::dragEnterEvent( QDragEnterEvent * event )
{
	if ( mimeData_2_fileName( event->mimeData() ).isEmpty() == false )
	{
		event->acceptProposedAction();
	}
}
void shape::dropEvent( QDropEvent * event )
{
	QString filename = mimeData_2_fileName( event->mimeData() );//qDebug() << filename;
	if ( filename.isEmpty() == false )
	{
		std::string filestr = (const char*)filename.toLatin1(); // do it like this because microscoft stl and qt stl clash
		std::cout<<"Opening: "<<filestr<<std::endl;
		if(!(filestr.length()==0))	{this->import_shape_from_disc(filestr);	}
	}
}
void shape::import_shape_from_disc(std::string filestr)
{			
	_s->read_trunc(filestr.c_str());		// read_trunc causes truncation of series or filling with zeros to keep L_max the same
	_s->center_to_zero();
	//_s->xc(0) = 0.0;_s->yc(0) = 0.0;_s->zc(0) = 0.0;
	cscale = 1.0/cscale_calc();
	_axes_actor->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
	_s->needsUpdating = YES;
	if(_s->b->L_max<6){	_s->set_L_max(6);}


	synchronize_shape2sliders();
	update_comboBox_sf();
	this->comboBox_sf->setCurrentIndex(0);
	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = 0.0;}	// set the scalar field to zero
	//
	this->update();
	update_vtkwindow();
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
}
void shape::on_actionExit_triggered() {
  qApp->exit();
}

void shape::on_actionSave_triggered()
{
	QString filename = QFileDialog::getSaveFileName(vtkWidget1, "Save shape to file ...",  "*.shp3");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0)) {
		_s->write(filestr.c_str());
	}
}
// Action to be taken upon Export
void shape::on_actionExport_triggered()
{
    /*
		QString filename = QFileDialog::getSaveFileName(vtkWidget1, "",  "*.obj");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0))
	{
		int tmp = 0;
		// prepare the scene by removing the axis actor if present
		if(_axis_on)
		{
			
			_ren1->RemoveAllViewProps(); // remove the existing actors
			_ren1->ResetCamera();
			_ren1->GetRenderWindow()->Render();
			_ren1->AddActor(Actor);		//  add the axis actor
		    _ren1->ResetCamera();
			_ren1->GetRenderWindow()->Render();
		}

		vtkOBJExporter *obj = vtkOBJExporter::New();
		obj->SetInput(_ren1->GetRenderWindow());
		obj->SetFilePrefix((filename.toStdString()).c_str());
		obj->Write();
		if(_axis_on)
		{
			_ren1->AddActor(_axes_actor);		//  add the axis actor
		    _ren1->ResetCamera();
			_ren1->GetRenderWindow()->Render();
		}
	}*/
}
void shape::on_actionExport_to_STL_triggered()
{
    /*
		QString filename = QFileDialog::getSaveFileName(vtkWidget1, "",  "*.stl");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0))
	{
		vtkPolyData*		polydata	= vtkPolyData::New();
		shape2polydata(_s,polydata);
		vtkSTLWriter *w = vtkSTLWriter::New();
		w->SetInputData(polydata);
		w->SetFileName((filename.toStdString()).c_str());
		w->Write();
	}
     */
}
void shape::on_actionExport_to_PLY_triggered()
{
    /*
	QString filename = QFileDialog::getSaveFileName(vtkWidget1, "",  "*.ply");
	std::string filestr = (const char*)filename.toLatin1();
	if(!(filestr.length()==0))
	{
		vtkPolyData*		polydata	= vtkPolyData::New();
		shape2polydata(_s,polydata);
		vtkPLYWriter *ply = vtkPLYWriter::New();
		ply->SetInputData(polydata);
		ply->SetFileName((filename.toStdString()).c_str());
		ply->Write();
	}
     */
}



void shape::shape2polydata(const shp_surface *s, vtkPolyData* polydata){
// communicating the shape to vtk goes here
// surface triangulation is stored in a vtk polydata object
	vtkPoints* points = vtkPoints::New();
	points->SetDataTypeToDouble();
	points->SetNumberOfPoints(_s->sm->n_points);
	for (int i=0; i<_s->sm->n_points; i++){	points->SetPoint(i, s->sm->X[i][0], s->sm->X[i][1], s->sm->X[i][2]);}
	polydata->SetPoints(points);points->Delete();// Assign points and cells
	polydata->SetPolys(_polys);//polys->Delete();
	if(this->comboBox_sf->currentIndex())
	{
		// add the scalar field 
		vtkFloatArray *scalars = vtkFloatArray::New();
		float tmp;
		_minH = 0.0;
		_maxH = 0.0;
		for (int i=0; i<s->sm->n_points; i++){
	 		if((this->_s->sm->sf[i])>_maxH){_maxH = this->_s->sm->sf[i];}
			if((this->_s->sm->sf[i])<_minH){_minH = this->_s->sm->sf[i];}
		}
		//std::cout<<"found range"<<_minH<<" "<<_maxH<<"  ... rescaling."<<std::endl;
		double range = _maxH-_minH;
		//double lo = 0.0;
		//double hi = 250;
		//double fac = 255.0/(range);
		double fac = _maxH/range;
		double val = 999;
		for (int i=0; i<s->sm->n_points; i++)
		{
			
			//val = (this->_s->sm->sf[i]) * fac;
			//val = val-_minH + lo + (val-_minH) * (hi-lo)/range;
			val = this->_s->sm->sf[i];
			scalars->InsertNextValue(val);
		}
		polydata->GetPointData()->SetScalars(scalars);scalars->Delete();//uncomment if there are curvature values associated with polydata
	}
}

void shape::update_vtkwindow(){
	//gets called when the shp_surface object needs to be repainted
	_s->needsUpdating = true;
		_s->update();
		_s->surfaceGen();
		_s->update_tri();
		_s->needsUpdating = NO;
	this->_update_gui();

/*	_ren2->RemoveActor(_ren2->GetActors()->GetLastItem());
	_ren2->AddActor(Actor);
	//_ren->ResetCamera();
	_ren2->GetRenderWindow()->Render();
*/
	/*	//let's try to use vtk to show the curvature
	vtkCleanPolyData* cleaner = vtkCleanPolyData::New();
	cleaner->SetInput(polydata);

	vtkCurvatures* curve = vtkCurvatures::New();
	curve->SetInput(cleaner->GetOutput());
	curve->SetCurvatureTypeToMean();

	vtkLookupTable* lut2 = vtkLookupTable::New();
	lut2->SetNumberOfColors(256);
	lut2->SetHueRange(0.15,1.0);
    	lut2->SetSaturationRange (1.0,1.0);
    	lut2->SetValueRange (1.0,1.0);
    	lut2->SetAlphaRange(1.0,1.0);
    	lut2->SetRange(0,4);

	Mapper->SetInputConnection(curve->GetOutputPort());
	Mapper->SetLookupTable(lut2);
	Mapper->SetUseLookupTableScalarRange(1);

*/
	/*
	//let's try to use vtk to show the curvature
	vtkCurvatures* curve = vtkCurvatures::New();
	curve->SetInput(polydata);
	curve->SetCurvatureTypeToMean();

	vtkLookupTable* lut2 = vtkLookupTable::New();
	lut2->SetNumberOfColors(256);
	lut2->SetHueRange(0.15,1.0);
    	lut2->SetSaturationRange (1.0,1.0);
    	lut2->SetValueRange (1.0,1.0);
    	lut2->SetAlphaRange(1.0,1.0);
    	lut2->SetRange(0,4);

	Mapper->SetInputConnection(curve->GetOutputPort());
	Mapper->SetLookupTable(lut2);
	Mapper->SetUseLookupTableScalarRange(1);


	Mapper->SetInput(polydata);
	Actor->SetMapper( Mapper);
	
	//
	_ren1->RemoveActor(_ren1->GetActors()->GetLastItem());
	_ren1->AddActor(Actor);
	_ren1->GetRenderWindow()->Render();
*/


	


/*
	_ren2->RemoveActor(_ren2->GetActors()->GetLastItem());
	_ren2->AddActor(Actor);
	_ren2->GetRenderWindow()->Render();
*/	
	/*
	//let's try to use vtk to show the curvature
	vtkCurvatures* curve = vtkCurvatures::New();
	curve->SetInput(polydata);
	curve->SetCurvatureTypeToMean();

	vtkLookupTable* lut2 = vtkLookupTable::New();
	lut2->SetNumberOfColors(256);
	lut2->SetHueRange(0.15,1.0);
    	lut2->SetSaturationRange (1.0,1.0);
    	lut2->SetValueRange (1.0,1.0);
    	lut2->SetAlphaRange(1.0,1.0);
    	lut2->SetRange(0,4);

	Mapper->SetInputConnection(curve->GetOutputPort());
	Mapper->SetLookupTable(lut2);
	Mapper->SetUseLookupTableScalarRange(1);


	Mapper->SetInput(polydata);	
	Actor->SetMapper( Mapper);
	_ren1->RemoveActor(_ren1->GetActors()->GetLastItem());
	_ren1->AddActor(Actor);
	_ren1->ResetCamera();
	_ren1->GetRenderWindow()->Render();
*/
}
void shape::update_vtkwindow(int coordix, int cix, double old_val, double new_val)
{
	//gets called when the shp_surface object needs to be repainted
	if(_s->needsUpdating)
	{
		_s->update(coordix, cix, old_val, new_val);
		//_s->update();
		_s->surfaceGen(coordix, cix, old_val, new_val);
		_s->update_tri();
		_s->needsUpdating = NO;
	}
	this->_update_gui();
	/*
	this->update_sf(this->comboBox_sf->currentIndex());
	vtkPolyData*		polydata	= vtkPolyData::New();
	shape2polydata(_s,polydata);
	if(this->check_self_intersection)
	{
		if(_s->self_intersect())
		{
			std::cout<<"Self intersection detected"<<std::endl;
		}
	}

	////////////////
	//lcdNumberLmax->display(_s->b->L_max);
	lineEdit_Lmax->setText(QString::number(uint(_s->b->L_max)));
	lcdNumberArea->display(_s->A);
	lcdNumberVolume->display(_s->V);
	lcdNumberReducedVolume->display(_s->v);
	lcdNumberEb->display(_s->Eb);
	lcdNumberh->display(_s->h);
	lcdNumberWb->display(_s->T);

	lcdNumberArea_2->display(_s->sm->A);
	lcdNumberVolume_2->display(_s->sm->V);
	lcdNumberReducedVolume_2->display(_s->sm->v);
	lcdNumberEb_2->display(_s->sm->Eb);
	lcdNumberh_2->display(_s->sm->h);
	
	////////////////
	Mapper->SetInput(polydata);
	Mapper->SetScalarRange(_minH, _maxH);
	if(this->comboBox_sf->currentIndex()>9)
	{
		//Mapper->SetScalarRange(20, 250);
	}
	Actor->SetMapper( Mapper);
	_ren1->RemoveActor(_ren1->GetActors()->GetLastItem());
	_ren1->AddActor(Actor);
		//sosi
	if(this->scalar_bar_vis)
	{
		if(this->comboBox_sf->currentIndex())
		{
			
			scalarBar->SetLookupTable(Mapper->GetLookupTable());
			scalarBar->SetVisibility(true);
			_ren1->GetRenderWindow()->Render();
			//scalarBar->SetTitle("Scalar field");
		}
		else
		{
			scalarBar->SetVisibility(false);
		}
	}
	_ren1->GetRenderWindow()->Render();
	*/
}
void shape::_update_gui()
{
	this->update_sf(this->comboBox_sf->currentIndex());
	vtkPolyData*		polydata	= vtkPolyData::New();
	shape2polydata(_s,polydata);
	if(this->check_self_intersection)
	{
		if(_s->self_intersect())
		{
			std::cout<<"Self intersection detected"<<std::endl;
		}
	}

	////////////////
	//lcdNumberLmax->display(_s->b->L_max);
	lineEdit_Lmax->setText(QString::number(uint(_s->b->L_max)));
	lcdNumberArea->display(_s->A);
	lcdNumberVolume->display(_s->V);
	lcdNumberReducedVolume->display(_s->v);
	lcdNumberEb->display(_s->Eb);
	lcdNumberh->display(_s->h);
	lcdNumberWb->display(_s->T);

	lcdNumberArea_2->display(_s->sm->A);
	lcdNumberVolume_2->display(_s->sm->V);
	lcdNumberReducedVolume_2->display(_s->sm->v);
	lcdNumberEb_2->display(_s->sm->Eb);
	lcdNumberh_2->display(_s->sm->h);
	
	////////////////
	//cout<<"minimum scalar field value: "<<_minH<<std::endl;
	//cout<<"maximum scalar field value: "<<_maxH<<std::endl;
	
	//this->lut->SetTableRange(_minH, _maxH);
	//this->lut->Build();
	Mapper->SetInputData(polydata);
	Mapper->SetScalarRange(_minH, _maxH);
	//Mapper->SetLookupTable(lut);
	Actor->SetMapper( Mapper);
	_ren1->RemoveActor(_ren1->GetActors()->GetLastItem());
	_ren1->AddActor(Actor);
	if(this->cube_axis_actor1_on)
	{
		this->_cube_axes_actor1->SetBounds(Actor->GetBounds());
		this->_cube_axes_actor1->SetCamera(_ren1->GetActiveCamera());
		_ren1->AddActor(this->_cube_axes_actor1);
	}
	if(this->scalar_bar_vis)
	{
		if(this->comboBox_sf->currentIndex())
		{
			double range[2];
			polydata->GetPointData()->GetScalars()->GetRange(range);
			//cout<<range[0]<<std::endl<<range[1]<<std::endl;
			lut->SetRange(range);
			lut->Build();
			scalarBar->SetLookupTable(Mapper->GetLookupTable());
			//scalarBar->SetLookupTable(lut);
			scalarBar->SetVisibility(true);
			//_ren1->GetRenderWindow()->Render();
			//scalarBar->SetTitle("Scalar field");
		}
		else
		{
			scalarBar->SetVisibility(false);
		}
	}

	_ren1->GetRenderWindow()->Render();
}
void shape::on_lineEdit_Lmax_returnPressed()
{
	const QString text  = this->lineEdit_Lmax->text();
	int Lmax_new = text.toInt();
	if(Lmax_new<6){Lmax_new = 6;};
	if(Lmax_new>80){Lmax_new = 80;};
	if(KKDEBUG){std::cout<<"Setting new basis to "<<Lmax_new<<" ... ";}
	//this->_s->set_new_basis(Lmax_new, this->_s->b->dim);
	this->_s->set_L_max(Lmax_new);
	if(KKDEBUG){std::cout<<"Done ! "<<std::endl;}
	//this->on_initializeButton_clicked();
	cscale = 1.0/cscale_calc();
	_axes_actor->SetTotalLength(_scale_fac/cscale, _scale_fac/cscale, _scale_fac/cscale);
	this->update_sf(this->comboBox_sf->currentIndex());		// forces a recalculation of the scalar field
	synchronize_shape2sliders();
	update_vtkwindow();
}
void shape::update_comboBox_sf()
{
	if(_s->sf_tags.size())
	{
		this->comboBox_sf->blockSignals(true);
		this->comboBox_sf->clear();
		this->comboBox_sf->addItem("none");
		this->comboBox_sf->addItem("theta");
		this->comboBox_sf->addItem("phi");
		this->comboBox_sf->addItem("mean curvature (analytical)");
		this->comboBox_sf->addItem("mean curvature (mesh approx.)");
		this->comboBox_sf->addItem("Gaussian curvature (analytical)");
		for(int i=0;i<_s->sf_tags.size();i++)
		{
			this->comboBox_sf->addItem(_s->sf_tags[i].c_str());
		}
		this->comboBox_sf->blockSignals(false);
	}
}
double shape::cscale_calc()
{
	// find maximum value among  the L = 1 coefficients
	double mc = 0.0; //mc = max.coeff.
	for(int i=1;i<4;i++){if(std::abs(_s->xc(i))>mc){mc = std::abs(_s->xc(i));}}
	for(int i=1;i<4;i++){if(std::abs(_s->yc(i))>mc){mc = std::abs(_s->yc(i));}}
	for(int i=1;i<4;i++){if(std::abs(_s->zc(i))>mc){mc = std::abs(_s->zc(i));}}
	return abs(mc);	// dividing by the largest value
}
void shape::on_comboBox_gdim_currentIndexChanged(int val)
{
	int new_dim = 40;		// corresponds to case 5
	if(val==0){new_dim = 5;}
	if(val==1){new_dim = 10;}
	if(val==2){new_dim = 15;}
	if(val==3){new_dim = 20;}
	if(val==4){new_dim = 30;}
	if(val==5){new_dim = 40;}
	if(val==6){new_dim = 50;}
	if(val==7){new_dim = 60;}
	if(val==8){new_dim = 80;}
	if(val==9){new_dim = 120;}
	if(val==10){new_dim = 240;}
	if(val==11){new_dim = 320;}
	if(val>11){new_dim = 20;if(KKDEBUG){std::cout<<"Warning: Setting to default!"<<std::endl;}}

	if(KKDEBUG){std::cout<<"Changing basis to L_max = "<<_s->b->L_max<<" and dim :"<<new_dim<<" ..... ";}
	_s->set_new_basis(_s->b->L_max, new_dim);
	if(KKDEBUG){std::cout<<"Done !"<<std::endl;}
	_s->needsUpdating = YES;
	update_vtkwindow();

}
void shape::update_sf(int val)
{
	if(val<=0){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = 1.0;}}
	if(val==1){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->t(i);}	}
	if(val==2){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->p(i);}	}

	if(val==3){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->Han[i];}	}
	if(val==4){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->H[i];}	}
	if(val==5){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->KGan[i];}	}

	if(val==6){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->X[i][0];}}
	if(val==7){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->X[i][1];}}	
	if(val==8){	for(int i=0;i<_s->sm->n_points;i++)	{_s->sm->sf[i] = _s->sm->X[i][2];}}	
	if(val>8){this->_s->sfGen(val-9);}
}
void shape::on_comboBox_sf_currentIndexChanged(int val)
{
	this->update_vtkwindow();
}
void shape::on_checkBox_self_intersection_stateChanged(int val)
{
	this->check_self_intersection = bool(val);
	this->update_vtkwindow();
}
void shape::on_comboBox_spherical_mesh_currentIndexChanged(int val)
{
	int new_n = 3;		// corresponds to default
	if(val==0){new_n = 1;}
	if(val==1){new_n = 2;}
	if(val==2){new_n = 3;}
	if(val==3){new_n = 4;}
	if(val==4){new_n = 5;}
	if(val==5){new_n = 6;}

	if(KKDEBUG){std::cout<<"Changing mesh representation..... ";}
	_s->set_new_spherical_mesh(_s->b->L_max, new_n);
	_polys = vtkCellArray::New();
	for (int i=0; i<_s->sm->n_faces; i++)
	{
		_polys->InsertNextCell(3);
		_polys->InsertCellPoint(_s->sm->f0[i]-1);
		_polys->InsertCellPoint(_s->sm->f1[i]-1);
		_polys->InsertCellPoint(_s->sm->f2[i]-1);
	}
	if(KKDEBUG){std::cout<<"Done !"<<std::endl;}
	
	if(new_n>2)
	{
		if(KKDEBUG){std::cout<<".... self-intersection test disabled! "<<std::endl;}
		checkBox_self_intersection->setChecked(false);
		this->checkBox_self_intersection->setEnabled(false);
		this->check_self_intersection= NO;
	}
	else
	{
		this->checkBox_self_intersection->setEnabled(true);
		if(this->checkBox_self_intersection->checkState() ){this->check_self_intersection=YES;}
	}

	_s->needsUpdating = YES;
	update_vtkwindow();
}





// Slider value change handling --------- z ------------
void shape::synchronize_shape2sliders(){
/// block all signals, so that set value doesn't trigger a valueChanged
verticalSliderx01->blockSignals(true);
verticalSliderx02->blockSignals(true);
verticalSliderx03->blockSignals(true);
verticalSliderx04->blockSignals(true);
verticalSliderx05->blockSignals(true);
verticalSliderx06->blockSignals(true);
verticalSliderx07->blockSignals(true);
verticalSliderx08->blockSignals(true);
verticalSliderx09->blockSignals(true);
verticalSliderx10->blockSignals(true);
verticalSliderx11->blockSignals(true);
verticalSliderx12->blockSignals(true);
verticalSliderx13->blockSignals(true);
verticalSliderx14->blockSignals(true);
verticalSliderx15->blockSignals(true);
verticalSliderx16->blockSignals(true);
verticalSliderx17->blockSignals(true);
verticalSliderx18->blockSignals(true);
verticalSliderx19->blockSignals(true);
verticalSliderx20->blockSignals(true);
verticalSliderx21->blockSignals(true);
verticalSliderx22->blockSignals(true);
verticalSliderx23->blockSignals(true);
verticalSliderx24->blockSignals(true);
verticalSliderx25->blockSignals(true);
verticalSliderx26->blockSignals(true);
verticalSliderx27->blockSignals(true);
verticalSliderx28->blockSignals(true);
verticalSliderx29->blockSignals(true);
verticalSliderx30->blockSignals(true);
verticalSliderx31->blockSignals(true);
verticalSliderx32->blockSignals(true);
verticalSliderx33->blockSignals(true);
verticalSliderx34->blockSignals(true);
verticalSliderx35->blockSignals(true);
verticalSliderx36->blockSignals(true);
verticalSliderx37->blockSignals(true);
verticalSliderx38->blockSignals(true);
verticalSliderx39->blockSignals(true);
verticalSliderx40->blockSignals(true);
verticalSliderx41->blockSignals(true);
verticalSliderx42->blockSignals(true);
verticalSliderx43->blockSignals(true);
verticalSliderx44->blockSignals(true);
verticalSliderx45->blockSignals(true);
verticalSliderx46->blockSignals(true);
verticalSliderx47->blockSignals(true);
verticalSliderx48->blockSignals(true);
verticalSliderx49->blockSignals(true);


verticalSliderx01xx->blockSignals(true);
verticalSliderx02xx->blockSignals(true);
verticalSliderx03xx->blockSignals(true);
verticalSliderx04xx->blockSignals(true);
verticalSliderx05xx->blockSignals(true);
verticalSliderx06xx->blockSignals(true);
verticalSliderx07xx->blockSignals(true);
verticalSliderx08xx->blockSignals(true);
verticalSliderx09xx->blockSignals(true);
verticalSliderx10xx->blockSignals(true);
verticalSliderx11xx->blockSignals(true);
verticalSliderx12xx->blockSignals(true);
verticalSliderx13xx->blockSignals(true);
verticalSliderx14xx->blockSignals(true);
verticalSliderx15xx->blockSignals(true);
verticalSliderx16xx->blockSignals(true);
verticalSliderx17xx->blockSignals(true);
verticalSliderx18xx->blockSignals(true);
verticalSliderx19xx->blockSignals(true);
verticalSliderx20xx->blockSignals(true);
verticalSliderx21xx->blockSignals(true);
verticalSliderx22xx->blockSignals(true);
verticalSliderx23xx->blockSignals(true);
verticalSliderx24xx->blockSignals(true);
verticalSliderx25xx->blockSignals(true);
verticalSliderx26xx->blockSignals(true);
verticalSliderx27xx->blockSignals(true);
verticalSliderx28xx->blockSignals(true);
verticalSliderx29xx->blockSignals(true);
verticalSliderx30xx->blockSignals(true);
verticalSliderx31xx->blockSignals(true);
verticalSliderx32xx->blockSignals(true);
verticalSliderx33xx->blockSignals(true);
verticalSliderx34xx->blockSignals(true);
verticalSliderx35xx->blockSignals(true);
verticalSliderx36xx->blockSignals(true);
verticalSliderx37xx->blockSignals(true);
verticalSliderx38xx->blockSignals(true);
verticalSliderx39xx->blockSignals(true);
verticalSliderx40xx->blockSignals(true);
verticalSliderx41xx->blockSignals(true);
verticalSliderx42xx->blockSignals(true);
verticalSliderx43xx->blockSignals(true);
verticalSliderx44xx->blockSignals(true);
verticalSliderx45xx->blockSignals(true);
verticalSliderx46xx->blockSignals(true);
verticalSliderx47xx->blockSignals(true);
verticalSliderx48xx->blockSignals(true);
verticalSliderx49xx->blockSignals(true);

verticalSliderx01yy->blockSignals(true);
verticalSliderx02yy->blockSignals(true);
verticalSliderx03yy->blockSignals(true);
verticalSliderx04yy->blockSignals(true);
verticalSliderx05yy->blockSignals(true);
verticalSliderx06yy->blockSignals(true);
verticalSliderx07yy->blockSignals(true);
verticalSliderx08yy->blockSignals(true);
verticalSliderx09yy->blockSignals(true);
verticalSliderx10yy->blockSignals(true);
verticalSliderx11yy->blockSignals(true);
verticalSliderx12yy->blockSignals(true);
verticalSliderx13yy->blockSignals(true);
verticalSliderx14yy->blockSignals(true);
verticalSliderx15yy->blockSignals(true);
verticalSliderx16yy->blockSignals(true);
verticalSliderx17yy->blockSignals(true);
verticalSliderx18yy->blockSignals(true);
verticalSliderx19yy->blockSignals(true);
verticalSliderx20yy->blockSignals(true);
verticalSliderx21yy->blockSignals(true);
verticalSliderx22yy->blockSignals(true);
verticalSliderx23yy->blockSignals(true);
verticalSliderx24yy->blockSignals(true);
verticalSliderx25yy->blockSignals(true);
verticalSliderx26yy->blockSignals(true);
verticalSliderx27yy->blockSignals(true);
verticalSliderx28yy->blockSignals(true);
verticalSliderx29yy->blockSignals(true);
verticalSliderx30yy->blockSignals(true);
verticalSliderx31yy->blockSignals(true);
verticalSliderx32yy->blockSignals(true);
verticalSliderx33yy->blockSignals(true);
verticalSliderx34yy->blockSignals(true);
verticalSliderx35yy->blockSignals(true);
verticalSliderx36yy->blockSignals(true);
verticalSliderx37yy->blockSignals(true);
verticalSliderx38yy->blockSignals(true);
verticalSliderx39yy->blockSignals(true);
verticalSliderx40yy->blockSignals(true);
verticalSliderx41yy->blockSignals(true);
verticalSliderx42yy->blockSignals(true);
verticalSliderx43yy->blockSignals(true);
verticalSliderx44yy->blockSignals(true);
verticalSliderx45yy->blockSignals(true);
verticalSliderx46yy->blockSignals(true);
verticalSliderx47yy->blockSignals(true);
verticalSliderx48yy->blockSignals(true);
verticalSliderx49yy->blockSignals(true);

 verticalSliderx01->setValue(( _s->zc(0) *49*cscale)+49);
 verticalSliderx02->setValue(( _s->zc(1) *49*cscale)+49);
 verticalSliderx03->setValue(( _s->zc(2) *49*cscale)+49);
 verticalSliderx04->setValue(( _s->zc(3) *49*cscale)+49);
 verticalSliderx05->setValue(( _s->zc(4) *49*cscale)+49);
 verticalSliderx06->setValue(( _s->zc(5) *49*cscale)+49);
 verticalSliderx07->setValue(( _s->zc(6) *49*cscale)+49);
 verticalSliderx08->setValue(( _s->zc(7) *49*cscale)+49);
 verticalSliderx09->setValue(( _s->zc(8) *49*cscale)+49);
 verticalSliderx10->setValue(( _s->zc(9) *49*cscale)+49);
 verticalSliderx11->setValue(( _s->zc(10) *49*cscale)+49);
 verticalSliderx12->setValue(( _s->zc(11) *49*cscale)+49);
 verticalSliderx13->setValue(( _s->zc(12) *49*cscale)+49);
 verticalSliderx14->setValue(( _s->zc(13) *49*cscale)+49);
 verticalSliderx15->setValue(( _s->zc(14) *49*cscale)+49);
 verticalSliderx16->setValue(( _s->zc(15) *49*cscale)+49);
 verticalSliderx17->setValue(( _s->zc(16) *49*cscale)+49);
 verticalSliderx18->setValue(( _s->zc(17) *49*cscale)+49);
 verticalSliderx19->setValue(( _s->zc(18) *49*cscale)+49);
 verticalSliderx20->setValue(( _s->zc(19) *49*cscale)+49);
 verticalSliderx21->setValue(( _s->zc(20) *49*cscale)+49);
 verticalSliderx22->setValue(( _s->zc(21) *49*cscale)+49);
 verticalSliderx23->setValue(( _s->zc(22) *49*cscale)+49);
 verticalSliderx24->setValue(( _s->zc(23) *49*cscale)+49);
 verticalSliderx25->setValue(( _s->zc(24) *49*cscale)+49);
 verticalSliderx26->setValue(( _s->zc(25) *49*cscale)+49);
 verticalSliderx27->setValue(( _s->zc(26) *49*cscale)+49);
 verticalSliderx28->setValue(( _s->zc(27) *49*cscale)+49);
 verticalSliderx29->setValue(( _s->zc(28) *49*cscale)+49);
 verticalSliderx30->setValue(( _s->zc(29) *49*cscale)+49);
 verticalSliderx31->setValue(( _s->zc(30) *49*cscale)+49);
 verticalSliderx32->setValue(( _s->zc(31) *49*cscale)+49);
 verticalSliderx33->setValue(( _s->zc(32) *49*cscale)+49);
 verticalSliderx34->setValue(( _s->zc(33) *49*cscale)+49);
 verticalSliderx35->setValue(( _s->zc(34) *49*cscale)+49);
 verticalSliderx36->setValue(( _s->zc(35) *49*cscale)+49);
 verticalSliderx37->setValue(( _s->zc(36) *49*cscale)+49);
 verticalSliderx38->setValue(( _s->zc(37) *49*cscale)+49);
 verticalSliderx39->setValue(( _s->zc(38) *49*cscale)+49);
 verticalSliderx40->setValue(( _s->zc(39) *49*cscale)+49);
 verticalSliderx41->setValue(( _s->zc(40) *49*cscale)+49);
 verticalSliderx42->setValue(( _s->zc(41) *49*cscale)+49);
 verticalSliderx43->setValue(( _s->zc(42) *49*cscale)+49);
 verticalSliderx44->setValue(( _s->zc(43) *49*cscale)+49);
 verticalSliderx45->setValue(( _s->zc(44) *49*cscale)+49);
 verticalSliderx46->setValue(( _s->zc(45) *49*cscale)+49);
 verticalSliderx47->setValue(( _s->zc(46) *49*cscale)+49);
 verticalSliderx48->setValue(( _s->zc(47) *49*cscale)+49);
 verticalSliderx49->setValue(( _s->zc(48) *49*cscale)+49);



// Slider value change handling --------- x ------------

 verticalSliderx01xx->setValue(( _s->xc(0) *49*cscale)+49);
 verticalSliderx02xx->setValue(( _s->xc(1) *49*cscale)+49);
 verticalSliderx03xx->setValue(( _s->xc(2) *49*cscale)+49);
 verticalSliderx04xx->setValue(( _s->xc(3) *49*cscale)+49);
 verticalSliderx05xx->setValue(( _s->xc(4) *49*cscale)+49);
 verticalSliderx06xx->setValue(( _s->xc(5) *49*cscale)+49);
 verticalSliderx07xx->setValue(( _s->xc(6) *49*cscale)+49);
 verticalSliderx08xx->setValue(( _s->xc(7) *49*cscale)+49);
 verticalSliderx09xx->setValue(( _s->xc(8) *49*cscale)+49);
 verticalSliderx10xx->setValue(( _s->xc(9) *49*cscale)+49);
 verticalSliderx11xx->setValue(( _s->xc(10) *49*cscale)+49);
 verticalSliderx12xx->setValue(( _s->xc(11) *49*cscale)+49);
 verticalSliderx13xx->setValue(( _s->xc(12) *49*cscale)+49);
 verticalSliderx14xx->setValue(( _s->xc(13) *49*cscale)+49);
 verticalSliderx15xx->setValue(( _s->xc(14) *49*cscale)+49);
 verticalSliderx16xx->setValue(( _s->xc(15) *49*cscale)+49);
 verticalSliderx17xx->setValue(( _s->xc(16) *49*cscale)+49);
 verticalSliderx18xx->setValue(( _s->xc(17) *49*cscale)+49);
 verticalSliderx19xx->setValue(( _s->xc(18) *49*cscale)+49);
 verticalSliderx20xx->setValue(( _s->xc(19) *49*cscale)+49);
 verticalSliderx21xx->setValue(( _s->xc(20) *49*cscale)+49);
 verticalSliderx22xx->setValue(( _s->xc(21) *49*cscale)+49);
 verticalSliderx23xx->setValue(( _s->xc(22) *49*cscale)+49);
 verticalSliderx24xx->setValue(( _s->xc(23) *49*cscale)+49);
 verticalSliderx25xx->setValue(( _s->xc(24) *49*cscale)+49);
 verticalSliderx26xx->setValue(( _s->xc(25) *49*cscale)+49);
 verticalSliderx27xx->setValue(( _s->xc(26) *49*cscale)+49);
 verticalSliderx28xx->setValue(( _s->xc(27) *49*cscale)+49);
 verticalSliderx29xx->setValue(( _s->xc(28) *49*cscale)+49);
 verticalSliderx30xx->setValue(( _s->xc(29) *49*cscale)+49);
 verticalSliderx31xx->setValue(( _s->xc(30) *49*cscale)+49);
 verticalSliderx32xx->setValue(( _s->xc(31) *49*cscale)+49);
 verticalSliderx33xx->setValue(( _s->xc(32) *49*cscale)+49);
 verticalSliderx34xx->setValue(( _s->xc(33) *49*cscale)+49);
 verticalSliderx35xx->setValue(( _s->xc(34) *49*cscale)+49);
 verticalSliderx36xx->setValue(( _s->xc(35) *49*cscale)+49);
 verticalSliderx37xx->setValue(( _s->xc(36) *49*cscale)+49);
 verticalSliderx38xx->setValue(( _s->xc(37) *49*cscale)+49);
 verticalSliderx39xx->setValue(( _s->xc(38) *49*cscale)+49);
 verticalSliderx40xx->setValue(( _s->xc(39) *49*cscale)+49);
 verticalSliderx41xx->setValue(( _s->xc(40) *49*cscale)+49);
 verticalSliderx42xx->setValue(( _s->xc(41) *49*cscale)+49);
 verticalSliderx43xx->setValue(( _s->xc(42) *49*cscale)+49);
 verticalSliderx44xx->setValue(( _s->xc(43) *49*cscale)+49);
 verticalSliderx45xx->setValue(( _s->xc(44) *49*cscale)+49);
 verticalSliderx46xx->setValue(( _s->xc(45) *49*cscale)+49);
 verticalSliderx47xx->setValue(( _s->xc(46) *49*cscale)+49);
 verticalSliderx48xx->setValue(( _s->xc(47) *49*cscale)+49);
 verticalSliderx49xx->setValue(( _s->xc(48) *49*cscale)+49);

// Slider value change handling --------- y ------------

 verticalSliderx01yy->setValue(( _s->yc(0) *49*cscale)+49);
 verticalSliderx02yy->setValue(( _s->yc(1) *49*cscale)+49);
 verticalSliderx03yy->setValue(( _s->yc(2) *49*cscale)+49);
 verticalSliderx04yy->setValue(( _s->yc(3) *49*cscale)+49);
 verticalSliderx05yy->setValue(( _s->yc(4) *49*cscale)+49);
 verticalSliderx06yy->setValue(( _s->yc(5) *49*cscale)+49);
 verticalSliderx07yy->setValue(( _s->yc(6) *49*cscale)+49);
 verticalSliderx08yy->setValue(( _s->yc(7) *49*cscale)+49);
 verticalSliderx09yy->setValue(( _s->yc(8) *49*cscale)+49);
 verticalSliderx10yy->setValue(( _s->yc(9) *49*cscale)+49);
 verticalSliderx11yy->setValue(( _s->yc(10) *49*cscale)+49);
 verticalSliderx12yy->setValue(( _s->yc(11) *49*cscale)+49);
 verticalSliderx13yy->setValue(( _s->yc(12) *49*cscale)+49);
 verticalSliderx14yy->setValue(( _s->yc(13) *49*cscale)+49);
 verticalSliderx15yy->setValue(( _s->yc(14) *49*cscale)+49);
 verticalSliderx16yy->setValue(( _s->yc(15) *49*cscale)+49);
 verticalSliderx17yy->setValue(( _s->yc(16) *49*cscale)+49);
 verticalSliderx18yy->setValue(( _s->yc(17) *49*cscale)+49);
 verticalSliderx19yy->setValue(( _s->yc(18) *49*cscale)+49);
 verticalSliderx20yy->setValue(( _s->yc(19) *49*cscale)+49);
 verticalSliderx21yy->setValue(( _s->yc(20) *49*cscale)+49);
 verticalSliderx22yy->setValue(( _s->yc(21) *49*cscale)+49);
 verticalSliderx23yy->setValue(( _s->yc(22) *49*cscale)+49);
 verticalSliderx24yy->setValue(( _s->yc(23) *49*cscale)+49);
 verticalSliderx25yy->setValue(( _s->yc(24) *49*cscale)+49);
 verticalSliderx26yy->setValue(( _s->yc(25) *49*cscale)+49);
 verticalSliderx27yy->setValue(( _s->yc(26) *49*cscale)+49);
 verticalSliderx28yy->setValue(( _s->yc(27) *49*cscale)+49);
 verticalSliderx29yy->setValue(( _s->yc(28) *49*cscale)+49);
 verticalSliderx30yy->setValue(( _s->yc(29) *49*cscale)+49);
 verticalSliderx31yy->setValue(( _s->yc(30) *49*cscale)+49);
 verticalSliderx32yy->setValue(( _s->yc(31) *49*cscale)+49);
 verticalSliderx33yy->setValue(( _s->yc(32) *49*cscale)+49);
 verticalSliderx34yy->setValue(( _s->yc(33) *49*cscale)+49);
 verticalSliderx35yy->setValue(( _s->yc(34) *49*cscale)+49);
 verticalSliderx36yy->setValue(( _s->yc(35) *49*cscale)+49);
 verticalSliderx37yy->setValue(( _s->yc(36) *49*cscale)+49);
 verticalSliderx38yy->setValue(( _s->yc(37) *49*cscale)+49);
 verticalSliderx39yy->setValue(( _s->yc(38) *49*cscale)+49);
 verticalSliderx40yy->setValue(( _s->yc(39) *49*cscale)+49);
 verticalSliderx41yy->setValue(( _s->yc(40) *49*cscale)+49);
 verticalSliderx42yy->setValue(( _s->yc(41) *49*cscale)+49);
 verticalSliderx43yy->setValue(( _s->yc(42) *49*cscale)+49);
 verticalSliderx44yy->setValue(( _s->yc(43) *49*cscale)+49);
 verticalSliderx45yy->setValue(( _s->yc(44) *49*cscale)+49);
 verticalSliderx46yy->setValue(( _s->yc(45) *49*cscale)+49);
 verticalSliderx47yy->setValue(( _s->yc(46) *49*cscale)+49);
 verticalSliderx48yy->setValue(( _s->yc(47) *49*cscale)+49);
 verticalSliderx49yy->setValue(( _s->yc(48) *49*cscale)+49);


/// release block all signals
verticalSliderx01->blockSignals(false);
verticalSliderx02->blockSignals(false);
verticalSliderx03->blockSignals(false);
verticalSliderx04->blockSignals(false);
verticalSliderx05->blockSignals(false);
verticalSliderx06->blockSignals(false);
verticalSliderx07->blockSignals(false);
verticalSliderx08->blockSignals(false);
verticalSliderx09->blockSignals(false);
verticalSliderx10->blockSignals(false);
verticalSliderx11->blockSignals(false);
verticalSliderx12->blockSignals(false);
verticalSliderx13->blockSignals(false);
verticalSliderx14->blockSignals(false);
verticalSliderx15->blockSignals(false);
verticalSliderx16->blockSignals(false);
verticalSliderx17->blockSignals(false);
verticalSliderx18->blockSignals(false);
verticalSliderx19->blockSignals(false);
verticalSliderx20->blockSignals(false);
verticalSliderx21->blockSignals(false);
verticalSliderx22->blockSignals(false);
verticalSliderx23->blockSignals(false);
verticalSliderx24->blockSignals(false);
verticalSliderx25->blockSignals(false);
verticalSliderx26->blockSignals(false);
verticalSliderx27->blockSignals(false);
verticalSliderx28->blockSignals(false);
verticalSliderx29->blockSignals(false);
verticalSliderx30->blockSignals(false);
verticalSliderx31->blockSignals(false);
verticalSliderx32->blockSignals(false);
verticalSliderx33->blockSignals(false);
verticalSliderx34->blockSignals(false);
verticalSliderx35->blockSignals(false);
verticalSliderx36->blockSignals(false);
verticalSliderx37->blockSignals(false);
verticalSliderx38->blockSignals(false);
verticalSliderx39->blockSignals(false);
verticalSliderx40->blockSignals(false);
verticalSliderx41->blockSignals(false);
verticalSliderx42->blockSignals(false);
verticalSliderx43->blockSignals(false);
verticalSliderx44->blockSignals(false);
verticalSliderx45->blockSignals(false);
verticalSliderx46->blockSignals(false);
verticalSliderx47->blockSignals(false);
verticalSliderx48->blockSignals(false);
verticalSliderx49->blockSignals(false);


verticalSliderx01xx->blockSignals(false);
verticalSliderx02xx->blockSignals(false);
verticalSliderx03xx->blockSignals(false);
verticalSliderx04xx->blockSignals(false);
verticalSliderx05xx->blockSignals(false);
verticalSliderx06xx->blockSignals(false);
verticalSliderx07xx->blockSignals(false);
verticalSliderx08xx->blockSignals(false);
verticalSliderx09xx->blockSignals(false);
verticalSliderx10xx->blockSignals(false);
verticalSliderx11xx->blockSignals(false);
verticalSliderx12xx->blockSignals(false);
verticalSliderx13xx->blockSignals(false);
verticalSliderx14xx->blockSignals(false);
verticalSliderx15xx->blockSignals(false);
verticalSliderx16xx->blockSignals(false);
verticalSliderx17xx->blockSignals(false);
verticalSliderx18xx->blockSignals(false);
verticalSliderx19xx->blockSignals(false);
verticalSliderx20xx->blockSignals(false);
verticalSliderx21xx->blockSignals(false);
verticalSliderx22xx->blockSignals(false);
verticalSliderx23xx->blockSignals(false);
verticalSliderx24xx->blockSignals(false);
verticalSliderx25xx->blockSignals(false);
verticalSliderx26xx->blockSignals(false);
verticalSliderx27xx->blockSignals(false);
verticalSliderx28xx->blockSignals(false);
verticalSliderx29xx->blockSignals(false);
verticalSliderx30xx->blockSignals(false);
verticalSliderx31xx->blockSignals(false);
verticalSliderx32xx->blockSignals(false);
verticalSliderx33xx->blockSignals(false);
verticalSliderx34xx->blockSignals(false);
verticalSliderx35xx->blockSignals(false);
verticalSliderx36xx->blockSignals(false);
verticalSliderx37xx->blockSignals(false);
verticalSliderx38xx->blockSignals(false);
verticalSliderx39xx->blockSignals(false);
verticalSliderx40xx->blockSignals(false);
verticalSliderx41xx->blockSignals(false);
verticalSliderx42xx->blockSignals(false);
verticalSliderx43xx->blockSignals(false);
verticalSliderx44xx->blockSignals(false);
verticalSliderx45xx->blockSignals(false);
verticalSliderx46xx->blockSignals(false);
verticalSliderx47xx->blockSignals(false);
verticalSliderx48xx->blockSignals(false);
verticalSliderx49xx->blockSignals(false);

verticalSliderx01yy->blockSignals(false);
verticalSliderx02yy->blockSignals(false);
verticalSliderx03yy->blockSignals(false);
verticalSliderx04yy->blockSignals(false);
verticalSliderx05yy->blockSignals(false);
verticalSliderx06yy->blockSignals(false);
verticalSliderx07yy->blockSignals(false);
verticalSliderx08yy->blockSignals(false);
verticalSliderx09yy->blockSignals(false);
verticalSliderx10yy->blockSignals(false);
verticalSliderx11yy->blockSignals(false);
verticalSliderx12yy->blockSignals(false);
verticalSliderx13yy->blockSignals(false);
verticalSliderx14yy->blockSignals(false);
verticalSliderx15yy->blockSignals(false);
verticalSliderx16yy->blockSignals(false);
verticalSliderx17yy->blockSignals(false);
verticalSliderx18yy->blockSignals(false);
verticalSliderx19yy->blockSignals(false);
verticalSliderx20yy->blockSignals(false);
verticalSliderx21yy->blockSignals(false);
verticalSliderx22yy->blockSignals(false);
verticalSliderx23yy->blockSignals(false);
verticalSliderx24yy->blockSignals(false);
verticalSliderx25yy->blockSignals(false);
verticalSliderx26yy->blockSignals(false);
verticalSliderx27yy->blockSignals(false);
verticalSliderx28yy->blockSignals(false);
verticalSliderx29yy->blockSignals(false);
verticalSliderx30yy->blockSignals(false);
verticalSliderx31yy->blockSignals(false);
verticalSliderx32yy->blockSignals(false);
verticalSliderx33yy->blockSignals(false);
verticalSliderx34yy->blockSignals(false);
verticalSliderx35yy->blockSignals(false);
verticalSliderx36yy->blockSignals(false);
verticalSliderx37yy->blockSignals(false);
verticalSliderx38yy->blockSignals(false);
verticalSliderx39yy->blockSignals(false);
verticalSliderx40yy->blockSignals(false);
verticalSliderx41yy->blockSignals(false);
verticalSliderx42yy->blockSignals(false);
verticalSliderx43yy->blockSignals(false);
verticalSliderx44yy->blockSignals(false);
verticalSliderx45yy->blockSignals(false);
verticalSliderx46yy->blockSignals(false);
verticalSliderx47yy->blockSignals(false);
verticalSliderx48yy->blockSignals(false);
verticalSliderx49yy->blockSignals(false);
}




void shape::on_verticalSliderx01_valueChanged(int val){int coordix = 2; int cix = 0; double old_val = _s->zc(0);_s->zc(0) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx02_valueChanged(int val){int coordix = 2; int cix = 1; double old_val = _s->zc(1);_s->zc(1) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx03_valueChanged(int val){int coordix = 2; int cix = 2; double old_val = _s->zc(2);_s->zc(2) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx04_valueChanged(int val){int coordix = 2; int cix = 3; double old_val = _s->zc(3);_s->zc(3) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx05_valueChanged(int val){int coordix = 2; int cix = 4; double old_val = _s->zc(4);_s->zc(4) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx06_valueChanged(int val){int coordix = 2; int cix = 5; double old_val = _s->zc(5);_s->zc(5) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx07_valueChanged(int val){int coordix = 2; int cix = 6; double old_val = _s->zc(6);_s->zc(6) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx08_valueChanged(int val){int coordix = 2; int cix = 7; double old_val = _s->zc(7);_s->zc(7) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx09_valueChanged(int val){int coordix = 2; int cix = 8; double old_val = _s->zc(8);_s->zc(8) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx10_valueChanged(int val){int coordix = 2; int cix = 9; double old_val = _s->zc(9);_s->zc(9) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx11_valueChanged(int val){int coordix = 2; int cix = 10; double old_val = _s->zc(10);_s->zc(10) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx12_valueChanged(int val){int coordix = 2; int cix = 11; double old_val = _s->zc(11);_s->zc(11) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx13_valueChanged(int val){int coordix = 2; int cix = 12; double old_val = _s->zc(12);_s->zc(12) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx14_valueChanged(int val){int coordix = 2; int cix = 13; double old_val = _s->zc(13);_s->zc(13) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx15_valueChanged(int val){int coordix = 2; int cix = 14; double old_val = _s->zc(14);_s->zc(14) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx16_valueChanged(int val){int coordix = 2; int cix = 15; double old_val = _s->zc(15);_s->zc(15) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx17_valueChanged(int val){int coordix = 2; int cix = 16; double old_val = _s->zc(16);_s->zc(16) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx18_valueChanged(int val){int coordix = 2; int cix = 17; double old_val = _s->zc(17);_s->zc(17) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx19_valueChanged(int val){int coordix = 2; int cix = 18; double old_val = _s->zc(18);_s->zc(18) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx20_valueChanged(int val){int coordix = 2; int cix = 19; double old_val = _s->zc(19);_s->zc(19) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx21_valueChanged(int val){int coordix = 2; int cix = 20; double old_val = _s->zc(20);_s->zc(20) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx22_valueChanged(int val){int coordix = 2; int cix = 21; double old_val = _s->zc(21);_s->zc(21) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx23_valueChanged(int val){int coordix = 2; int cix = 22; double old_val = _s->zc(22);_s->zc(22) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx24_valueChanged(int val){int coordix = 2; int cix = 23; double old_val = _s->zc(23);_s->zc(23) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx25_valueChanged(int val){int coordix = 2; int cix = 24; double old_val = _s->zc(24);_s->zc(24) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx26_valueChanged(int val){int coordix = 2; int cix = 25; double old_val = _s->zc(25);_s->zc(25) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx27_valueChanged(int val){int coordix = 2; int cix = 26; double old_val = _s->zc(26);_s->zc(26) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx28_valueChanged(int val){int coordix = 2; int cix = 27; double old_val = _s->zc(27);_s->zc(27) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx29_valueChanged(int val){int coordix = 2; int cix = 28; double old_val = _s->zc(28);_s->zc(28) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx30_valueChanged(int val){int coordix = 2; int cix = 29; double old_val = _s->zc(29);_s->zc(29) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx31_valueChanged(int val){int coordix = 2; int cix = 30; double old_val = _s->zc(30);_s->zc(30) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx32_valueChanged(int val){int coordix = 2; int cix = 31; double old_val = _s->zc(31);_s->zc(31) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx33_valueChanged(int val){int coordix = 2; int cix = 32; double old_val = _s->zc(32);_s->zc(32) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx34_valueChanged(int val){int coordix = 2; int cix = 33; double old_val = _s->zc(33);_s->zc(33) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx35_valueChanged(int val){int coordix = 2; int cix = 34; double old_val = _s->zc(34);_s->zc(34) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx36_valueChanged(int val){int coordix = 2; int cix = 35; double old_val = _s->zc(35);_s->zc(35) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx37_valueChanged(int val){int coordix = 2; int cix = 36; double old_val = _s->zc(36);_s->zc(36) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx38_valueChanged(int val){int coordix = 2; int cix = 37; double old_val = _s->zc(37);_s->zc(37) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx39_valueChanged(int val){int coordix = 2; int cix = 38; double old_val = _s->zc(38);_s->zc(38) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx40_valueChanged(int val){int coordix = 2; int cix = 39; double old_val = _s->zc(39);_s->zc(39) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx41_valueChanged(int val){int coordix = 2; int cix = 40; double old_val = _s->zc(40);_s->zc(40) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx42_valueChanged(int val){int coordix = 2; int cix = 41; double old_val = _s->zc(41);_s->zc(41) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx43_valueChanged(int val){int coordix = 2; int cix = 42; double old_val = _s->zc(42);_s->zc(42) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx44_valueChanged(int val){int coordix = 2; int cix = 43; double old_val = _s->zc(43);_s->zc(43) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx45_valueChanged(int val){int coordix = 2; int cix = 44; double old_val = _s->zc(44);_s->zc(44) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx46_valueChanged(int val){int coordix = 2; int cix = 45; double old_val = _s->zc(45);_s->zc(45) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx47_valueChanged(int val){int coordix = 2; int cix = 46; double old_val = _s->zc(46);_s->zc(46) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx48_valueChanged(int val){int coordix = 2; int cix = 47; double old_val = _s->zc(47);_s->zc(47) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx49_valueChanged(int val){int coordix = 2; int cix = 48; double old_val = _s->zc(48);_s->zc(48) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
///

// Slider value change handling --------- x ------------

void shape::on_verticalSliderx01xx_valueChanged(int val){int coordix = 0; int cix = 0; double old_val = _s->xc(0);_s->xc(0) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx02xx_valueChanged(int val){int coordix = 0; int cix = 1; double old_val = _s->xc(1);_s->xc(1) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx03xx_valueChanged(int val){int coordix = 0; int cix = 2; double old_val = _s->xc(2);_s->xc(2) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx04xx_valueChanged(int val){int coordix = 0; int cix = 3; double old_val = _s->xc(3);_s->xc(3) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx05xx_valueChanged(int val){int coordix = 0; int cix = 4; double old_val = _s->xc(4);_s->xc(4) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx06xx_valueChanged(int val){int coordix = 0; int cix = 5; double old_val = _s->xc(5);_s->xc(5) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx07xx_valueChanged(int val){int coordix = 0; int cix = 6; double old_val = _s->xc(6);_s->xc(6) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx08xx_valueChanged(int val){int coordix = 0; int cix = 7; double old_val = _s->xc(7);_s->xc(7) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx09xx_valueChanged(int val){int coordix = 0; int cix = 8; double old_val = _s->xc(8);_s->xc(8) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx10xx_valueChanged(int val){int coordix = 0; int cix = 9; double old_val = _s->xc(9);_s->xc(9) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx11xx_valueChanged(int val){int coordix = 0; int cix = 10; double old_val = _s->xc(10);_s->xc(10) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx12xx_valueChanged(int val){int coordix = 0; int cix = 11; double old_val = _s->xc(11);_s->xc(11) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx13xx_valueChanged(int val){int coordix = 0; int cix = 12; double old_val = _s->xc(12);_s->xc(12) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx14xx_valueChanged(int val){int coordix = 0; int cix = 13; double old_val = _s->xc(13);_s->xc(13) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx15xx_valueChanged(int val){int coordix = 0; int cix = 14; double old_val = _s->xc(14);_s->xc(14) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx16xx_valueChanged(int val){int coordix = 0; int cix = 15; double old_val = _s->xc(15);_s->xc(15) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx17xx_valueChanged(int val){int coordix = 0; int cix = 16; double old_val = _s->xc(16);_s->xc(16) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx18xx_valueChanged(int val){int coordix = 0; int cix = 17; double old_val = _s->xc(17);_s->xc(17) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx19xx_valueChanged(int val){int coordix = 0; int cix = 18; double old_val = _s->xc(18);_s->xc(18) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx20xx_valueChanged(int val){int coordix = 0; int cix = 19; double old_val = _s->xc(19);_s->xc(19) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx21xx_valueChanged(int val){int coordix = 0; int cix = 20; double old_val = _s->xc(20);_s->xc(20) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx22xx_valueChanged(int val){int coordix = 0; int cix = 21; double old_val = _s->xc(21);_s->xc(21) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx23xx_valueChanged(int val){int coordix = 0; int cix = 22; double old_val = _s->xc(22);_s->xc(22) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx24xx_valueChanged(int val){int coordix = 0; int cix = 23; double old_val = _s->xc(23);_s->xc(23) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx25xx_valueChanged(int val){int coordix = 0; int cix = 24; double old_val = _s->xc(24);_s->xc(24) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx26xx_valueChanged(int val){int coordix = 0; int cix = 25; double old_val = _s->xc(25);_s->xc(25) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx27xx_valueChanged(int val){int coordix = 0; int cix = 26; double old_val = _s->xc(26);_s->xc(26) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx28xx_valueChanged(int val){int coordix = 0; int cix = 27; double old_val = _s->xc(27);_s->xc(27) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx29xx_valueChanged(int val){int coordix = 0; int cix = 28; double old_val = _s->xc(28);_s->xc(28) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx30xx_valueChanged(int val){int coordix = 0; int cix = 29; double old_val = _s->xc(29);_s->xc(29) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx31xx_valueChanged(int val){int coordix = 0; int cix = 30; double old_val = _s->xc(30);_s->xc(30) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx32xx_valueChanged(int val){int coordix = 0; int cix = 31; double old_val = _s->xc(31);_s->xc(31) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx33xx_valueChanged(int val){int coordix = 0; int cix = 32; double old_val = _s->xc(32);_s->xc(32) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx34xx_valueChanged(int val){int coordix = 0; int cix = 33; double old_val = _s->xc(33);_s->xc(33) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx35xx_valueChanged(int val){int coordix = 0; int cix = 34; double old_val = _s->xc(34);_s->xc(34) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx36xx_valueChanged(int val){int coordix = 0; int cix = 35; double old_val = _s->xc(35);_s->xc(35) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx37xx_valueChanged(int val){int coordix = 0; int cix = 36; double old_val = _s->xc(36);_s->xc(36) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx38xx_valueChanged(int val){int coordix = 0; int cix = 37; double old_val = _s->xc(37);_s->xc(37) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx39xx_valueChanged(int val){int coordix = 0; int cix = 38; double old_val = _s->xc(38);_s->xc(38) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx40xx_valueChanged(int val){int coordix = 0; int cix = 39; double old_val = _s->xc(39);_s->xc(39) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx41xx_valueChanged(int val){int coordix = 0; int cix = 40; double old_val = _s->xc(40);_s->xc(40) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx42xx_valueChanged(int val){int coordix = 0; int cix = 41; double old_val = _s->xc(41);_s->xc(41) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx43xx_valueChanged(int val){int coordix = 0; int cix = 42; double old_val = _s->xc(42);_s->xc(42) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx44xx_valueChanged(int val){int coordix = 0; int cix = 43; double old_val = _s->xc(43);_s->xc(43) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx45xx_valueChanged(int val){int coordix = 0; int cix = 44; double old_val = _s->xc(44);_s->xc(44) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx46xx_valueChanged(int val){int coordix = 0; int cix = 45; double old_val = _s->xc(45);_s->xc(45) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx47xx_valueChanged(int val){int coordix = 0; int cix = 46; double old_val = _s->xc(46);_s->xc(46) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx48xx_valueChanged(int val){int coordix = 0; int cix = 47; double old_val = _s->xc(47);_s->xc(47) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx49xx_valueChanged(int val){int coordix = 0; int cix = 48; double old_val = _s->xc(48);_s->xc(48) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}

// Slider value change handling --------- y ------------

void shape::on_verticalSliderx01yy_valueChanged(int val){int coordix = 1; int cix = 0; double old_val = _s->yc(0);_s->yc(0) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx02yy_valueChanged(int val){int coordix = 1; int cix = 1; double old_val = _s->yc(1);_s->yc(1) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx03yy_valueChanged(int val){int coordix = 1; int cix = 2; double old_val = _s->yc(2);_s->yc(2) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx04yy_valueChanged(int val){int coordix = 1; int cix = 3; double old_val = _s->yc(3);_s->yc(3) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx05yy_valueChanged(int val){int coordix = 1; int cix = 4; double old_val = _s->yc(4);_s->yc(4) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx06yy_valueChanged(int val){int coordix = 1; int cix = 5; double old_val = _s->yc(5);_s->yc(5) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx07yy_valueChanged(int val){int coordix = 1; int cix = 6; double old_val = _s->yc(6);_s->yc(6) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx08yy_valueChanged(int val){int coordix = 1; int cix = 7; double old_val = _s->yc(7);_s->yc(7) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx09yy_valueChanged(int val){int coordix = 1; int cix = 8; double old_val = _s->yc(8);_s->yc(8) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx10yy_valueChanged(int val){int coordix = 1; int cix = 9; double old_val = _s->yc(9);_s->yc(9) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx11yy_valueChanged(int val){int coordix = 1; int cix = 10; double old_val = _s->yc(10);_s->yc(10) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx12yy_valueChanged(int val){int coordix = 1; int cix = 11; double old_val = _s->yc(11);_s->yc(11) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx13yy_valueChanged(int val){int coordix = 1; int cix = 12; double old_val = _s->yc(12);_s->yc(12) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx14yy_valueChanged(int val){int coordix = 1; int cix = 13; double old_val = _s->yc(13);_s->yc(13) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx15yy_valueChanged(int val){int coordix = 1; int cix = 14; double old_val = _s->yc(14);_s->yc(14) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx16yy_valueChanged(int val){int coordix = 1; int cix = 15; double old_val = _s->yc(15);_s->yc(15) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx17yy_valueChanged(int val){int coordix = 1; int cix = 16; double old_val = _s->yc(16);_s->yc(16) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx18yy_valueChanged(int val){int coordix = 1; int cix = 17; double old_val = _s->yc(17);_s->yc(17) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx19yy_valueChanged(int val){int coordix = 1; int cix = 18; double old_val = _s->yc(18);_s->yc(18) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx20yy_valueChanged(int val){int coordix = 1; int cix = 19; double old_val = _s->yc(19);_s->yc(19) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx21yy_valueChanged(int val){int coordix = 1; int cix = 20; double old_val = _s->yc(20);_s->yc(20) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx22yy_valueChanged(int val){int coordix = 1; int cix = 21; double old_val = _s->yc(21);_s->yc(21) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx23yy_valueChanged(int val){int coordix = 1; int cix = 22; double old_val = _s->yc(22);_s->yc(22) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx24yy_valueChanged(int val){int coordix = 1; int cix = 23; double old_val = _s->yc(23);_s->yc(23) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx25yy_valueChanged(int val){int coordix = 1; int cix = 24; double old_val = _s->yc(24);_s->yc(24) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx26yy_valueChanged(int val){int coordix = 1; int cix = 25; double old_val = _s->yc(25);_s->yc(25) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx27yy_valueChanged(int val){int coordix = 1; int cix = 26; double old_val = _s->yc(26);_s->yc(26) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx28yy_valueChanged(int val){int coordix = 1; int cix = 27; double old_val = _s->yc(27);_s->yc(27) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx29yy_valueChanged(int val){int coordix = 1; int cix = 28; double old_val = _s->yc(28);_s->yc(28) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx30yy_valueChanged(int val){int coordix = 1; int cix = 29; double old_val = _s->yc(29);_s->yc(29) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx31yy_valueChanged(int val){int coordix = 1; int cix = 30; double old_val = _s->yc(30);_s->yc(30) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx32yy_valueChanged(int val){int coordix = 1; int cix = 31; double old_val = _s->yc(31);_s->yc(31) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx33yy_valueChanged(int val){int coordix = 1; int cix = 32; double old_val = _s->yc(32);_s->yc(32) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx34yy_valueChanged(int val){int coordix = 1; int cix = 33; double old_val = _s->yc(33);_s->yc(33) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx35yy_valueChanged(int val){int coordix = 1; int cix = 34; double old_val = _s->yc(34);_s->yc(34) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx36yy_valueChanged(int val){int coordix = 1; int cix = 35; double old_val = _s->yc(35);_s->yc(35) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx37yy_valueChanged(int val){int coordix = 1; int cix = 36; double old_val = _s->yc(36);_s->yc(36) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx38yy_valueChanged(int val){int coordix = 1; int cix = 37; double old_val = _s->yc(37);_s->yc(37) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx39yy_valueChanged(int val){int coordix = 1; int cix = 38; double old_val = _s->yc(38);_s->yc(38) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx40yy_valueChanged(int val){int coordix = 1; int cix = 39; double old_val = _s->yc(39);_s->yc(39) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx41yy_valueChanged(int val){int coordix = 1; int cix = 40; double old_val = _s->yc(40);_s->yc(40) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx42yy_valueChanged(int val){int coordix = 1; int cix = 41; double old_val = _s->yc(41);_s->yc(41) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx43yy_valueChanged(int val){int coordix = 1; int cix = 42; double old_val = _s->yc(42);_s->yc(42) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx44yy_valueChanged(int val){int coordix = 1; int cix = 43; double old_val = _s->yc(43);_s->yc(43) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx45yy_valueChanged(int val){int coordix = 1; int cix = 44; double old_val = _s->yc(44);_s->yc(44) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx46yy_valueChanged(int val){int coordix = 1; int cix = 45; double old_val = _s->yc(45);_s->yc(45) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx47yy_valueChanged(int val){int coordix = 1; int cix = 46; double old_val = _s->yc(46);_s->yc(46) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx48yy_valueChanged(int val){int coordix = 1; int cix = 47; double old_val = _s->yc(47);_s->yc(47) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
void shape::on_verticalSliderx49yy_valueChanged(int val){int coordix = 1; int cix = 48; double old_val = _s->yc(48);_s->yc(48) = ((double)(val-49))/49/cscale ;_s->needsUpdating = YES;update_vtkwindow(coordix, cix, old_val,((double)(val-49))/49/cscale);}
