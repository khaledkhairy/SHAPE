#ifndef SHAPE_H
#define SHAPE_H

#include <QMainWindow>
#include "ui_shape.h"


// Forward class declarations
class vtkDataSetMapper;
class vtkPolyDataMapper;
class vtkPolyData;
class vtkActor;
class vtkRenderer;
class vtkCellArray;
class shp_surface;
class vtkAxesActor;
class vtkCubeAxesActor;
class vtkScalarBarActor;
class vtkLookupTable;

class shape : public QMainWindow, private Ui_MainWindow
{
	Q_OBJECT

public:
	//shape(QWidget *parent = 0, Qt::WFlags flags = 0);
	shape(QWidget *parent = 0);
	~shape();

	private:
	bool cube_axis_actor1_on;
	bool scalar_bar_vis;
 	double _minH;
	double _maxH;
    vtkRenderer* _ren1;
	vtkCellArray*_polys;
    vtkAxesActor* _axes_actor;
	vtkPolyDataMapper* Mapper;
	vtkActor* Actor;
	vtkCubeAxesActor *_cube_axes_actor1;
	vtkScalarBarActor* scalarBar;
	vtkLookupTable *lut;
	shp_surface *_s;
	bool check_self_intersection;
	bool _axis_on;
	double cscale, _scale_fac, _bgrdr, _bgrdg, _bgrdb;		// scale factor for clks necessary when importing shapes and default background color
	void shape2polydata(const shp_surface *s, vtkPolyData*);
	void update_vtkwindow();
	void update_vtkwindow(int coordix, int cix, double old_val, double new_val);
	void synchronize_shape2sliders();
	double cscale_calc();
	QString mimeData_2_fileName( const QMimeData * mimeData );
	void dragEnterEvent( QDragEnterEvent * event );
	void dropEvent( QDropEvent * event );
	void update_comboBox_sf();
	void update_sf(int);
	void import_shape_from_disc(std::string filestr);
	void _update_gui();
//private:
	//Ui::MainWindow ui;

	public slots:
	virtual void on_lineEdit_Lmax_returnPressed();
	virtual	void on_initializeButton_clicked();
	virtual void on_axes_toggle_Button_clicked();
	virtual void on_pushButton_cube_axis_toggle_clicked();
	virtual void on_actionOpen_triggered();
	virtual void on_actionSave_triggered(); 
	virtual void on_actionExport_triggered();
	virtual void on_actionExport_to_STL_triggered();
	virtual void on_actionExport_to_PLY_triggered();
	virtual	void on_actionExit_triggered();
	virtual void on_comboBox_gdim_currentIndexChanged(int);
	virtual void on_comboBox_spherical_mesh_currentIndexChanged(int);
	virtual void on_comboBox_sf_currentIndexChanged(int);
	virtual void on_checkBox_self_intersection_stateChanged(int);
	 
/// Slider declarations ---- z -----

	void on_verticalSliderx01_valueChanged(int);
   	void on_verticalSliderx02_valueChanged(int);
	void on_verticalSliderx03_valueChanged(int);
   	void on_verticalSliderx04_valueChanged(int);
	void on_verticalSliderx05_valueChanged(int);
   	void on_verticalSliderx06_valueChanged(int);
	void on_verticalSliderx07_valueChanged(int);
   	void on_verticalSliderx08_valueChanged(int);
	void on_verticalSliderx09_valueChanged(int);
   	void on_verticalSliderx10_valueChanged(int);
	void on_verticalSliderx11_valueChanged(int);
   	void on_verticalSliderx12_valueChanged(int);
	void on_verticalSliderx13_valueChanged(int);
   	void on_verticalSliderx14_valueChanged(int);
   	void on_verticalSliderx15_valueChanged(int);
	void on_verticalSliderx16_valueChanged(int);
   	void on_verticalSliderx17_valueChanged(int);
	void on_verticalSliderx18_valueChanged(int);
   	void on_verticalSliderx19_valueChanged(int);
   	void on_verticalSliderx20_valueChanged(int);
	void on_verticalSliderx21_valueChanged(int);
   	void on_verticalSliderx22_valueChanged(int);
	void on_verticalSliderx23_valueChanged(int);
   	void on_verticalSliderx24_valueChanged(int);
	void on_verticalSliderx25_valueChanged(int);
   	void on_verticalSliderx26_valueChanged(int);
   	void on_verticalSliderx27_valueChanged(int);
   	void on_verticalSliderx28_valueChanged(int);
	void on_verticalSliderx29_valueChanged(int);
   	void on_verticalSliderx30_valueChanged(int);
	void on_verticalSliderx31_valueChanged(int);
   	void on_verticalSliderx32_valueChanged(int);
	void on_verticalSliderx33_valueChanged(int);
   	void on_verticalSliderx34_valueChanged(int);
   	void on_verticalSliderx35_valueChanged(int);
   	void on_verticalSliderx36_valueChanged(int);
	void on_verticalSliderx37_valueChanged(int);
   	void on_verticalSliderx38_valueChanged(int);
	void on_verticalSliderx39_valueChanged(int);
   	void on_verticalSliderx40_valueChanged(int);
	void on_verticalSliderx41_valueChanged(int);
   	void on_verticalSliderx42_valueChanged(int);
   	void on_verticalSliderx43_valueChanged(int);
   	void on_verticalSliderx44_valueChanged(int);
	void on_verticalSliderx45_valueChanged(int);
   	void on_verticalSliderx46_valueChanged(int);
   	void on_verticalSliderx47_valueChanged(int);
   	void on_verticalSliderx48_valueChanged(int);
   	void on_verticalSliderx49_valueChanged(int);

/// Slider declarations ---- x -----

	void on_verticalSliderx01xx_valueChanged(int);
   	void on_verticalSliderx02xx_valueChanged(int);
	void on_verticalSliderx03xx_valueChanged(int);
   	void on_verticalSliderx04xx_valueChanged(int);
	void on_verticalSliderx05xx_valueChanged(int);
   	void on_verticalSliderx06xx_valueChanged(int);
	void on_verticalSliderx07xx_valueChanged(int);
   	void on_verticalSliderx08xx_valueChanged(int);
	void on_verticalSliderx09xx_valueChanged(int);
   	void on_verticalSliderx10xx_valueChanged(int);
	void on_verticalSliderx11xx_valueChanged(int);
   	void on_verticalSliderx12xx_valueChanged(int);
	void on_verticalSliderx13xx_valueChanged(int);
   	void on_verticalSliderx14xx_valueChanged(int);
   	void on_verticalSliderx15xx_valueChanged(int);
	void on_verticalSliderx16xx_valueChanged(int);
   	void on_verticalSliderx17xx_valueChanged(int);
	void on_verticalSliderx18xx_valueChanged(int);
   	void on_verticalSliderx19xx_valueChanged(int);
   	void on_verticalSliderx20xx_valueChanged(int);
	void on_verticalSliderx21xx_valueChanged(int);
   	void on_verticalSliderx22xx_valueChanged(int);
	void on_verticalSliderx23xx_valueChanged(int);
   	void on_verticalSliderx24xx_valueChanged(int);
	void on_verticalSliderx25xx_valueChanged(int);
   	void on_verticalSliderx26xx_valueChanged(int);
   	void on_verticalSliderx27xx_valueChanged(int);
   	void on_verticalSliderx28xx_valueChanged(int);
	void on_verticalSliderx29xx_valueChanged(int);
   	void on_verticalSliderx30xx_valueChanged(int);
	void on_verticalSliderx31xx_valueChanged(int);
   	void on_verticalSliderx32xx_valueChanged(int);
	void on_verticalSliderx33xx_valueChanged(int);
   	void on_verticalSliderx34xx_valueChanged(int);
   	void on_verticalSliderx35xx_valueChanged(int);
   	void on_verticalSliderx36xx_valueChanged(int);
	void on_verticalSliderx37xx_valueChanged(int);
   	void on_verticalSliderx38xx_valueChanged(int);
	void on_verticalSliderx39xx_valueChanged(int);
   	void on_verticalSliderx40xx_valueChanged(int);
	void on_verticalSliderx41xx_valueChanged(int);
   	void on_verticalSliderx42xx_valueChanged(int);
   	void on_verticalSliderx43xx_valueChanged(int);
   	void on_verticalSliderx44xx_valueChanged(int);
	void on_verticalSliderx45xx_valueChanged(int);
   	void on_verticalSliderx46xx_valueChanged(int);
   	void on_verticalSliderx47xx_valueChanged(int);
   	void on_verticalSliderx48xx_valueChanged(int);
   	void on_verticalSliderx49xx_valueChanged(int);

/// Slider declarations ---- y -----

	void on_verticalSliderx01yy_valueChanged(int);
   	void on_verticalSliderx02yy_valueChanged(int);
	void on_verticalSliderx03yy_valueChanged(int);
   	void on_verticalSliderx04yy_valueChanged(int);
	void on_verticalSliderx05yy_valueChanged(int);
   	void on_verticalSliderx06yy_valueChanged(int);
	void on_verticalSliderx07yy_valueChanged(int);
   	void on_verticalSliderx08yy_valueChanged(int);
	void on_verticalSliderx09yy_valueChanged(int);
   	void on_verticalSliderx10yy_valueChanged(int);
	void on_verticalSliderx11yy_valueChanged(int);
   	void on_verticalSliderx12yy_valueChanged(int);
	void on_verticalSliderx13yy_valueChanged(int);
   	void on_verticalSliderx14yy_valueChanged(int);
   	void on_verticalSliderx15yy_valueChanged(int);
	void on_verticalSliderx16yy_valueChanged(int);
   	void on_verticalSliderx17yy_valueChanged(int);
	void on_verticalSliderx18yy_valueChanged(int);
   	void on_verticalSliderx19yy_valueChanged(int);
   	void on_verticalSliderx20yy_valueChanged(int);
	void on_verticalSliderx21yy_valueChanged(int);
   	void on_verticalSliderx22yy_valueChanged(int);
	void on_verticalSliderx23yy_valueChanged(int);
   	void on_verticalSliderx24yy_valueChanged(int);
	void on_verticalSliderx25yy_valueChanged(int);
   	void on_verticalSliderx26yy_valueChanged(int);
   	void on_verticalSliderx27yy_valueChanged(int);
   	void on_verticalSliderx28yy_valueChanged(int);
	void on_verticalSliderx29yy_valueChanged(int);
   	void on_verticalSliderx30yy_valueChanged(int);
	void on_verticalSliderx31yy_valueChanged(int);
   	void on_verticalSliderx32yy_valueChanged(int);
	void on_verticalSliderx33yy_valueChanged(int);
   	void on_verticalSliderx34yy_valueChanged(int);
   	void on_verticalSliderx35yy_valueChanged(int);
   	void on_verticalSliderx36yy_valueChanged(int);
	void on_verticalSliderx37yy_valueChanged(int);
   	void on_verticalSliderx38yy_valueChanged(int);
	void on_verticalSliderx39yy_valueChanged(int);
   	void on_verticalSliderx40yy_valueChanged(int);
	void on_verticalSliderx41yy_valueChanged(int);
   	void on_verticalSliderx42yy_valueChanged(int);
   	void on_verticalSliderx43yy_valueChanged(int);
   	void on_verticalSliderx44yy_valueChanged(int);
	void on_verticalSliderx45yy_valueChanged(int);
   	void on_verticalSliderx46yy_valueChanged(int);
   	void on_verticalSliderx47yy_valueChanged(int);
   	void on_verticalSliderx48yy_valueChanged(int);
   	void on_verticalSliderx49yy_valueChanged(int);




private slots:


};

#endif // SHAPE_H
