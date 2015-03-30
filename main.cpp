#include "shape.h"
#include <QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	shape w;
	w.show();
	return a.exec();
}
