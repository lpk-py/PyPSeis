#include <QCoreApplication>
#include <stdio.h>
#include "calculator.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    Calculator Calculator;


    Calculator.ReadData("E:\\Microseismic\\MSET\\MyWorks\\MSET_DEMO\\DEMO_FORWARD4\\demo4_49g_500us.sgy");

    Calculator.SetPseudoWell();
    Calculator.GetGrid();

    QString geoFile = "E:\\Microseismic\\MSET\\MyWorks\\MSET_DEMO\\DEMO_FORWARD4\\demo4_geo-coor.prn";
    Calculator.ReadGeophones(geoFile);

    QString velFile = "E:\\Microseismic\\MSET\\MyWorks\\MSET_DEMO\\DEMO_FORWARD4\\velocity_SS103H_true_new.txt";
    Calculator.LoadVModel(velFile);

    //Calculator.CalTTT();
    Calculator.CalSET();

    return a.exec();
}
