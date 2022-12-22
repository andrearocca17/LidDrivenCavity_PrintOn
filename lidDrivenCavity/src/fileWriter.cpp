#include "fileWriter.h"
#include <iostream>
#include <string>
fileWriter::fileWriter()
{
}

fileWriter::~fileWriter()
{
}

void fileWriter::writeUVP(string& name,int time_, Grid& Grid_,Fields::vectorField& Utemp, Fields::vectorField& Vtemp,  Fields::vectorField& Ptemp)
  {
    string  myfileType = ".dat";
    int ttime = time_;
    string timename ;
    ostringstream temp;
    temp << ttime;
    timename = temp.str();
    string newnames_ = name;
    string newfilename= newnames_.append(timename);
    string name2 = newfilename.append(myfileType);
    int size = name.size() + 10;
    char* cstr = new char();
    strcpy_s(cstr, size, newfilename.c_str());
    FILE* outfile=nullptr;
    if (outfile == NULL) perror("Error opening file");  
    else {
        fopen_s(&outfile, cstr, "w+t");
        int NXtemp = Utemp.size();
        int NYtemp = Utemp[0].size();
        fprintf(outfile, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
        fprintf(outfile, "ZONE  F=POINT\n");
        fprintf(outfile, "I=%d, J=%d\n", NXtemp, NYtemp);
        double xpos, ypos, UU, VV, PP;
        for (int i = 0; i < Utemp.size(); i++)
        {
            for (int j = 0; j < Utemp[0].size(); j++)
            {
                xpos = Grid_.XC[i];
                ypos = Grid_.YC[j];
                UU = Utemp[i][j].value;
                VV = Vtemp[i][j].value;
                PP = Ptemp[i][j].value;

                //%5.8lf\t
                fprintf(outfile, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, UU, VV, PP);
            }
        }

        fclose(outfile);
    }
  }

